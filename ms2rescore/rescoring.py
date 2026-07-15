"""
Rescore PSMs with ristretto.

:py:mod:`ristretto` is a lean, dependency-light (numpy/scikit-learn/pandas only)
reimplementation of the Percolator/Käll 2007 semi-supervised rescoring algorithm, purpose-built
to serve MS²Rescore. It takes a plain :py:class:`pandas.DataFrame` in and returns a
:py:class:`~ristretto.result.RescoreResult` of plain DataFrames out.

"""

import logging
import re
from concurrent.futures import BrokenExecutor
from dataclasses import replace
from pathlib import Path
from typing import Dict, Optional, Set, Tuple

import numpy as np
import pandas as pd
import ristretto
from psm_utils import PSMList
from ristretto import RescoreResult

from ms2rescore.constants import CHARGE_PATTERN
from ms2rescore.exceptions import RescoringError
from ms2rescore.parse_psms import infer_score_direction

logger = logging.getLogger(__name__)


def build_features_dataframe(
    psm_list: PSMList,
    feature_names: Set[str],
    lower_score_is_better: bool,
) -> pd.DataFrame:
    """
    Build the plain identifier + feature DataFrame that ristretto expects.

    Identifier columns are the real ``spectrum_id`` (for true spectrum-grouped CV, not an
    artificial per-row index), a charge-stripped ``peptidoform`` string, a bare (no
    modifications, no charge) ``peptide`` string, and a semicolon-joined ``protein`` string
    (only added if every PSM has a non-empty ``protein_list``). ``score`` is negated if
    ``lower_score_is_better``, so ristretto always sees "higher is better".

    """
    psm_df = psm_list.to_dataframe().reset_index(drop=True)

    features_df = pd.DataFrame(
        {
            "spectrum_id": psm_df["spectrum_id"],
            "is_decoy": psm_df["is_decoy"],
            "peptidoform": psm_df["peptidoform"]
            .astype(str)
            .str.replace(CHARGE_PATTERN, "", n=1, regex=True),
            "peptide": psm_df["peptidoform"].apply(lambda p: p.sequence),
            "score": -psm_df["score"] if lower_score_is_better else psm_df["score"],
        }
    )

    if psm_df["protein_list"].apply(lambda p: bool(p)).all():
        features_df["protein"] = psm_df["protein_list"].apply(";".join)

    if feature_names:
        feature_df = pd.DataFrame(list(psm_df["rescoring_features"]))[sorted(feature_names)]
        feature_df = feature_df.astype(float).fillna(0.0)
        features_df = pd.concat([features_df, feature_df], axis=1)

    return features_df


def evaluate_before(psm_list: PSMList, config: Dict) -> RescoreResult:
    """
    Evaluate the PSMs' current (pre-rescoring) score with ristretto, for report baselines.

    Called right after :py:func:`~ms2rescore.parse_psms.parse_psms`, before any feature
    generator runs. Always competes to one best PSM per spectrum, regardless of
    ``max_psm_rank_output``: the report is always a rank-1, one-row-per-spectrum view, so
    this is directly comparable to :py:func:`rescore`'s ``report_result``.

    """
    train_fdr = config["rescoring"]["train_fdr"] if config["rescoring"] else 0.01
    features_df = build_features_dataframe(psm_list, set(), infer_score_direction(psm_list, train_fdr))
    return ristretto.evaluate(
        features_df,
        peptide_col="peptide",
        protein_col="protein" if "protein" in features_df.columns else None,
        decoy_pattern=config["id_decoy_pattern"],
        multi_rank_rescoring=False,
    )


def rescore(
    psm_list: PSMList, config: Dict, output_file_root: str
) -> Tuple[PSMList, RescoreResult, RescoreResult]:
    """
    Rescore PSMs with ristretto and write the new scores, q-values, and PEPs back to ``psm_list``.

    Always scores every candidate rank with ristretto's semi-supervised model
    (``multi_rank_rescoring=True``), then runs the actual competition/FDR/PEP step
    separately via :py:func:`ristretto.evaluate`, competing down to one PSM per spectrum
    unless ``max_psm_rank_output > 1`` (mirroring the existing rank-based config options).

    Returns
    -------
    tuple[PSMList, RescoreResult, RescoreResult]
        The (possibly filtered) PSMList with updated scores; the final competition/FDR
        result respecting ``max_psm_rank_output`` (matches ``psm_list`` and the written
        output files); and a report-view result always competed to one PSM per spectrum,
        for the HTML report and identification counts, regardless of
        ``max_psm_rank_output``. The two results are the same object when
        ``max_psm_rank_output == 1``.

    """
    feature_names = {f for psm in psm_list for f in psm.rescoring_features}
    lower_score_is_better = infer_score_direction(psm_list, config["rescoring"]["train_fdr"])
    features_df = build_features_dataframe(psm_list, feature_names, lower_score_is_better)

    peptide_col = "peptide"
    protein_col = "protein" if "protein" in features_df.columns else None
    decoy_pattern = config["id_decoy_pattern"]

    try:
        ml_result = ristretto.rescore(
            features_df,
            peptide_col=peptide_col,
            protein_col=protein_col,
            feature_cols=sorted(feature_names),
            decoy_pattern=decoy_pattern,
            train_fdr=config["rescoring"]["train_fdr"],
            n_jobs=int(config["processes"]),
            multi_rank_rescoring=True,
        )
        # For max_psm_rank_output > 1, trim to the top-N ranks (by ML score) *before* the
        # final FDR calculation, so q-values/PEP are computed over exactly the population
        # that ends up in the output -- not recomputed after the fact on a subset.
        ml_psms = ml_result.psms
        if config["max_psm_rank_output"] > 1:
            output_rank = ml_psms.groupby("spectrum_id")["score"].rank(
                method="first", ascending=False
            )
            output_ml_psms = ml_psms[output_rank <= config["max_psm_rank_output"]]
        else:
            output_ml_psms = ml_psms

        final_result = ristretto.evaluate(
            output_ml_psms,
            peptide_col=peptide_col,
            protein_col=protein_col,
            decoy_pattern=decoy_pattern,
            multi_rank_rescoring=config["max_psm_rank_output"] != 1,
        )
        report_result = None
        if config["max_psm_rank_output"] > 1:
            # Separate, always-rank-1 competition on the full (pre-trim) ML scores, purely
            # for the report -- independent of how many ranks per spectrum the user wants
            # in the actual output.
            report_result = ristretto.evaluate(
                ml_psms,
                peptide_col=peptide_col,
                protein_col=protein_col,
                decoy_pattern=decoy_pattern,
                multi_rank_rescoring=False,
            )
    except (RuntimeError, BrokenExecutor, ValueError) as e:
        raise RescoringError("Ristretto could not be run. Please check the input data.") from e

    keep_mask = np.zeros(len(psm_list), dtype=bool)
    keep_mask[final_result.psms.index.to_numpy()] = True
    psm_list = psm_list[keep_mask]

    psm_list["score"] = final_result.psms["score"].to_numpy()
    psm_list["qvalue"] = final_result.psms["qvalue"].to_numpy()
    psm_list["pep"] = final_result.psms["pep"].to_numpy()
    psm_list.set_ranks(lower_score_better=False)

    psm_list, final_result = fix_constant_pep(psm_list, final_result)
    if report_result is None:
        # max_psm_rank_output == 1: the report view is the same competition as the final
        # result, so it must reflect the same fix_constant_pep filtering.
        report_result = final_result

    _write_group_metadata(psm_list, final_result.psms, final_result.peptidoforms, "peptidoform")
    if final_result.peptides is not None:
        _write_group_metadata(psm_list, final_result.psms, final_result.peptides, "peptide")
    if final_result.proteins is not None:
        if decoy_pattern:
            logger.debug("Protein-level FDR used picked-protein competition (decoy_pattern set).")
        else:
            logger.warning(
                "id_decoy_pattern is not set: protein-level FDR falls back to a plain rollup "
                "instead of picked-protein competition. Protein q-values/PEPs may be less "
                "accurate. Set id_decoy_pattern to enable picked-protein competition."
            )
        _write_group_metadata(
            psm_list, final_result.psms, final_result.proteins, "protein", decoy_pattern
        )

    return psm_list, final_result, report_result


def _write_group_metadata(
    psm_list: PSMList,
    psms: pd.DataFrame,
    rollup: pd.DataFrame,
    group_col: str,
    decoy_pattern: Optional[str] = None,
) -> None:
    """
    Write a rollup's score/qvalue/pep onto each PSM's metadata, keyed by ``group_col``.

    For protein-level rollups computed with ``decoy_pattern`` set, ristretto's
    ``picked_rollup`` strips the decoy pattern from each group key (picked-protein
    competition, Savitski et al. 2015), so it must be stripped from ``psms[group_col]``
    the same way before the lookup.

    """
    lookup = rollup.set_index(group_col)[["score", "qvalue", "pep"]].to_dict(orient="index")
    strip_pattern = re.compile(decoy_pattern) if group_col == "protein" and decoy_pattern else None
    for psm, group_key in zip(psm_list, psms[group_col]):
        if strip_pattern is not None:
            group_key = strip_pattern.sub("", group_key)
        values = lookup[group_key]
        psm.metadata.update(
            {
                f"{group_col}_score": values["score"],
                f"{group_col}_qvalue": values["qvalue"],
                f"{group_col}_pep": values["pep"],
            }
        )


def write_rescoring_tables(
    before: RescoreResult, after: RescoreResult, output_file_root: str
) -> None:
    """Write before/after rescoring result tables to Parquet files."""
    _write_result_tables(before, output_file_root, "before")
    _write_result_tables(after, output_file_root, "after")
    after.feature_weights.to_parquet(f"{output_file_root}.ristretto.weights.parquet")


def _write_result_tables(result: RescoreResult, output_file_root: str, suffix: str) -> None:
    """Write one RescoreResult's psms/peptidoforms/peptides/proteins tables to Parquet."""
    result.psms.to_parquet(Path(f"{output_file_root}.ristretto.psms_{suffix}.parquet"))
    result.peptidoforms.to_parquet(
        Path(f"{output_file_root}.ristretto.peptidoforms_{suffix}.parquet")
    )
    if result.peptides is not None:
        result.peptides.to_parquet(Path(f"{output_file_root}.ristretto.peptides_{suffix}.parquet"))
    if result.proteins is not None:
        result.proteins.to_parquet(Path(f"{output_file_root}.ristretto.proteins_{suffix}.parquet"))


def fix_constant_pep(
    psm_list: PSMList, result: RescoreResult
) -> Tuple[PSMList, RescoreResult]:
    """
    Workaround for broken PEP calculation if the best-scoring PSM is a decoy.

    Filters ``psm_list`` and ``result.psms`` together (same length, same row order) so they
    never desync -- callers (parquet tables, report, 1%-FDR counts) all keep seeing the same
    PSM population as the returned ``psm_list``.

    """
    if not all(psm_list["pep"] == 1.0):
        return psm_list, result

    logger.warning(
        "Attempting to fix constant PEP values by removing decoy PSMs that score higher than the "
        "best target PSM."
    )
    max_target_score = psm_list["score"][~psm_list["is_decoy"]].max()
    higher_scoring_decoys = psm_list["is_decoy"] & (psm_list["score"] > max_target_score)

    if not higher_scoring_decoys.any():
        logger.warning("No decoys scoring higher than the best target found. Skipping fix.")
        return psm_list, result

    keep = ~higher_scoring_decoys
    logger.warning(f"Removed {higher_scoring_decoys.sum()} decoy PSMs.")
    return psm_list[keep], replace(result, psms=result.psms[keep])
