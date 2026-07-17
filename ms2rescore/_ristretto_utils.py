"""
Package-internal helpers around ristretto, not part of the public API.

The public entry point is :py:func:`ms2rescore.rescoring.rescore`. Everything here backs
that function (and is imported by it), plus a few standalone evaluation/reporting helpers
used elsewhere in the package (:py:mod:`ms2rescore.core`, :py:mod:`ms2rescore.report.data`) --
none of it is meant to be called directly by external users.

"""

import logging
import re
from dataclasses import replace
from typing import Dict, Optional, Set, Tuple

import numpy as np
import pandas as pd
import ristretto
from psm_utils import PSMList
from ristretto import RescoreResult

from ms2rescore.constants import CHARGE_PATTERN
from ms2rescore.parse_psms import infer_score_direction

logger = logging.getLogger(__name__)


def _build_features_dataframe(
    psm_list: PSMList,
    feature_names: Set[str],
    lower_score_is_better: bool,
) -> pd.DataFrame:
    """
    Build the plain identifier + feature DataFrame that ristretto expects.

    Identifier columns are the real ``spectrum_id`` and ``run`` (for true spectrum-grouped
    CV and competition, not an artificial per-row index -- ``run`` guards against
    ``spectrum_id`` collisions across multiple input files), a charge-stripped
    ``peptidoform`` string, a bare (no modifications, no charge) ``peptide`` string, and a
    semicolon-joined ``protein`` string (only added if every PSM has a non-empty
    ``protein_list``). ``score`` is negated if ``lower_score_is_better``, so ristretto
    always sees "higher is better".

    """
    psm_df = psm_list.to_dataframe().reset_index(drop=True)

    features_df = pd.DataFrame(
        {
            "spectrum_id": psm_df["spectrum_id"],
            "run": psm_df["run"],
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


def _trim_and_evaluate(
    features_df: pd.DataFrame,
    max_rank: int,
    *,
    run_col: str,
    peptide_col: Optional[str],
    protein_col: Optional[str],
    decoy_pattern: Optional[str],
) -> RescoreResult:
    """
    Compete to at most ``max_rank`` PSMs per spectrum, then compute q-values/PEP/rollups.

    ``max_rank == 1``: proper target-decoy spectrum competition (ristretto's
    ``multi_rank_rescoring=False``), matching Percolator/mokapot/Sage default behavior.

    ``max_rank > 1``: keep the top ``max_rank`` PSMs per ``(run, spectrum_id)`` by score
    first (a plain rank cut, not a competition -- ristretto has no "top-N" concept, only
    "best-1" or "keep-everything"), then compute q-values/PEP over that whole population as
    independent rows (``multi_rank_rescoring=True``). This is *not* statistically rigorous
    FDR control -- ``max_psm_rank_output > 1`` is meant for surfacing ambiguous results
    (e.g. multiple candidate peptidoforms per spectrum from mumble), not for a corrected
    identification count.

    """
    if max_rank == 1:
        trimmed = features_df
        multi_rank = False
    else:
        rank = features_df.groupby([run_col, "spectrum_id"])["score"].rank(
            method="first", ascending=False
        )
        trimmed = features_df[rank <= max_rank]
        multi_rank = True

    return ristretto.evaluate(
        trimmed,
        run_col=run_col,
        peptide_col=peptide_col,
        protein_col=protein_col,
        decoy_pattern=decoy_pattern,
        multi_rank_rescoring=multi_rank,
    )


def _is_original_psm(psm) -> bool:
    """
    Whether a PSM is original (not mumble-generated), robust to metadata's round-trip.

    ``psm_utils.PSM.metadata`` is typed ``dict[str, str]``, so a value written to a TSV file
    and read back is always a string (e.g. ``"False"``), not the ``bool`` mumble may have set
    in memory. Plain truthiness would treat that string as truthy (any non-empty string is),
    silently including mumble candidates in "before" for the standalone-regen path. Only the
    literal string ``"false"`` (any case) counts as excluded; anything else, including an
    absent key, is treated as original.

    """
    value = psm.metadata.get("original_psm", True)
    if isinstance(value, str):
        return value.strip().lower() != "false"
    return bool(value)


def evaluate_before(psm_list: PSMList, config: Dict) -> RescoreResult:
    """
    Evaluate the PSMs' current (pre-rescoring) score with ristretto, for report baselines.

    Called right after :py:func:`~ms2rescore.parse_psms.parse_psms`, before any feature
    generator runs. Excludes mumble-generated alternate candidates (``metadata
    ["original_psm"] is False``) entirely -- they're ms2rescore's own addition, not
    something the original search engine reported, so they shouldn't inflate the "before"
    baseline. Trimmed to ``max_psm_rank_output`` the same way :py:func:`~ms2rescore.rescoring.
    rescore` trims ``after``, so before/after counts stay comparable at the same rank depth.

    """
    original_mask = np.array([_is_original_psm(psm) for psm in psm_list])
    psm_list = psm_list[original_mask]

    train_fdr = config["rescoring"].get("train_fdr", 0.01)
    features_df = _build_features_dataframe(
        psm_list, set(), infer_score_direction(psm_list, train_fdr)
    )
    return _trim_and_evaluate(
        features_df,
        config["max_psm_rank_output"],
        run_col="run",
        peptide_col="peptide",
        protein_col="protein" if "protein" in features_df.columns else None,
        decoy_pattern=config["id_decoy_pattern"],
    )


def evaluate_before_from_provenance(psm_list: PSMList, config: Dict) -> RescoreResult:
    """
    Rebuild the "before" ``RescoreResult`` for standalone report regeneration.

    Uses each PSM's pre-rescoring score, stashed in ``provenance_data`` by
    :py:func:`~ms2rescore.parse_psms.parse_psms`, and otherwise delegates to
    :py:func:`evaluate_before` -- which already excludes mumble-generated candidates and
    trims to ``max_psm_rank_output`` -- so this reconstructs the same population a live run
    would have evaluated, without needing a separately persisted table.

    """
    before_scores = np.array(
        [float(psm.provenance_data["before_rescoring_score"]) for psm in psm_list]
    )
    # Boolean-mask indexing on PSMList makes a new list of the *same* PSM object references,
    # not copies -- mutating .score below would otherwise corrupt the caller's psm_list too
    # (e.g. ReportData.from_files uses the same psm_list again afterwards, for "after").
    psm_list = PSMList(psm_list=[psm.model_copy(deep=True) for psm in psm_list])
    psm_list["score"] = before_scores
    return evaluate_before(psm_list, config)


def evaluate_after_from_psm_list(psm_list: PSMList, config: Dict) -> RescoreResult:
    """
    Rebuild the "after" ``RescoreResult`` for standalone report regeneration.

    Re-competes ``psm_list``'s current (already-rescored) score and identity columns from
    scratch via ristretto -- reproduces the same result :py:func:`~ms2rescore.rescoring.
    rescore` computed during the live run, without needing a separately persisted rollup
    table. ``psm_list``'s score is already in ristretto's "higher is better" convention
    (written back by ``rescore``), and its rows are already trimmed to
    ``max_psm_rank_output``, so :py:func:`_trim_and_evaluate` applies that trim again as a
    no-op.

    """
    features_df = _build_features_dataframe(psm_list, set(), lower_score_is_better=False)
    return _trim_and_evaluate(
        features_df,
        config["max_psm_rank_output"],
        run_col="run",
        peptide_col="peptide",
        protein_col="protein" if "protein" in features_df.columns else None,
        decoy_pattern=config["id_decoy_pattern"],
    )


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


def write_rescoring_tables(after: RescoreResult, output_file_root: str) -> None:
    """
    Write the rescoring result tables to TSV.

    The main user-facing rescoring deliverable -- separate from the full PSM list output
    (which includes rescoring features, provenance data, etc). Only ``after`` is written:
    "before" is fully reconstructable from the PSM list's provenance data (see
    :py:func:`evaluate_before_from_provenance`), so it's never persisted.

    """
    after.psms.to_csv(f"{output_file_root}.psms.tsv", sep="\t", index=False)
    after.peptidoforms.to_csv(f"{output_file_root}.peptidoforms.tsv", sep="\t", index=False)
    if after.peptides is not None:
        after.peptides.to_csv(f"{output_file_root}.peptides.tsv", sep="\t", index=False)
    if after.proteins is not None:
        after.proteins.to_csv(f"{output_file_root}.proteins.tsv", sep="\t", index=False)
    after.feature_weights.to_csv(f"{output_file_root}.weights.tsv", sep="\t")


def _fix_constant_pep_result(
    result: RescoreResult,
    peptide_col: Optional[str],
    protein_col: Optional[str],
    decoy_pattern: Optional[str],
) -> Tuple[RescoreResult, Optional[np.ndarray]]:
    """
    Detect and fix constant PEP (all 1.0) on a single ``RescoreResult``.

    Removes decoy PSMs that score higher than the best target, then recomputes q-values,
    PEPs, and the peptidoform/peptide/protein rollups from scratch on the corrected
    population -- the original per-row values (and rollups derived from them) are products
    of the broken calibration and must not be carried forward.

    Returns
    -------
    tuple[RescoreResult, numpy.ndarray | None]
        The (possibly fixed) result, and the boolean keep-mask applied to ``result.psms`` --
        or ``None`` if no fix was needed/possible, in which case ``result`` is unchanged.

    """
    psms = result.psms
    if not (psms["pep"] == 1.0).all():
        return result, None

    logger.warning(
        "Attempting to fix constant PEP values by removing decoy PSMs that score higher than the "
        "best target PSM."
    )
    max_target_score = psms["score"][~psms["is_decoy"]].max()
    higher_scoring_decoys = psms["is_decoy"].to_numpy() & (
        psms["score"].to_numpy() > max_target_score
    )

    if not higher_scoring_decoys.any():
        logger.warning("No decoys scoring higher than the best target found. Skipping fix.")
        return result, None

    keep = ~higher_scoring_decoys
    logger.warning(f"Removed {higher_scoring_decoys.sum()} decoy PSMs.")
    filtered_psms = psms[keep]
    if filtered_psms["is_decoy"].nunique() < 2:
        # No decoys (or no targets) remain to recompute a competition against -- an extreme
        # edge case of an already-degenerate input. Keep the filtered rows as-is rather than
        # erroring on ristretto.evaluate()'s target/decoy requirement.
        logger.warning(
            "Cannot recompute q-values/PEP/rollups after removing decoys: no decoys remain "
            "in the corrected population. Keeping filtered rows with their prior values."
        )
        return replace(result, psms=filtered_psms), keep

    fixed = ristretto.evaluate(
        filtered_psms,
        run_col="run" if "run" in filtered_psms.columns else None,
        peptide_col=peptide_col,
        protein_col=protein_col,
        decoy_pattern=decoy_pattern,
        multi_rank_rescoring=True,
    )
    # ristretto.evaluate() always returns empty feature_weights/n_iterations (no training
    # occurs there) -- carry over the real values already present on the input result.
    fixed = replace(
        fixed, feature_weights=result.feature_weights, n_iterations=result.n_iterations
    )
    return fixed, keep


def _fix_constant_pep(
    psm_list: PSMList,
    result: RescoreResult,
    peptide_col: Optional[str] = None,
    protein_col: Optional[str] = None,
    decoy_pattern: Optional[str] = None,
) -> Tuple[PSMList, RescoreResult]:
    """
    Workaround for broken PEP calculation if the best-scoring PSM is a decoy.

    Filters ``psm_list`` and ``result`` together (same length, same row order) so they never
    desync -- callers (result tables, report, identification counts) all keep seeing the same
    PSM population as the returned ``psm_list``, with the same freshly recomputed rollups.

    """
    fixed_result, keep = _fix_constant_pep_result(result, peptide_col, protein_col, decoy_pattern)
    if keep is None:
        return psm_list, fixed_result

    psm_list = psm_list[keep]
    psm_list["score"] = fixed_result.psms["score"].to_numpy()
    psm_list["qvalue"] = fixed_result.psms["qvalue"].to_numpy()
    psm_list["pep"] = fixed_result.psms["pep"].to_numpy()
    return psm_list, fixed_result
