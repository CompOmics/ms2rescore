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

import numpy as np
import pandas as pd
import ristretto
from psm_utils import PSMList
from ristretto import RescoreResult

from ms2rescore._ristretto_utils import _build_features_dataframe, _trim_and_evaluate
from ms2rescore.exceptions import RescoringError
from ms2rescore.parse_psms import infer_score_direction

logger = logging.getLogger(__name__)


def rescore(psm_list: PSMList, config: dict, output_file_root: str) -> tuple[PSMList, RescoreResult]:
    """
    Rescore PSMs with ristretto and write the new scores, q-values, and PEPs back to ``psm_list``.

    Always scores every candidate rank with ristretto's semi-supervised model
    (``multi_rank_rescoring=True``), then runs the actual competition/FDR/PEP step
    separately via :py:func:`~ms2rescore._ristretto_utils._trim_and_evaluate`, trimming to
    ``max_psm_rank_output`` PSMs per spectrum.

    Returns
    -------
    tuple[PSMList, RescoreResult]
        The (possibly filtered) PSMList with updated scores, and the competition/FDR result.

    """
    feature_names = {f for psm in psm_list for f in psm.rescoring_features}
    lower_score_is_better = infer_score_direction(
        psm_list, config["rescoring"].get("train_fdr", 0.01)
    )
    features_df = _build_features_dataframe(psm_list, feature_names, lower_score_is_better)

    peptide_col = "peptide"
    protein_col = "protein" if "protein" in features_df.columns else None
    decoy_pattern = config["id_decoy_pattern"]

    try:
        ml_result = ristretto.rescore(
            features_df,
            run_col="run",
            peptide_col=peptide_col,
            protein_col=protein_col,
            feature_cols=sorted(feature_names),
            decoy_pattern=decoy_pattern,
            n_jobs=int(config["processes"]),
            multi_rank_rescoring=True,
            # train_fdr/model; ristretto's own kwarg defaults apply to any key a partial
            # user-provided rescoring dict omitted.
            **config["rescoring"],
        )
        final_result = _trim_and_evaluate(
            ml_result.psms,
            config["max_psm_rank_output"],
            run_col="run",
            peptide_col=peptide_col,
            protein_col=protein_col,
            decoy_pattern=decoy_pattern,
        )
        # _trim_and_evaluate() calls ristretto.evaluate(), which always returns empty
        # feature_weights/n_iterations (no training occurs there) -- carry over the real
        # values from ml_result, the actual ristretto.rescore() call above.
        final_result = replace(
            final_result,
            feature_weights=ml_result.feature_weights,
            n_iterations=ml_result.n_iterations,
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

    psm_list, final_result = _fix_constant_pep(
        psm_list, final_result, peptide_col, protein_col, decoy_pattern
    )

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

    return psm_list, final_result


def _write_group_metadata(
    psm_list: PSMList,
    psms: pd.DataFrame,
    rollup: pd.DataFrame,
    group_col: str,
    decoy_pattern: str | None = None,
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


def _fix_constant_pep_result(
    result: RescoreResult,
    peptide_col: str | None,
    protein_col: str | None,
    decoy_pattern: str | None,
) -> tuple[RescoreResult, np.ndarray | None]:
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
    peptide_col: str | None = None,
    protein_col: str | None = None,
    decoy_pattern: str | None = None,
) -> tuple[PSMList, RescoreResult]:
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
