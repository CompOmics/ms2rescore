"""
Rescore PSMs with ristretto.

:py:mod:`ristretto` is a lean, dependency-light (numpy/scikit-learn/pandas only)
reimplementation of the Percolator/Käll 2007 semi-supervised rescoring algorithm, purpose-built
to serve MS²Rescore. It takes a plain :py:class:`pandas.DataFrame` in and returns a
:py:class:`~ristretto.result.RescoreResult` of plain DataFrames out.

"""

import logging
from concurrent.futures import BrokenExecutor
from dataclasses import replace
from typing import Dict, Tuple

import numpy as np
import ristretto
from psm_utils import PSMList
from ristretto import RescoreResult

from ms2rescore._ristretto_utils import (
    _build_features_dataframe,
    _fix_constant_pep,
    _trim_and_evaluate,
    _write_group_metadata,
)
from ms2rescore.exceptions import RescoringError
from ms2rescore.parse_psms import infer_score_direction

logger = logging.getLogger(__name__)


def rescore(psm_list: PSMList, config: Dict, output_file_root: str) -> Tuple[PSMList, RescoreResult]:
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
