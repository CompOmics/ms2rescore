"""
Mokapot integration for MS²Rescore.

:py:mod:`mokapot` is a full-Python implementation of the semi-supervised learning algorithms
introduced with Percolator. It builds upon the flexible scikit-learn package, which makes it
highly efficient for routine applications, but also customizable for experimental research
settings. Using Mokapot through MS²Rescore brings several advantages over Percolator: It can be
easily installed in the same Python environment, and it is generally faster as the communication
between the tools happens completely within Python, without the need to write and read files
or communicate through the command line. See
`mokapot.readthedocs.io <https://mokapot.readthedocs.io/>`_ for more information.

If you use Mokapot through MS²Rescore, please cite:

.. epigraph::
   Fondrie W. E. & Noble W. S. mokapot: Fast and Flexible Semisupervised
   Learning for Peptide Detection. *J Proteome Res* (2021).
   `doi:10.1021/acs.jproteome.0c01010 <https://doi.org/10.1021/acs.jproteome.0c01010>`_

"""

import logging
from concurrent.futures import BrokenExecutor
from typing import Any, Dict, List, Optional, Tuple

import mokapot
import numpy as np
import pandas as pd
import psm_utils
from mokapot.brew import brew
from mokapot.dataset import LinearPsmDataset
from mokapot.model import PercolatorModel
from pyteomics.mass import nist_mass

from ms2rescore.constants import CHARGE_PATTERN
from ms2rescore.exceptions import RescoringError

logger = logging.getLogger(__name__)
logging.getLogger("numba").setLevel(logging.WARNING)


def rescore(
    psm_list: psm_utils.PSMList,
    output_file_root: str = "ms2rescore",
    fasta_file: Optional[str] = None,
    train_fdr: float = 0.01,
    write_weights: bool = False,
    write_txt: bool = False,
    protein_kwargs: Optional[Dict[str, Any]] = None,
    **kwargs: Any,
) -> None:
    """
    Rescore PSMs with Mokapot.

    The function provides a high-level interface to use Mokapot within MS²Rescore. It first
    converts the :py:class:`~psm_utils.psm_list.PSMList` to a
    :py:class:`~mokapot.dataset.LinearPsmDataset`, and then optionally adds protein information
    from a FASTA file. The dataset is then passed to the :py:func:`~mokapot.brew` function, which
    returns the new scores, q-values, and PEPs. These are then written back to the original
    :py:class:`~psm_utils.psm_list.PSMList`.

    Parameters
    ----------
    psm_list
        PSMs to be rescored.
    output_file_root
        Root of output file names. Defaults to ``"ms2rescore"``.
    fasta_file
        Path to FASTA file with protein sequences to use for protein inference. Defaults to
        ``None``.
    train_fdr
        FDR to use for training the Mokapot model. Defaults to ``0.01``.
    write_weights
        Write model weights to a text file. Defaults to ``False``.
    write_txt
        Write Mokapot results to a text file. Defaults to ``False``.
    protein_kwargs
        Keyword arguments to pass to the :py:meth:`~mokapot.dataset.LinearPsmDataset.add_proteins`
        method.
    **kwargs
        Additional keyword arguments are passed to the Mokapot :py:func:`~mokapot.brew` function.

    """
    _set_log_levels()

    if "write_flashlfq" in kwargs:
        _ = kwargs.pop("write_flashlfq")
        logger.warning(
            "The `write_flashlfq` argument has moved. To write FlashLFQ generic TSV, use the "
            "MS²Rescore-level `write_flashlfq` option instead."
        )

    # Convert PSMList to Mokapot dataset
    lin_psm_data = convert_psm_list(psm_list)
    feature_names = list(lin_psm_data.features.columns)

    # Add proteins
    if fasta_file:
        logger.debug(f"Adding protein info from {fasta_file} with options: `{protein_kwargs}`")
        lin_psm_data.add_proteins(fasta_file, **protein_kwargs)

    # Rescore
    logger.debug(f"Mokapot brew options: `{kwargs}`")
    try:
        confidence_results, models = brew(
            lin_psm_data, model=PercolatorModel(train_fdr=train_fdr), rng=8, **kwargs
        )
    except (RuntimeError, BrokenExecutor) as e:
        raise RescoringError("Mokapot could not be run. Please check the input data.") from e

    add_psm_confidence(psm_list, confidence_results)
    add_peptide_confidence(psm_list, confidence_results)

    # Write results
    if write_weights:
        try:
            save_model_weights(models, feature_names, output_file_root)
        except AttributeError:
            logger.warning(
                "Could not extract Mokapot model weights with the `coef_` attribute. Most likely, "
                "a model type different from the default (linear SVM) was used. No weights will "
                "be saved."
            )
    if write_txt:
        confidence_results.to_txt(file_root=output_file_root, decoys=True)


def convert_psm_list(
    psm_list: psm_utils.PSMList,
    feature_names: Optional[List[str]] = None,
) -> LinearPsmDataset:
    """
    Convert a PSM list to a Mokapot dataset.

    Parameters
    ----------
    psm_list
        PSMList to rescore.
    feature_names
        List of feature names to use. Items must be keys in the PSM `rescoring_features` dict.

    """
    psm_df = psm_list.to_dataframe()
    psm_df = psm_df.reset_index(drop=True).reset_index()

    psm_df["peptide"] = (
        psm_df["peptidoform"].astype(str).str.replace(CHARGE_PATTERN, "", n=1, regex=True)
    )
    psm_df["is_target"] = ~psm_df["is_decoy"]
    psm_df["charge"] = psm_df["peptidoform"].apply(lambda x: x.precursor_charge)
    psm_df["calcmass"] = psm_df["peptidoform"].apply(lambda x: x.theoretical_mass)
    psm_df["expmass"] = _mz_to_mass(psm_df["precursor_mz"], psm_df["charge"])

    required_columns = [
        "index",
        "spectrum_id",
        "peptide",
        "is_target",
        "protein_list",
        "run",
        "calcmass",
        "expmass",
        "retention_time",
        "charge",
    ]
    feature_df = pd.DataFrame(list(psm_df["rescoring_features"])).astype(float).fillna(0.0)
    feature_df.columns = [f"feature:{f}" for f in feature_df.columns]
    combined_df = pd.concat([psm_df[required_columns], feature_df], axis=1)

    feature_names = [f"feature:{f}" for f in feature_names] if feature_names else None

    lin_psm_data = LinearPsmDataset(
        psms=combined_df,
        target_column="is_target",
        spectrum_columns="index",  # Use artificial index to allow multi-rank rescoring
        peptide_column="peptide",
        protein_column="protein_list",
        feature_columns=feature_names or list(feature_df.columns),
        filename_column="run",
        scan_column="spectrum_id",  # Keep as spectrum_id?
        calcmass_column="calcmass",
        expmass_column="expmass",
        rt_column="retention_time",
        charge_column="charge",
    )

    return lin_psm_data


def save_model_weights(
    models: Tuple[mokapot.model.Model], feature_names: List[str], output_file_root: str
):
    """
    Save model weights to a file.

    Parameters
    ----------
    models
        Tuple of Mokapot models (one for each fold) to save.
    feature_names
        List of feature names that were used to train the models.
    output_file_root
        Root of output file names.

    """
    try:
        coefficients = np.stack([m.estimator.coef_[0] for m in models])
    except AttributeError as e:
        raise AttributeError(
            "Could not extract Mokapot model weights with the `coef_` attribute. Most likely, "
            "a model type different from the default (linear SVM) was used."
        ) from e

    pd.DataFrame(coefficients, columns=list(feature_names)).to_csv(
        output_file_root + ".mokapot.weights.tsv", sep="\t", index=False
    )


def add_psm_confidence(
    psm_list: psm_utils.PSMList, confidence_results: mokapot.confidence.Confidence
) -> None:
    """Add PSM-level confidence estimates to PSM list, updating score, qvalue, pep, and rank."""
    # Reshape confidence estimates to match PSMList
    keys = ["mokapot score", "mokapot q-value", "mokapot PEP"]
    mokapot_values_targets = (
        confidence_results.confidence_estimates["psms"].set_index("index").sort_index()[keys]
    )
    mokapot_values_decoys = (
        confidence_results.decoy_confidence_estimates["psms"].set_index("index").sort_index()[keys]
    )
    q = np.full((len(psm_list), 3), np.nan)
    q[mokapot_values_targets.index] = mokapot_values_targets.values
    q[mokapot_values_decoys.index] = mokapot_values_decoys.values

    # Add Mokapot results to PSMList
    psm_list["score"] = q[:, 0]
    psm_list["qvalue"] = q[:, 1]
    psm_list["pep"] = q[:, 2]

    # Reset ranks to match new scores
    psm_list.set_ranks(lower_score_better=False)


def add_peptide_confidence(
    psm_list: psm_utils.PSMList, confidence_results: mokapot.confidence.Confidence
) -> None:
    """Add Mokapot peptide-level confidence estimates to PSM list."""
    keys = ["mokapot score", "mokapot q-value", "mokapot PEP"]
    peptide_info = pd.concat(
        [
            confidence_results.confidence_estimates["peptides"].set_index("peptide")[keys],
            confidence_results.decoy_confidence_estimates["peptides"].set_index("peptide")[keys],
        ],
        axis=0,
    ).to_dict(orient="index")

    # Add peptide-level scores to PSM metadata
    # run_key = "na" if not all(psm.run for psm in psm_list) else None
    for psm in psm_list:
        peptide_scores = peptide_info[(CHARGE_PATTERN.sub("", str(psm.peptidoform), 1))]
        psm.metadata.update(
            {
                "peptide_score": peptide_scores["mokapot score"],
                "peptide_qvalue": peptide_scores["mokapot q-value"],
                "peptide_pep": peptide_scores["mokapot PEP"],
            }
        )


def _mz_to_mass(mz: float, charge: int) -> float:
    """Convert m/z to mass."""
    return mz * charge - charge * nist_mass["H"][1][0]


def _set_log_levels() -> None:
    """Set log levels for Mokapot and Numba to avoid too-high verbosity."""
    # Set mokapot logging to WARNING if not in debug mode
    if logger.getEffectiveLevel() > logging.DEBUG:
        logging.getLogger("mokapot").setLevel(logging.WARNING)

    # Keep Numba logging to INFO or higher
    if logger.getEffectiveLevel() < logging.INFO:
        logging.getLogger("numba").setLevel(logging.INFO)
