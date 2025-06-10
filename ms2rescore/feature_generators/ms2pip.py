"""
MS²PIP fragmentation intensity-based feature generator.

MS²PIP is a machine learning tool that predicts the MS2 spectrum of a peptide given its sequence.
It is previously identified MS2 spectra and their corresponding peptide sequences. Because MS²PIP
uses the highly performant - but traditional - machine learning approach XGBoost, it can already
produce accurate predictions even if trained on smaller spectral libraries. This makes MS²PIP a
very flexible platform to train new models on custom datasets. Nevertheless, MS²PIP comes with
several pre-trained models. See
`github.com/compomics/ms2pip <https://github.com/compomics/ms2pip>`_ for more information.

Because traditional proteomics search engines do not fully consider MS2 peak intensities in their
scoring functions, adding rescoring features derived from spectrum prediction tools has proved to
be a very effective way to further improve the sensitivity of peptide-spectrum matching.

If you use MS²PIP through MS²Rescore, please cite:

.. epigraph::
    Declercq, A., Bouwmeester, R., Chiva, C., Sabidó, E., Hirschler, A., Carapito, C., Martens, L.,
    Degroeve, S., Gabriels, R. Updated MS²PIP web server supports cutting-edge proteomics
    applications. *Nucleic Acids Research* (2023)
    `doi:10.1093/nar/gkad335 <https://doi.org/10.1093/nar/gkad335>`_

"""

import logging
import multiprocessing
import os
import warnings
from itertools import chain
from typing import List, Optional, Union

import numpy as np
import pandas as pd
from ms2pip import correlate
from ms2pip.exceptions import NoMatchingSpectraFound
from ms2pip.result import ProcessingResult
from psm_utils import PSMList
from rich.progress import track

from ms2rescore.feature_generators.base import FeatureGeneratorBase, FeatureGeneratorException
from ms2rescore.parse_spectra import MSDataType
from ms2rescore.utils import infer_spectrum_path

logger = logging.getLogger(__name__)


class MS2PIPFeatureGenerator(FeatureGeneratorBase):
    """Generate MS²PIP-based features."""

    required_ms_data = {MSDataType.ms2_spectra}

    def __init__(
        self,
        *args,
        model: str = "HCD",
        ms2_tolerance: float = 0.02,
        spectrum_path: Optional[str] = None,
        spectrum_id_pattern: str = "(.*)",
        model_dir: Optional[str] = None,
        processes: 1,
        **kwargs,
    ) -> None:
        """
        Generate MS²PIP-based features.

        Parameters
        ----------
        model
            MS²PIP prediction model to use. Defaults to :py:const:`HCD`.
        ms2_tolerance
            MS2 mass tolerance in Da. Defaults to :py:const:`0.02`.
        spectrum_path
            Path to spectrum file or directory with spectrum files. If None, inferred from ``run``
            field in PSMs. Defaults to :py:const:`None`.
        spectrum_id_pattern : str, optional
            Regular expression pattern to extract spectrum ID from spectrum file. Defaults to
            :py:const:`.*`.
        model_dir
            Directory containing MS²PIP models. Defaults to :py:const:`None` (use MS²PIP default).
        processes : int, optional
            Number of processes to use. Defaults to 1.

        Attributes
        ----------
        feature_names: list[str]
            Names of the features that will be added to the PSMs.

        """
        super().__init__(*args, **kwargs)
        self.model = model
        self.ms2_tolerance = ms2_tolerance
        self.spectrum_path = spectrum_path
        self.spectrum_id_pattern = spectrum_id_pattern
        self.model_dir = model_dir
        self.processes = processes

    @property
    def feature_names(self):
        return [
            "spec_pearson_norm",
            "ionb_pearson_norm",
            "iony_pearson_norm",
            "spec_mse_norm",
            "ionb_mse_norm",
            "iony_mse_norm",
            "min_abs_diff_norm",
            "max_abs_diff_norm",
            "abs_diff_Q1_norm",
            "abs_diff_Q2_norm",
            "abs_diff_Q3_norm",
            "mean_abs_diff_norm",
            "std_abs_diff_norm",
            "ionb_min_abs_diff_norm",
            "ionb_max_abs_diff_norm",
            "ionb_abs_diff_Q1_norm",
            "ionb_abs_diff_Q2_norm",
            "ionb_abs_diff_Q3_norm",
            "ionb_mean_abs_diff_norm",
            "ionb_std_abs_diff_norm",
            "iony_min_abs_diff_norm",
            "iony_max_abs_diff_norm",
            "iony_abs_diff_Q1_norm",
            "iony_abs_diff_Q2_norm",
            "iony_abs_diff_Q3_norm",
            "iony_mean_abs_diff_norm",
            "iony_std_abs_diff_norm",
            "dotprod_norm",
            "dotprod_ionb_norm",
            "dotprod_iony_norm",
            "cos_norm",
            "cos_ionb_norm",
            "cos_iony_norm",
            "spec_pearson",
            "ionb_pearson",
            "iony_pearson",
            "spec_spearman",
            "ionb_spearman",
            "iony_spearman",
            "spec_mse",
            "ionb_mse",
            "iony_mse",
            "min_abs_diff_iontype",
            "max_abs_diff_iontype",
            "min_abs_diff",
            "max_abs_diff",
            "abs_diff_Q1",
            "abs_diff_Q2",
            "abs_diff_Q3",
            "mean_abs_diff",
            "std_abs_diff",
            "ionb_min_abs_diff",
            "ionb_max_abs_diff",
            "ionb_abs_diff_Q1",
            "ionb_abs_diff_Q2",
            "ionb_abs_diff_Q3",
            "ionb_mean_abs_diff",
            "ionb_std_abs_diff",
            "iony_min_abs_diff",
            "iony_max_abs_diff",
            "iony_abs_diff_Q1",
            "iony_abs_diff_Q2",
            "iony_abs_diff_Q3",
            "iony_mean_abs_diff",
            "iony_std_abs_diff",
            "dotprod",
            "dotprod_ionb",
            "dotprod_iony",
            "cos",
            "cos_ionb",
            "cos_iony",
        ]

    def add_features(self, psm_list: PSMList) -> None:
        """
        Add MS²PIP-derived features to PSMs.

        Parameters
        ----------
        psm_list
            PSMs to add features to.

        """
        logger.info("Adding MS²PIP-derived features to PSMs.")
        psm_dict = psm_list.get_psm_dict()
        current_run = 1
        total_runs = sum(len(runs) for runs in psm_dict.values())

        for runs in psm_dict.values():
            for run, psms in runs.items():
                logger.info(
                    f"Running MS²PIP for PSMs from run ({current_run}/{total_runs}) `{run}`..."
                )
                psm_list_run = PSMList(psm_list=list(chain.from_iterable(psms.values())))
                spectrum_filename = infer_spectrum_path(self.spectrum_path, run)
                logger.debug(f"Using spectrum file `{spectrum_filename}`")
                try:
                    os.environ.pop("CUDA_VISIBLE_DEVICES", None)
                    ms2pip_results = correlate(
                        psms=psm_list_run,
                        spectrum_file=str(spectrum_filename),
                        spectrum_id_pattern=self.spectrum_id_pattern,
                        model=self.model,
                        ms2_tolerance=self.ms2_tolerance,
                        compute_correlations=False,
                        model_dir=self.model_dir,
                        processes=self.processes,
                    )
                except NoMatchingSpectraFound as e:
                    raise FeatureGeneratorException(
                        f"Could not find any matching spectra for PSMs from run `{run}`. "
                        "Please check that the `spectrum_id_pattern` and `psm_id_pattern` "
                        "options are configured correctly. See "
                        "https://ms2rescore.readthedocs.io/en/latest/userguide/configuration/#mapping-psms-to-spectra"
                        " for more information."
                    ) from e
                self._calculate_features(psm_list_run, ms2pip_results)
                current_run += 1

    def _calculate_features(
        self, psm_list: PSMList, ms2pip_results: List[ProcessingResult]
    ) -> None:
        """Calculate features from all MS²PIP results and add to PSMs."""
        logger.debug("Calculating features from predicted spectra")
        with multiprocessing.Pool(int(self.processes)) as pool:
            # Use imap, so we can use a progress bar
            counts_failed = 0
            for result, features in zip(
                ms2pip_results,
                track(
                    pool.imap(self._calculate_features_single, ms2pip_results, chunksize=1000),
                    total=len(ms2pip_results),
                    description="Calculating features...",
                    transient=True,
                ),
            ):
                if features:
                    # Cannot use result.psm directly, as it is a copy from MS²PIP multiprocessing
                    try:
                        psm_list[result.psm_index]["rescoring_features"].update(features)
                    except (AttributeError, TypeError):
                        psm_list[result.psm_index]["rescoring_features"] = features
                else:
                    counts_failed += 1

        if counts_failed > 0:
            logger.warning(f"Failed to calculate features for {counts_failed} PSMs")

    def _calculate_features_single(self, processing_result: ProcessingResult) -> Union[dict, None]:
        """Calculate MS²PIP-based features for single PSM."""
        if (
            processing_result.observed_intensity is None
            or processing_result.predicted_intensity is None
        ):
            return None

        # Suppress RuntimeWarnings about invalid values
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            # Convert intensities to arrays
            target_b = processing_result.predicted_intensity["b"].clip(np.log2(0.001))
            target_y = processing_result.predicted_intensity["y"].clip(np.log2(0.001))
            target_all = np.concatenate([target_b, target_y])
            prediction_b = processing_result.observed_intensity["b"].clip(np.log2(0.001))
            prediction_y = processing_result.observed_intensity["y"].clip(np.log2(0.001))
            prediction_all = np.concatenate([prediction_b, prediction_y])

            # Prepare 'unlogged' intensity arrays
            target_b_unlog = 2**target_b - 0.001
            target_y_unlog = 2**target_y - 0.001
            target_all_unlog = 2**target_all - 0.001
            prediction_b_unlog = 2**prediction_b - 0.001
            prediction_y_unlog = 2**prediction_y - 0.001
            prediction_all_unlog = 2**prediction_all - 0.001

            # Calculate absolute differences
            abs_diff_b = np.abs(target_b - prediction_b)
            abs_diff_y = np.abs(target_y - prediction_y)
            abs_diff_all = np.abs(target_all - prediction_all)
            abs_diff_b_unlog = np.abs(target_b_unlog - prediction_b_unlog)
            abs_diff_y_unlog = np.abs(target_y_unlog - prediction_y_unlog)
            abs_diff_all_unlog = np.abs(target_all_unlog - prediction_all_unlog)

            # Compute features
            feature_values = [
                # Features between spectra in log space
                np.corrcoef(target_all, prediction_all)[0][1],  # Pearson all ions
                np.corrcoef(target_b, prediction_b)[0][1],  # Pearson b ions
                np.corrcoef(target_y, prediction_y)[0][1],  # Pearson y ions
                _mse(target_all, prediction_all),  # MSE all ions
                _mse(target_b, prediction_b),  # MSE b ions
                _mse(target_y, prediction_y),  # MSE y ions
                np.min(abs_diff_all),  # min_abs_diff_norm
                np.max(abs_diff_all),  # max_abs_diff_norm
                np.quantile(abs_diff_all, 0.25),  # abs_diff_Q1_norm
                np.quantile(abs_diff_all, 0.5),  # abs_diff_Q2_norm
                np.quantile(abs_diff_all, 0.75),  # abs_diff_Q3_norm
                np.mean(abs_diff_all),  # mean_abs_diff_norm
                np.std(abs_diff_all),  # std_abs_diff_norm
                np.min(abs_diff_b),  # ionb_min_abs_diff_norm
                np.max(abs_diff_b),  # ionb_max_abs_diff_norm
                np.quantile(abs_diff_b, 0.25),  # ionb_abs_diff_Q1_norm
                np.quantile(abs_diff_b, 0.5),  # ionb_abs_diff_Q2_norm
                np.quantile(abs_diff_b, 0.75),  # ionb_abs_diff_Q3_norm
                np.mean(abs_diff_b),  # ionb_mean_abs_diff_norm
                np.std(abs_diff_b),  # ionb_std_abs_diff_norm
                np.min(abs_diff_y),  # iony_min_abs_diff_norm
                np.max(abs_diff_y),  # iony_max_abs_diff_norm
                np.quantile(abs_diff_y, 0.25),  # iony_abs_diff_Q1_norm
                np.quantile(abs_diff_y, 0.5),  # iony_abs_diff_Q2_norm
                np.quantile(abs_diff_y, 0.75),  # iony_abs_diff_Q3_norm
                np.mean(abs_diff_y),  # iony_mean_abs_diff_norm
                np.std(abs_diff_y),  # iony_std_abs_diff_norm
                np.dot(target_all, prediction_all),  # Dot product all ions
                np.dot(target_b, prediction_b),  # Dot product b ions
                np.dot(target_y, prediction_y),  # Dot product y ions
                _cosine_similarity(target_all, prediction_all),  # Cos similarity all ions
                _cosine_similarity(target_b, prediction_b),  # Cos similarity b ions
                _cosine_similarity(target_y, prediction_y),  # Cos similarity y ions
                # Same features in normal space
                np.corrcoef(target_all_unlog, prediction_all_unlog)[0][1],  # Pearson all
                np.corrcoef(target_b_unlog, prediction_b_unlog)[0][1],  # Pearson b
                np.corrcoef(target_y_unlog, prediction_y_unlog)[0][1],  # Pearson y
                _spearman(target_all_unlog, prediction_all_unlog),  # Spearman all ions
                _spearman(target_b_unlog, prediction_b_unlog),  # Spearman b ions
                _spearman(target_y_unlog, prediction_y_unlog),  # Spearman y ions
                _mse(target_all_unlog, prediction_all_unlog),  # MSE all ions
                _mse(target_b_unlog, prediction_b_unlog),  # MSE b ions
                _mse(target_y_unlog, prediction_y_unlog),  # MSE y ions,
                # Ion type with min absolute difference
                0 if np.min(abs_diff_b_unlog) <= np.min(abs_diff_y_unlog) else 1,
                # Ion type with max absolute difference
                0 if np.max(abs_diff_b_unlog) >= np.max(abs_diff_y_unlog) else 1,
                np.min(abs_diff_all_unlog),  # min_abs_diff
                np.max(abs_diff_all_unlog),  # max_abs_diff
                np.quantile(abs_diff_all_unlog, 0.25),  # abs_diff_Q1
                np.quantile(abs_diff_all_unlog, 0.5),  # abs_diff_Q2
                np.quantile(abs_diff_all_unlog, 0.75),  # abs_diff_Q3
                np.mean(abs_diff_all_unlog),  # mean_abs_diff
                np.std(abs_diff_all_unlog),  # std_abs_diff
                np.min(abs_diff_b_unlog),  # ionb_min_abs_diff
                np.max(abs_diff_b_unlog),  # ionb_max_abs_diff_norm
                np.quantile(abs_diff_b_unlog, 0.25),  # ionb_abs_diff_Q1
                np.quantile(abs_diff_b_unlog, 0.5),  # ionb_abs_diff_Q2
                np.quantile(abs_diff_b_unlog, 0.75),  # ionb_abs_diff_Q3
                np.mean(abs_diff_b_unlog),  # ionb_mean_abs_diff
                np.std(abs_diff_b_unlog),  # ionb_std_abs_diff
                np.min(abs_diff_y_unlog),  # iony_min_abs_diff
                np.max(abs_diff_y_unlog),  # iony_max_abs_diff
                np.quantile(abs_diff_y_unlog, 0.25),  # iony_abs_diff_Q1
                np.quantile(abs_diff_y_unlog, 0.5),  # iony_abs_diff_Q2
                np.quantile(abs_diff_y_unlog, 0.75),  # iony_abs_diff_Q3
                np.mean(abs_diff_y_unlog),  # iony_mean_abs_diff
                np.std(abs_diff_y_unlog),  # iony_std_abs_diff
                np.dot(target_all_unlog, prediction_all_unlog),  # Dot product all ions
                np.dot(target_b_unlog, prediction_b_unlog),  # Dot product b ions
                np.dot(target_y_unlog, prediction_y_unlog),  # Dot product y ions
                _cosine_similarity(target_all_unlog, prediction_all_unlog),  # Cos similarity all
                _cosine_similarity(target_b_unlog, prediction_b_unlog),  # Cos similarity b ions
                _cosine_similarity(target_y_unlog, prediction_y_unlog),  # Cos similarity y ions
            ]

        features = dict(
            zip(
                self.feature_names,
                [0.0 if np.isnan(ft) else ft for ft in feature_values],
            )
        )

        return features


def _spearman(x: np.ndarray, y: np.ndarray) -> float:
    """Spearman rank correlation."""
    x = np.array(x)
    y = np.array(y)
    x_rank = pd.Series(x).rank()
    y_rank = pd.Series(y).rank()
    return np.corrcoef(x_rank, y_rank)[0][1]


def _mse(x: np.ndarray, y: np.ndarray) -> float:
    """Mean squared error"""
    x = np.array(x)
    y = np.array(y)
    return np.mean((x - y) ** 2)


def _cosine_similarity(x: np.ndarray, y: np.ndarray) -> float:
    """Cosine similarity"""
    x = np.array(x)
    y = np.array(y)
    return np.dot(x, y) / (np.linalg.norm(x, 2) * np.linalg.norm(y, 2))
