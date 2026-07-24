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

MS²Rescore uses the unified ``ms2pip.correlate()`` API with spectra preloaded on
``psm.spectrum``. Both raw ``MS2Spectrum`` and pre-annotated ``AnnotatedMS2Spectrum`` objects
are accepted; ms2pip uses each PSM's own annotation directly, so multi-rank inputs are handled
correctly.

If you use MS²PIP through MS²Rescore, please cite:

.. epigraph::
    Declercq, A., Bouwmeester, R., Chiva, C., Sabidó, E., Hirschler, A., Carapito, C., Martens, L.,
    Degroeve, S., Gabriels, R. Updated MS²PIP web server supports cutting-edge proteomics
    applications. *Nucleic Acids Research* (2023)
    `doi:10.1093/nar/gkad335 <https://doi.org/10.1093/nar/gkad335>`_

"""

import logging

from ms2pip import correlate
from ms2rescore_rs import (
    AnnotatedMS2Spectrum,
    MS2Spectrum,
    ms2pip_features_from_prediction_peak_arrays,
)
from psm_utils import PSMList

from ms2rescore.feature_generators.base import FeatureGeneratorBase, FeatureGeneratorException
from ms2rescore.parse_spectra import MSDataType

logger = logging.getLogger(__name__)


class MS2PIPFeatureGenerator(FeatureGeneratorBase):
    """Generate MS²PIP-based features from spectra preloaded on the input PSMs."""

    required_ms_data = {MSDataType.ms2_spectra}

    def __init__(
        self,
        *args,
        model: str = "HCD",
        model_dir: str | None = None,
        processes: int = 1,
        **kwargs,
    ) -> None:
        """
        Generate MS²PIP-based features.

        Parameters
        ----------
        model
            MS²PIP prediction model to use. Defaults to :py:const:`HCD`.
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
        Add MS²PIP-derived features to PSMs with preloaded spectra.

        Parameters
        ----------
        psm_list
            PSMs to add features to. Each PSM must have a spectrum attached on
            ``psm.spectrum``. Both ``MS2Spectrum`` and ``AnnotatedMS2Spectrum`` are accepted;
            ms2pip uses each PSM's own annotation directly.

        """
        logger.info("Adding MS²PIP-derived features to PSMs.")

        if any(
            not isinstance(psm.spectrum, (MS2Spectrum, AnnotatedMS2Spectrum)) for psm in psm_list
        ):
            raise FeatureGeneratorException(
                "MS²PIP feature generation expects spectra to be preloaded on `psm.spectrum`. "
                "Parse and attach spectra before calling `MS2PIPFeatureGenerator.add_features()`."
            )

        ms2pip_results = correlate(
            psm_list,
            compute_correlations=False,
            model=self.model,
            model_dir=self.model_dir,
            processes=self.processes,
        )
        self._calculate_features(psm_list, ms2pip_results)

    def _calculate_features(self, psm_list, ms2pip_results):
        idx = []
        pred_b = []
        pred_y = []
        obs_b = []
        obs_y = []

        for r in ms2pip_results:
            if r.observed_intensity is None or r.predicted_intensity is None:
                continue
            idx.append(r.psm_index)
            pred_b.append(r.predicted_intensity["b"])
            pred_y.append(r.predicted_intensity["y"])
            obs_b.append(r.observed_intensity["b"])
            obs_y.append(r.observed_intensity["y"])

        results = ms2pip_features_from_prediction_peak_arrays(idx, pred_b, pred_y, obs_b, obs_y)

        for psm_index, feats in results:
            if feats:
                try:
                    psm_list[psm_index]["rescoring_features"].update(feats)
                except (AttributeError, TypeError):
                    psm_list[psm_index]["rescoring_features"] = feats
