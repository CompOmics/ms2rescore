"""
MS2-based feature generator.

"""

import logging

from typing import List, Optional

from psm_utils import PSMList
from ms2rescore_rs import ms2_features_from_ms2spectra
from ms2rescore.feature_generators.base import FeatureGeneratorBase

logger = logging.getLogger(__name__)


class MS2FeatureGenerator(FeatureGeneratorBase):
    """DeepLC retention time-based feature generator."""

    def __init__(
        self,
        *args,
        spectrum_path: Optional[str] = None,
        spectrum_id_pattern: str = "(.*)",
        fragmentation_model: str = "cidhcd",
        mass_mode: str = "monoisotopic",
        processes: int = 1,
        calculate_hyperscore: bool = False,
        **kwargs,
    ) -> None:
        """
        Generate MS2-based features for rescoring.

        Parameters
        ----------

        Attributes
        ----------
        feature_names: list[str]
            Names of the features that will be added to the PSMs.

        """
        super().__init__(*args, **kwargs)

        self.spectrum_path = spectrum_path
        self.spectrum_id_pattern = spectrum_id_pattern
        self.fragmentation_model = fragmentation_model.lower()
        self.mass_mode = mass_mode.lower()
        self.processes = processes
        self.calculate_hyperscore = calculate_hyperscore

    @property
    def feature_names(self) -> List[str]:
        return [
            "ln_explained_intensity",
            "ln_total_intensity",
            "ln_explained_intensity_ratio",
            "ln_explained_b_ion_ratio",
            "ln_explained_y_ion_ratio",
            "longest_b_ion_sequence",
            "longest_y_ion_sequence",
            "matched_b_ions",
            "matched_b_ions_pct",
            "matched_y_ions",
            "matched_y_ions_pct",
            "matched_ions_pct",
            "hyperscore",
        ]

    def add_features(self, psm_list: PSMList) -> None:
        logger.info("Adding MS2-derived features to PSMs.")

        spectra = psm_list["spectrum"]
        # Keep parity with your current behavior:
        proformas = [psm.peptidoform.proforma.split("/")[0] for psm in psm_list]
        seq_lens = [len(psm.peptidoform.sequence) for psm in psm_list]

        feature_dicts = ms2_features_from_ms2spectra(
            spectra=spectra,
            proformas=proformas,
            seq_lens=seq_lens,
            fragmentation_model=self.fragmentation_model,
            mass_mode=self.mass_mode,
            calculate_hyperscore=self.calculate_hyperscore,
        )

        for psm, feats in zip(psm_list, feature_dicts):
            psm.rescoring_features.update(feats)
