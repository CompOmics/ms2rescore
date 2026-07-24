"""
MS2-based feature generator.

"""

import logging
from typing import ClassVar

from ms2rescore_rs import score_ms2_spectra
from psm_utils import PSMList

from ms2rescore.feature_generators.base import FeatureGeneratorBase
from ms2rescore.parse_spectra import MSDataType

logger = logging.getLogger(__name__)

# Map fragmentation model to active ion series
ACTIVE_SERIES = {
    "cidhcd": ["a", "b", "y"],
    "etd": ["c", "y", "z"],
    "ethcd": ["a", "b", "c", "y", "z"],
    "all": ["a", "b", "c", "x", "y", "z"],
}


class MS2FeatureGenerator(FeatureGeneratorBase):
    """MS2 spectrum-based feature generator."""

    required_ms_data: ClassVar[set[MSDataType]] = {MSDataType.ms2_spectra}

    def __init__(
        self,
        *args,
        fragmentation_model: str = "cidhcd",
        **kwargs,
    ) -> None:
        """
        Generate MS2-based features for rescoring.

        Parameters
        ----------
        fragmentation_model
            Fragmentation model, used to determine active ion series for scoring. Defaults to
            :py:const:`cidhcd` (a, b, and y ions).

        Attributes
        ----------
        feature_names: list[str]
            Names of the features that will be added to the PSMs.

        """
        super().__init__(*args, **kwargs)

        self.fragmentation_model = fragmentation_model.lower()

    @property
    def feature_names(self) -> list[str]:
        return [
            "ln_explained_intensity",
            "ln_total_intensity",
            "ln_explained_intensity_ratio",
            "ln_explained_a_ion_ratio",
            "ln_explained_b_ion_ratio",
            "ln_explained_c_ion_ratio",
            "ln_explained_x_ion_ratio",
            "ln_explained_y_ion_ratio",
            "ln_explained_z_ion_ratio",
            "longest_a_ion_sequence",
            "longest_b_ion_sequence",
            "longest_c_ion_sequence",
            "longest_x_ion_sequence",
            "longest_y_ion_sequence",
            "longest_z_ion_sequence",
            "matched_a_ions",
            "matched_a_ions_pct",
            "matched_b_ions",
            "matched_b_ions_pct",
            "matched_c_ions",
            "matched_c_ions_pct",
            "matched_x_ions",
            "matched_x_ions_pct",
            "matched_y_ions",
            "matched_y_ions_pct",
            "matched_z_ions",
            "matched_z_ions_pct",
            "matched_ions_pct",
            "hyperscore",
        ]

    def add_features(self, psm_list: PSMList) -> None:
        logger.info("Adding MS2-derived features to PSMs.")

        seq_lens = [len(psm.peptidoform.sequence) for psm in psm_list]

        feature_dicts = score_ms2_spectra(
            spectra=list(psm_list["spectrum"]),
            seq_lens=seq_lens,
            active_ion_series=ACTIVE_SERIES[self.fragmentation_model],
            calculate_hyperscore=True,
        )

        for psm, feats in zip(psm_list, feature_dicts):
            psm.rescoring_features.update(feats)
