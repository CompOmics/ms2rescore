"""Generate basic features that can be extracted from any PSM list."""

from __future__ import annotations

import logging

import numpy as np
from psm_utils import PSMList

from ms2rescore.feature_generators.base import FeatureGeneratorBase

logger = logging.getLogger(__name__)


class BasicFeatureGenerator(FeatureGeneratorBase):
    def __init__(self, *args, **kwargs) -> None:
        """
        Generate basic features that can be extracted from any PSM list, including search engine
        score, charge state, and MS1 error.

        Parameters
        ----------
        *args
            Positional arguments passed to the base class.
        **kwargs
            Keyword arguments passed to the base class.

        Attributes
        ----------
        feature_names: list[str]
            Names of the features that will be added to the PSMs.

        """
        super().__init__(*args, **kwargs)

    @property
    def feature_names(self) -> list[str]:
        return [
            "charge_n",
            "charge_1",
            "charge_2",
            "charge_3",
            "charge_4",
            "charge_5",
            "charge_6",
            "abs_ms1_error_ppm",
            "search_engine_score",
            "theoretical_mass",
            "experimental_mass",
            "mass_error",
            "pep_len",
        ]

    def add_features(self, psm_list: PSMList) -> None:
        """
        Add basic features to a PSM list.

        All features listed in ``feature_names`` are always added. When the required data
        (charge, m/z, score) is not available, the corresponding features are set to 0.

        Parameters
        ----------
        psm_list
            PSM list to add features to.

        """
        logger.info("Adding basic features to PSMs.")
        n = len(psm_list)

        charge_states = np.array([psm.peptidoform.precursor_charge for psm in psm_list])
        precursor_mzs = psm_list["precursor_mz"]
        scores = psm_list["score"]
        peptide_lengths = np.array([len(psm.peptidoform.sequence) for psm in psm_list])

        has_charge = None not in charge_states
        # precursor_mz and score come back as float arrays where missing values are NaN, not None
        has_mz = not np.isnan(precursor_mzs).any() and has_charge
        has_score = not np.isnan(scores).any()

        if has_charge:
            charge_n = charge_states
            charge_one_hot, _ = _one_hot_encode_charge(charge_states)
        else:
            logger.warning("Charge states not available for all PSMs; charge features will be 0.")
            charge_n = np.zeros(n)
            charge_one_hot = [dict.fromkeys([f"charge_{i}" for i in range(1, 7)], 0) for _ in range(n)]

        if has_mz:  # Charge also required for theoretical m/z
            theo_mz = np.array([psm.peptidoform.theoretical_mz for psm in psm_list])
            abs_ms1_error_ppm = np.abs((precursor_mzs - theo_mz) / theo_mz * 10**6)
            experimental_mass = (precursor_mzs * charge_n) - (charge_n * 1.007276466812)
            theoretical_mass = (theo_mz * charge_n) - (charge_n * 1.007276466812)
            mass_error = experimental_mass - theoretical_mass
        else:
            logger.warning("Precursor m/z not available for all PSMs; m/z features will be 0.")
            abs_ms1_error_ppm = np.zeros(n)
            experimental_mass = np.zeros(n)
            theoretical_mass = np.zeros(n)
            mass_error = np.zeros(n)

        for i, psm in enumerate(psm_list):
            if psm.rescoring_features is None:
                psm.rescoring_features = {}
            psm.rescoring_features.update(
                {
                    "charge_n": charge_n[i],
                    **charge_one_hot[i],
                    "abs_ms1_error_ppm": abs_ms1_error_ppm[i],
                    "search_engine_score": scores[i] if has_score else 0,
                    "theoretical_mass": theoretical_mass[i],
                    "experimental_mass": experimental_mass[i],
                    "mass_error": mass_error[i],
                    "pep_len": peptide_lengths[i],
                }
            )


def _one_hot_encode_charge(
    charge_states: np.ndarray,
) -> tuple[list[dict[str, int]], list[str]]:
    """One-hot encode charge states with fixed range 1-6.

    Charge states outside the 1-6 range are encoded as all zeros.
    """
    n_entries = len(charge_states)
    heading = [f"charge_{i}" for i in range(1, 7)]

    # Create mask for charges 1-6
    mask = np.zeros((n_entries, 6), dtype=bool)

    # Set the appropriate charge position to 1 for each entry
    for i, charge in enumerate(charge_states):
        if charge is not None and 1 <= charge <= 6:
            mask[i, int(charge) - 1] = 1

    one_hot = mask.view("i1")

    return [dict(zip(heading, row)) for row in one_hot], heading
