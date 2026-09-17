"""
MS2-based feature generator.

"""

import logging
from collections import defaultdict
from copy import deepcopy
from typing import ClassVar

import numpy as np
from ms2pip._spectrum_processing import proforma_to_mass_shift
from ms2rescore_rs import MS2Spectrum, annotate_ms2_spectra, score_ms2_spectra
from psm_utils import Peptidoform, PSMList

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

MOD_FEATURE_NAMES = [
    "n_mods",
    "mod_loss_n_matched",
    "mod_loss_intensity_ratio",
    "precursor_mod_loss_ratio",
    "diagnostic_ion_ratio",
    "delta_hyperscore_unmod",
    "delta_hyperscore_vs_spectrum_best",
]


class MS2FeatureGenerator(FeatureGeneratorBase):
    """MS2 spectrum-based feature generator."""

    required_ms_data: ClassVar[set[MSDataType]] = {MSDataType.ms2_spectra}

    def __init__(
        self,
        *args,
        fragmentation_model: str = "cidhcd",
        add_mod_info: bool = False,
        tolerance_value: float = 0.02,
        tolerance_mode: str = "Da",
        **kwargs,
    ) -> None:
        """
        Generate MS2-based features for rescoring.

        Parameters
        ----------
        fragmentation_model
            Fragmentation model, used to determine active ion series for scoring. Defaults to
            :py:const:`cidhcd` (a, b, and y ions).
        add_mod_info
            Add modification-aware features (see :py:const:`MOD_FEATURE_NAMES`). Requires
            spectra annotated with ``extended=True`` (see
            :py:func:`ms2rescore.parse_spectra.annotate_spectra`).
        tolerance_value, tolerance_mode
            Fragment mass tolerance, only used for the extra annotation passes behind
            ``add_mod_info``.

        Attributes
        ----------
        feature_names: list[str]
            Names of the features that will be added to the PSMs.

        """
        super().__init__(*args, **kwargs)

        self.fragmentation_model = fragmentation_model.lower()
        self.add_mod_info = add_mod_info
        self.tolerance_value = tolerance_value
        self.tolerance_mode = tolerance_mode

    @property
    def feature_names(self) -> list[str]:
        names = [
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
        if self.add_mod_info:
            names += MOD_FEATURE_NAMES
        return names

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

        if self.add_mod_info:
            self._add_mod_features(psm_list, seq_lens)

    # ------------------------------------------------------------------ mod features

    def _annotate(self, spectra, proformas, extended):
        return annotate_ms2_spectra(
            spectra=spectra,
            proformas=proformas,
            fragmentation_model=self.fragmentation_model,
            mass_mode="monoisotopic",
            tolerance_value=self.tolerance_value,
            tolerance_mode=self.tolerance_mode,
            extended=extended,
        )

    def _hyperscores(self, annotated, seq_lens) -> list[float]:
        feats = score_ms2_spectra(
            spectra=annotated,
            seq_lens=seq_lens,
            active_ion_series=ACTIVE_SERIES[self.fragmentation_model],
            calculate_hyperscore=True,
        )
        return [f.get("hyperscore", 0.0) for f in feats]

    def _add_mod_features(self, psm_list: PSMList, seq_lens: list[int]) -> None:
        logger.info("Adding modification-aware MS2 features to PSMs.")
        spectra = list(psm_list["spectrum"])
        n = len(psm_list)
        n_mods = [_count_mods(psm.peptidoform) for psm in psm_list]
        modified = [i for i in range(n) if n_mods[i] > 0]
        feats = {name: np.zeros(n) for name in MOD_FEATURE_NAMES}
        feats["n_mods"] = np.asarray(n_mods, dtype=float)

        if modified and not any(spectra[i].extended_annotations for i in modified):
            logger.warning(
                "Spectra carry no extended annotations; modification-specific ion features "
                "are all zero. Annotate with extended=True to enable them."
            )

        # Reference pass: same peptidoforms as numeric mass shifts (no modification identity),
        # so that every extended annotation not present here is modification-specific.
        raw = [_raw_spectrum(spectra[i]) for i in modified]
        reference = self._annotate(
            raw, [proforma_to_mass_shift(psm_list[i].peptidoform) for i in modified], True
        )
        for i, ref in zip(modified, reference):
            spec = spectra[i]
            total = float(sum(spec.intensity)) or 1.0
            generic = {
                (k, a.series, a.position, a.charge, a.neutral_loss)
                for k, anns in enumerate(ref.extended_annotations)
                for a in anns
            }
            loss_peaks, precursor_peaks, diagnostic_peaks = set(), set(), set()
            for k, anns in enumerate(spec.extended_annotations):
                for a in anns:
                    if (k, a.series, a.position, a.charge, a.neutral_loss) in generic:
                        continue
                    if a.ion_type == "backbone" and a.neutral_loss:
                        loss_peaks.add(k)
                    elif a.ion_type == "precursor" and a.neutral_loss:
                        precursor_peaks.add(k)
                    elif a.ion_type in ("diagnostic", "immonium"):
                        diagnostic_peaks.add(k)
            feats["mod_loss_n_matched"][i] = len(loss_peaks)
            feats["mod_loss_intensity_ratio"][i] = (
                sum(spec.intensity[k] for k in loss_peaks) / total
            )
            feats["precursor_mod_loss_ratio"][i] = (
                sum(spec.intensity[k] for k in precursor_peaks) / total
            )
            feats["diagnostic_ion_ratio"][i] = (
                sum(spec.intensity[k] for k in diagnostic_peaks) / total
            )

        # Leave-one-out pass: hyperscore gain of the least supported modification.
        alt_spectra, alt_proformas, alt_owner = [], [], []
        for i in modified:
            for proforma in _leave_one_out_proformas(psm_list[i].peptidoform):
                alt_spectra.append(_raw_spectrum(spectra[i]))
                alt_proformas.append(proforma)
                alt_owner.append(i)
        if alt_owner:
            alt_hs = self._hyperscores(
                self._annotate(alt_spectra, alt_proformas, False), [seq_lens[i] for i in alt_owner]
            )
            full_hs = [psm.rescoring_features["hyperscore"] for psm in psm_list]
            delta = defaultdict(list)
            for i, hs in zip(alt_owner, alt_hs):
                delta[i].append(full_hs[i] - hs)
            for i, deltas in delta.items():
                feats["delta_hyperscore_unmod"][i] = min(deltas)

        feats["delta_hyperscore_vs_spectrum_best"] = _delta_vs_spectrum_best(psm_list)

        for i, psm in enumerate(psm_list):
            psm.rescoring_features.update(
                {name: float(feats[name][i]) for name in MOD_FEATURE_NAMES}
            )


def _count_mods(peptidoform: Peptidoform) -> int:
    n = sum(len(mods) for _, mods in peptidoform.parsed_sequence if mods)
    n += len(peptidoform.properties.get("n_term") or [])
    n += len(peptidoform.properties.get("c_term") or [])
    return n


def _raw_spectrum(spectrum) -> MS2Spectrum:
    return MS2Spectrum(
        identifier=spectrum.identifier,
        mz=spectrum.mz,
        intensity=spectrum.intensity,
        precursor=spectrum.precursor,
    )


def _leave_one_out_proformas(peptidoform: Peptidoform):
    """Yield mass-shift ProForma strings with one modification removed at a time."""
    for terminus in ("n_term", "c_term"):
        mods = peptidoform.properties.get(terminus) or []
        for k in range(len(mods)):
            variant = deepcopy(peptidoform)
            variant.properties[terminus] = mods[:k] + mods[k + 1 :] or None
            yield proforma_to_mass_shift(variant)
    for i, (aa, mods) in enumerate(peptidoform.parsed_sequence):
        for k in range(len(mods or [])):
            variant = deepcopy(peptidoform)
            variant.parsed_sequence[i] = (aa, mods[:k] + mods[k + 1 :] or None)
            yield proforma_to_mass_shift(variant)


def _delta_vs_spectrum_best(psm_list: PSMList) -> np.ndarray:
    """Hyperscore minus the best hyperscore of any other PSM on the same spectrum; 0 if alone."""
    hs = np.asarray([psm.rescoring_features["hyperscore"] for psm in psm_list], dtype=float)
    groups = defaultdict(list)
    for i, (run, spectrum_id) in enumerate(zip(psm_list["run"], psm_list["spectrum_id"])):
        groups[(str(run), str(spectrum_id))].append(i)
    delta = np.zeros(len(psm_list))
    for idx in groups.values():
        if len(idx) < 2:
            continue
        values = hs[idx]
        top = np.sort(values)[::-1]
        for i, v in zip(idx, values):
            other_best = (
                top[0] if (v < top[0] or (top[0] == v and (top == v).sum() > 1)) else top[1]
            )
            delta[i] = v - other_best
    return delta
