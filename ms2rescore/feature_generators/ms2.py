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

CHUNK_SIZE = 5000
WHOLE_MOD_TOL = 0.01  # Da, loss mass equal to a modification mass = whole-modification loss

# Kept deliberately small. On a phospho-enriched mumble dataset, n_mods and a matched-loss
# count pushed the model towards unmodified PSMs, and a hyperscore-minus-spectrum-best delta
# was learned with inverted sign (target spectra have one dominant candidate and many far
# losers, decoy spectra are flat). Diagnostic ions were absent for all modifications seen.
MOD_FEATURE_NAMES = [
    "mod_loss_intensity_ratio",
    "precursor_mod_loss_ratio",
    "delta_hyperscore_unmod",
    "mod_site_flank_matched",
    "mod_site_flank_intensity_ratio",
    "mod_mass_error_offset_ppm",
]
MIN_IONS_FOR_OFFSET = 2
N_TERM_SERIES = {"a", "b", "c"}


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
            Add modification-aware features (see :py:const:`MOD_FEATURE_NAMES`). Unmodified
            PSMs get 0 for all of them. Requires
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

        if modified and not any(spectra[i].extended_annotations for i in modified):
            logger.warning(
                "Spectra carry no extended annotations; modification-specific ion features "
                "are all zero. Annotate with extended=True to enable them."
            )

        # ponytail: chunked so transient annotated spectra stay bounded (~900k PSMs at once
        # exhausted 80 GB); chunk size only trades Rust call overhead against peak memory.
        for i in modified:
            matched, intensity = _flank_features(psm_list[i].peptidoform, spectra[i])
            feats["mod_site_flank_matched"][i] = matched
            feats["mod_site_flank_intensity_ratio"][i] = intensity
            feats["mod_mass_error_offset_ppm"][i] = _mass_error_offset(
                psm_list[i].peptidoform, spectra[i]
            )

        full_hs = [psm.rescoring_features["hyperscore"] for psm in psm_list]
        for start in range(0, len(modified), CHUNK_SIZE):
            chunk = modified[start : start + CHUNK_SIZE]
            self._loss_features(psm_list, spectra, chunk, feats)
            self._leave_one_out_delta(psm_list, spectra, chunk, seq_lens, full_hs, feats)

        for i, psm in enumerate(psm_list):
            psm.rescoring_features.update(
                {name: float(feats[name][i]) for name in MOD_FEATURE_NAMES}
            )

    def _loss_features(self, psm_list, spectra, chunk, feats) -> None:
        """Reference pass: same peptidoforms as numeric mass shifts (no modification identity),
        so that every extended annotation not present there is modification-specific."""
        reference = self._annotate(
            [_raw_spectrum(spectra[i]) for i in chunk],
            [proforma_to_mass_shift(psm_list[i].peptidoform) for i in chunk],
            True,
        )
        for i, ref in zip(chunk, reference):
            spec = spectra[i]
            total = float(sum(spec.intensity)) or 1.0
            # A loss of the whole modification (e.g. Sulfo -SO3) leaves the plain unmodified
            # fragment, which is present anyway for fragments not covering the true site.
            # Such ions carry no site or identity evidence and are ignored.
            whole_mod = _mod_masses(psm_list[i].peptidoform)
            generic = {
                (k, a.series, a.position, a.charge, a.neutral_loss)
                for k, anns in enumerate(ref.extended_annotations)
                for a in anns
            }
            loss_peaks, precursor_peaks = set(), set()
            for k, anns in enumerate(spec.extended_annotations):
                for a in anns:
                    if (k, a.series, a.position, a.charge, a.neutral_loss) in generic:
                        continue
                    if any(abs(a.loss_mass - m) < WHOLE_MOD_TOL for m in whole_mod):
                        continue
                    if a.ion_type == "backbone" and a.neutral_loss:
                        loss_peaks.add(k)
                    elif a.ion_type == "precursor" and a.neutral_loss:
                        precursor_peaks.add(k)
            feats["mod_loss_intensity_ratio"][i] = (
                sum(spec.intensity[k] for k in loss_peaks) / total
            )
            feats["precursor_mod_loss_ratio"][i] = (
                sum(spec.intensity[k] for k in precursor_peaks) / total
            )

    def _leave_one_out_delta(self, psm_list, spectra, chunk, seq_lens, full_hs, feats) -> None:
        """Hyperscore gain of the least supported modification over its removal."""
        alt_spectra, alt_proformas, alt_owner = [], [], []
        for i in chunk:
            for proforma in _leave_one_out_proformas(psm_list[i].peptidoform):
                alt_spectra.append(_raw_spectrum(spectra[i]))
                alt_proformas.append(proforma)
                alt_owner.append(i)
        if not alt_owner:
            return
        alt_hs = self._hyperscores(
            self._annotate(alt_spectra, alt_proformas, False), [seq_lens[i] for i in alt_owner]
        )
        delta = defaultdict(list)
        for i, hs in zip(alt_owner, alt_hs):
            delta[i].append(full_hs[i] - hs)
        for i, deltas in delta.items():
            feats["delta_hyperscore_unmod"][i] = min(deltas)


def _flank_features(peptidoform: Peptidoform, spectrum) -> tuple[float, float]:
    """Site-flanking backbone ions for each modification, minimum over modifications.

    For a modification on residue ``s`` (0-based, length ``L``) the flanking ions are the last
    N-terminal ion without the site (b_s), the first with it (b_s+1), and likewise y_L-s-1 and
    y_L-s. Terminal modifications have two. Returns (fraction matched, matched intensity over
    total intensity); a/b/c and x/y/z series are pooled by terminus.
    """
    length = len(peptidoform.parsed_sequence)
    peak_intensity: dict[tuple[str, int], float] = {}
    for k, anns in enumerate(spectrum.peak_annotations):
        for terminus in {("N" if a.series in N_TERM_SERIES else "C", a.position) for a in anns}:
            peak_intensity[terminus] = peak_intensity.get(terminus, 0.0) + spectrum.intensity[k]
    total = float(sum(spectrum.intensity)) or 1.0

    sites = [s for s, (_, mods) in enumerate(peptidoform.parsed_sequence) if mods]
    flanks = [[("N", s), ("N", s + 1), ("C", length - s - 1), ("C", length - s)] for s in sites]
    if peptidoform.properties.get("n_term"):
        flanks.append([("N", 1), ("C", length - 1)])
    if peptidoform.properties.get("c_term"):
        flanks.append([("N", length - 1), ("C", 1)])

    best = (1.0, 1.0)
    for ions in flanks:
        possible = [ion for ion in ions if 1 <= ion[1] <= length - 1]
        hit = [ion for ion in possible if ion in peak_intensity]
        candidate = (
            len(hit) / len(possible) if possible else 0.0,
            sum(peak_intensity[ion] for ion in hit) / total,
        )
        best = min(best, candidate)
    return best if flanks else (0.0, 0.0)


def _mass_error_offset(peptidoform: Peptidoform, spectrum) -> float:
    """|median ppm error of backbone ions containing a modification site - median of the rest|.

    A wrong modification identity of near-identical mass (Phospho vs Sulfo, 9.5 mDa) shifts
    every site-containing fragment by that mass while the other fragments stay put. 0 when
    either group has fewer than MIN_IONS_FOR_OFFSET ions.
    """
    length = len(peptidoform.parsed_sequence)
    sites = [s for s, (_, mods) in enumerate(peptidoform.parsed_sequence) if mods]
    n_min = min(sites) + 1 if sites else length  # first N-terminal ion containing a site
    c_min = length - max(sites) if sites else length  # first C-terminal ion containing a site
    if peptidoform.properties.get("n_term"):
        n_min = 1
    if peptidoform.properties.get("c_term"):
        c_min = 1
    containing, other = [], []
    for k, anns in enumerate(spectrum.peak_annotations):
        mz = spectrum.mz[k]
        for a in anns:
            ppm = a.mz_error / mz * 1e6
            n_series = a.series in N_TERM_SERIES
            has_site = a.position >= (n_min if n_series else c_min)
            (containing if has_site else other).append(ppm)
    if len(containing) < MIN_IONS_FOR_OFFSET or len(other) < MIN_IONS_FOR_OFFSET:
        return 0.0
    return abs(float(np.median(containing) - np.median(other)))


def _mod_masses(peptidoform: Peptidoform) -> list[float]:
    mods = [m for _, ms in peptidoform.parsed_sequence if ms for m in ms]
    mods += peptidoform.properties.get("n_term") or []
    mods += peptidoform.properties.get("c_term") or []
    return [float(m.mass) for m in mods]


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
