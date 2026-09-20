import math

from ms2rescore_rs import AnnotatedMS2Spectrum, FragmentAnnotation, Precursor
from psm_utils import PSM, Peptidoform, PSMList

from ms2rescore.feature_generators.ms2 import MS2FeatureGenerator
from ms2rescore.parse_spectra import MSDataType


def _make_annotated_spectrum(identifier: str = "scan=1") -> AnnotatedMS2Spectrum:
    """Three-peak spectrum with a b1 and y1 annotation and one unmatched peak."""
    peak_annotations = [
        [FragmentAnnotation(series="b", position=1, charge=1)],
        [FragmentAnnotation(series="y", position=1, charge=1)],
        [],
    ]
    return AnnotatedMS2Spectrum(
        identifier=identifier,
        mz=[100.0, 200.0, 300.0],
        intensity=[1000.0, 500.0, 250.0],
        precursor=Precursor(mz=475.14, charge=2, rt=51.2),
        peak_annotations=peak_annotations,
    )


def _make_psm_list(peptidoform: str = "PEPTIDE/2", identifier: str = "scan=1") -> PSMList:
    psm = PSM(peptidoform=peptidoform, spectrum_id=identifier, run="run1")
    psm.rescoring_features = {}
    psm.spectrum = _make_annotated_spectrum(identifier)
    return PSMList(psm_list=[psm])


def test_required_ms_data():
    assert MS2FeatureGenerator.required_ms_data == {MSDataType.ms2_spectra}


def test_feature_names_default():
    generator = MS2FeatureGenerator()
    names = generator.feature_names
    # No duplicates and hyperscore present
    assert len(names) == len(set(names))
    assert "hyperscore" in names
    assert "matched_ions_pct" in names
    # Per-series features exist for all six primary ion series
    for series in ["a", "b", "c", "x", "y", "z"]:
        assert f"ln_explained_{series}_ion_ratio" in names
        assert f"longest_{series}_ion_sequence" in names
        assert f"matched_{series}_ions" in names
        assert f"matched_{series}_ions_pct" in names


def test_add_features_populates_all_feature_names():
    psm_list = _make_psm_list()
    generator = MS2FeatureGenerator()

    generator.add_features(psm_list)

    features = psm_list[0].rescoring_features
    for name in generator.feature_names:
        assert name in features


def test_add_features_scores_matched_ions():
    """cidhcd activates a, b, y; the b1 and y1 peaks must be matched."""
    psm_list = _make_psm_list()
    generator = MS2FeatureGenerator(fragmentation_model="cidhcd")

    generator.add_features(psm_list)

    features = psm_list[0].rescoring_features
    assert features["matched_b_ions"] == 1.0
    assert features["matched_y_ions"] == 1.0
    assert features["hyperscore"] > 0
    # Inactive series for cidhcd yield NaN
    for series in ["c", "x", "z"]:
        assert math.isnan(features[f"matched_{series}_ions"])


def test_fragmentation_model_selects_active_ion_series():
    """etd activates c, y, z; a, b, x become NaN even though a b peak is present."""
    psm_list = _make_psm_list()
    generator = MS2FeatureGenerator(fragmentation_model="etd")

    generator.add_features(psm_list)

    features = psm_list[0].rescoring_features
    for series in ["c", "y", "z"]:
        assert not math.isnan(features[f"matched_{series}_ions"])
    for series in ["a", "b", "x"]:
        assert math.isnan(features[f"matched_{series}_ions"])


# ---------------------------------------------------------------- add_mod_info features
from ms2rescore_rs import MS2Spectrum

from ms2rescore.feature_generators.ms2 import (
    MOD_FEATURE_NAMES,
    _leave_one_out_proformas,
)
from ms2rescore.parse_spectra import annotate_spectra

# PEPS[Phospho]TIDE/2: b3 324.1554, b4-H3PO4 393.1769, precursor-H3PO4 435.1980 (2+), y5 644.2175,
# plus b2 227.1026 and y3 376.1714 (no site) and b4 491.1538 (site) for the mass-error offset
PHOSPHO_MZ = [227.1026, 324.1554, 376.1714, 393.1769, 435.1980, 491.1538, 644.2175, 999.0]
PHOSPHO_INTENSITY = [10.0, 10.0, 10.0, 20.0, 30.0, 10.0, 40.0, 70.0]  # total 200


def _phospho_psm_list() -> PSMList:
    def psm(peptidoform, spectrum_id):
        p = PSM(peptidoform=peptidoform, spectrum_id=spectrum_id, run="run1")
        p.rescoring_features = {}
        p.spectrum = MS2Spectrum(
            identifier=spectrum_id,
            mz=PHOSPHO_MZ,
            intensity=PHOSPHO_INTENSITY,
            precursor=Precursor(mz=484.19, charge=2),
        )
        return p

    return PSMList(
        psm_list=[
            psm("PEPS[Phospho]TIDE/2", "scan=1"),  # correct, modified
            psm("PEPSTIDE/2", "scan=1"),  # unmodified competitor on same spectrum
            psm("PEPS[Sulfo]TIDE/2", "scan=1"),  # isobaric wrong identity, 9.5 mDa off
        ]
    )


def test_mod_features_phospho():
    psm_list = _phospho_psm_list()
    annotate_spectra(psm_list, "cidhcd", 0.02, "Da", extended=True)
    generator = MS2FeatureGenerator(add_mod_info=True, tolerance_value=0.02, tolerance_mode="Da")
    generator.add_features(psm_list)
    mod, unmod, sulfo = (psm.rescoring_features for psm in psm_list)

    assert math.isclose(mod["mod_loss_intensity_ratio"], 20 / 200)
    assert math.isclose(mod["precursor_mod_loss_ratio"], 30 / 200)
    assert mod["delta_hyperscore_unmod"] > 0  # y5 carries the phospho, lost when removed
    # site S4 (0-based 3) of 8: flanking b3, b4, y4, y5; b3, b4 and y5 are present
    assert mod["mod_site_flank_matched"] == 0.75
    assert math.isclose(mod["mod_site_flank_intensity_ratio"], 60 / 200)
    # Phospho fragments sit on their theoretical m/z, Sulfo site fragments are 9.5 mDa off
    assert mod["mod_mass_error_offset_ppm"] < 2
    assert 12 < sulfo["mod_mass_error_offset_ppm"] < 22
    for name in MOD_FEATURE_NAMES:
        assert unmod[name] == 0.0, name
    # Sulfo -SO3 leaves the plain fragment: excluded as evidence, so no loss features
    assert sulfo["mod_loss_intensity_ratio"] == 0.0 and sulfo["precursor_mod_loss_ratio"] == 0.0
    assert math.isclose(sulfo["mod_site_flank_matched"], 0.75)  # same peaks match within 0.02 Da


def test_mod_flag_off_adds_nothing():
    """Without `add_mod_info` no modification feature is declared, computed or annotated."""
    assert not set(MOD_FEATURE_NAMES) & set(MS2FeatureGenerator().feature_names)
    psm_list = _phospho_psm_list()
    annotate_spectra(psm_list, "cidhcd", 20.0, "ppm")
    MS2FeatureGenerator().add_features(psm_list)
    assert not set(MOD_FEATURE_NAMES) & set(psm_list[0].rescoring_features)
    assert psm_list[0].spectrum.extended_annotations == []


def test_leave_one_out_proformas():
    variants = list(
        _leave_one_out_proformas(Peptidoform("[Acetyl]-PEPS[Phospho]TIDE-[Amidated]/2"))
    )
    assert variants == [
        "PEPS[+79.9663]TIDE-[-0.9840]/2",
        "[+42.0106]-PEPS[+79.9663]TIDE/2",
        "[+42.0106]-PEPSTIDE-[-0.9840]/2",
    ]
