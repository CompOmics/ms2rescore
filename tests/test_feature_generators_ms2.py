import math

from ms2rescore_rs import AnnotatedMS2Spectrum, FragmentAnnotation, Precursor
from psm_utils import PSM, PSMList

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
from ms2rescore_rs import MS2Spectrum  # noqa: E402

from ms2rescore.feature_generators.ms2 import (  # noqa: E402
    MOD_FEATURE_NAMES,
    _delta_vs_spectrum_best,
    _leave_one_out_proformas,
)
from ms2rescore.parse_spectra import annotate_spectra  # noqa: E402

# PEPS[Phospho]TIDE/2: b3 324.1554, b4-H3PO4 393.1769, precursor-H3PO4 435.1980 (2+), y5 644.2175
PHOSPHO_MZ = [324.1554, 393.1769, 435.1980, 644.2175, 999.0]
PHOSPHO_INTENSITY = [10.0, 20.0, 30.0, 40.0, 100.0]  # total 200


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
            psm("PEPSTIDE/2", "scan=2"),  # unmodified, alone on its spectrum
        ]
    )


def test_feature_names_mod_flag():
    assert not set(MOD_FEATURE_NAMES) & set(MS2FeatureGenerator().feature_names)
    names = MS2FeatureGenerator(add_mod_info=True).feature_names
    assert names[-len(MOD_FEATURE_NAMES) :] == MOD_FEATURE_NAMES
    assert len(names) == len(set(names))


def test_mod_features_phospho():
    psm_list = _phospho_psm_list()
    annotate_spectra(psm_list, "cidhcd", 20.0, "ppm", extended=True)
    generator = MS2FeatureGenerator(add_mod_info=True, tolerance_value=20.0, tolerance_mode="ppm")
    generator.add_features(psm_list)
    mod, unmod, alone = (psm.rescoring_features for psm in psm_list)

    assert mod["n_mods"] == 1 and unmod["n_mods"] == 0
    assert mod["mod_loss_n_matched"] == 1
    assert math.isclose(mod["mod_loss_intensity_ratio"], 20 / 200)
    assert math.isclose(mod["precursor_mod_loss_ratio"], 30 / 200)
    assert mod["diagnostic_ion_ratio"] == 0.0
    assert mod["delta_hyperscore_unmod"] > 0  # y5 carries the phospho, lost when removed
    for name in MOD_FEATURE_NAMES:
        if name != "delta_hyperscore_vs_spectrum_best":
            assert unmod[name] == 0.0, name
    # Competition on scan=1: modified PSM explains more (y5), unmodified loses by the same amount
    assert mod["delta_hyperscore_vs_spectrum_best"] > 0
    assert math.isclose(
        unmod["delta_hyperscore_vs_spectrum_best"], -mod["delta_hyperscore_vs_spectrum_best"]
    )
    assert alone["delta_hyperscore_vs_spectrum_best"] == 0.0


def test_mod_flag_off_adds_nothing():
    psm_list = _phospho_psm_list()
    annotate_spectra(psm_list, "cidhcd", 20.0, "ppm")
    MS2FeatureGenerator().add_features(psm_list)
    assert not set(MOD_FEATURE_NAMES) & set(psm_list[0].rescoring_features)
    assert psm_list[0].spectrum.extended_annotations == []


def test_mod_helpers():
    from psm_utils import Peptidoform

    variants = list(
        _leave_one_out_proformas(Peptidoform("[Acetyl]-PEPS[Phospho]TIDE-[Amidated]/2"))
    )
    assert variants == [
        "PEPS[+79.9663]TIDE-[-0.9840]/2",
        "[+42.0106]-PEPS[+79.9663]TIDE/2",
        "[+42.0106]-PEPSTIDE-[-0.9840]/2",
    ]

    psms = PSMList(
        psm_list=[PSM(peptidoform="PEPTIDE", spectrum_id="s1", run="r") for _ in range(3)]
    )
    for psm, hs in zip(psms, [5.0, 5.0, 2.0]):
        psm.rescoring_features = {"hyperscore": hs}
    assert list(_delta_vs_spectrum_best(psms)) == [0.0, 0.0, -3.0]  # tie at the top gives 0
