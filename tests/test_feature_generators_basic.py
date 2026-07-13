import math

from psm_utils import PSM, PSMList

from ms2rescore.feature_generators.basic import BasicFeatureGenerator

EXPECTED_FEATURE_NAMES = [
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


def _make_psm(peptidoform="PEPTIDEK/2", with_mz=True, with_score=True) -> PSM:
    psm = PSM(peptidoform=peptidoform, spectrum_id="scan=1", run="run1")
    psm.rescoring_features = {}
    if with_mz:
        # Set precursor m/z exactly to the theoretical value -> zero mass error
        psm.precursor_mz = psm.peptidoform.theoretical_mz
    if with_score:
        psm.score = 42.0
    return psm


def test_feature_names_unique_and_expected():
    names = BasicFeatureGenerator().feature_names
    assert names == EXPECTED_FEATURE_NAMES
    assert len(names) == len(set(names))


def test_add_features_populates_all_feature_names():
    psm_list = PSMList(psm_list=[_make_psm()])
    BasicFeatureGenerator().add_features(psm_list)

    features = psm_list[0].rescoring_features
    for name in EXPECTED_FEATURE_NAMES:
        assert name in features


def test_add_features_values():
    psm_list = PSMList(psm_list=[_make_psm(peptidoform="PEPTIDEK/2")])
    BasicFeatureGenerator().add_features(psm_list)
    features = psm_list[0].rescoring_features

    assert features["pep_len"] == 8  # PEPTIDEK
    assert features["charge_n"] == 2
    assert features["charge_2"] == 1
    assert sum(features[f"charge_{i}"] for i in range(1, 7)) == 1  # one-hot
    assert features["search_engine_score"] == 42.0
    # precursor_mz set to theoretical -> essentially no MS1 error
    assert features["abs_ms1_error_ppm"] < 1e-3
    assert abs(features["mass_error"]) < 1e-3
    assert not math.isnan(features["theoretical_mass"])


def test_charge_out_of_range_is_all_zero_one_hot():
    psm_list = PSMList(psm_list=[_make_psm(peptidoform="PEPTIDEK/8")])
    BasicFeatureGenerator().add_features(psm_list)
    features = psm_list[0].rescoring_features

    assert features["charge_n"] == 8
    assert sum(features[f"charge_{i}"] for i in range(1, 7)) == 0


def test_missing_mz_and_score_default_to_zero():
    psm_list = PSMList(psm_list=[_make_psm(with_mz=False, with_score=False)])
    BasicFeatureGenerator().add_features(psm_list)
    features = psm_list[0].rescoring_features

    # Charge comes from the peptidoform, so charge features remain populated
    assert features["charge_2"] == 1
    # m/z- and score-derived features fall back to 0
    assert features["abs_ms1_error_ppm"] == 0
    assert features["experimental_mass"] == 0
    assert features["theoretical_mass"] == 0
    assert features["mass_error"] == 0
    assert features["search_engine_score"] == 0
