import math

import pytest
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
