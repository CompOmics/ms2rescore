from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from ms2rescore_rs import AnnotatedMS2Spectrum, MS2Spectrum, Precursor, annotate_ms2_spectra
from psm_utils import PSM, PSMList

from ms2rescore.feature_generators.base import FeatureGeneratorException
from ms2rescore.feature_generators.ms2pip import MS2PIPFeatureGenerator
from ms2rescore.parse_spectra import annotate_spectra


def _make_ms2_spectrum(identifier: str = "scan=1") -> MS2Spectrum:
    return MS2Spectrum(
        identifier=identifier,
        mz=[100.0, 200.0],
        intensity=[1000.0, 500.0],
        precursor=Precursor(mz=475.14, charge=2, rt=51.2),
    )


def _make_psm_list(with_spectrum: bool = True) -> PSMList:
    psm = PSM(peptidoform="PEPTIDE/2", spectrum_id="scan=1", run="run1")
    psm.rescoring_features = {}
    if with_spectrum:
        psm.spectrum = _make_ms2_spectrum()
    return PSMList(psm_list=[psm])


def test_ms2pip_feature_generator_uses_unified_correlate(monkeypatch):
    psm_list = _make_psm_list()
    captured = {}

    def fake_correlate(psms, **kwargs):
        captured["psms"] = psms
        captured["kwargs"] = kwargs
        return [
            SimpleNamespace(
                psm_index=0,
                predicted_intensity={"b": np.array([1.0]), "y": np.array([2.0])},
                observed_intensity={"b": np.array([3.0]), "y": np.array([4.0])},
            )
        ]

    def fake_feature_calculation(idx, pred_b, pred_y, obs_b, obs_y):
        captured["feature_inputs"] = (idx, pred_b, pred_y, obs_b, obs_y)
        return [(0, {"spec_pearson_norm": 0.91})]

    monkeypatch.setattr("ms2rescore.feature_generators.ms2pip.correlate", fake_correlate)
    monkeypatch.setattr(
        "ms2rescore.feature_generators.ms2pip.ms2pip_features_from_prediction_peak_arrays",
        fake_feature_calculation,
    )

    feature_generator = MS2PIPFeatureGenerator(model="HCD2021", processes=4)
    feature_generator.add_features(psm_list)

    assert captured["psms"] is psm_list
    assert "spectrum_file" not in captured["kwargs"]
    assert captured["kwargs"]["compute_correlations"] is False
    assert captured["kwargs"]["model"] == "HCD2021"
    assert "ms2_tolerance" not in captured["kwargs"]
    assert captured["kwargs"]["processes"] == 4
    assert psm_list[0].rescoring_features["spec_pearson_norm"] == 0.91


def test_ms2pip_feature_generator_requires_preloaded_spectra(monkeypatch):
    psm_list = _make_psm_list(with_spectrum=False)

    with pytest.raises(FeatureGeneratorException, match="preloaded on `psm.spectrum`"):
        MS2PIPFeatureGenerator().add_features(psm_list)


def _patch_ms2pip_model(monkeypatch):
    """Stub only ms2pip model lookup and XGBoost inference so real correlate/annotation runs."""
    monkeypatch.setattr("ms2pip.core.validate_model", lambda model, model_dir: Path("."))

    def fake_predict_intensities(features, num_ions, model_params, model_dir, **kwargs):
        return [
            {
                "b": np.linspace(0.1, 0.9, n, dtype=np.float32),
                "y": np.linspace(0.9, 0.1, n, dtype=np.float32),
            }
            for n in num_ions
        ]

    monkeypatch.setattr("ms2pip.core.predict_intensities", fake_predict_intensities)


@pytest.mark.parametrize("mode", ["raw", "annotate_spectra", "preannotated_none"])
def test_ms2pip_feature_generator_missing_precursor_matches_controls(monkeypatch, mode):
    """MS2PIPFeatureGenerator handles precursor=None across raw and annotated spectrum inputs."""
    _patch_ms2pip_model(monkeypatch)

    # PEPTIDE/2 matching b2 (227.1026), b3 (324.1554), b4 (425.2031), y3 (362.1558), y4 (463.2035)
    mz = [227.1026, 324.1554, 362.1558, 425.2031, 463.2035]
    intensity = [200.0, 500.0, 350.0, 400.0, 600.0]

    def build_psm_list(precursor):
        spec = MS2Spectrum(identifier="scan=1", mz=mz, intensity=intensity, precursor=precursor)
        if mode == "preannotated_none":
            spec = annotate_ms2_spectra(
                spectra=[spec],
                proformas=["PEPTIDE/2"],
                fragmentation_model="cidhcd",
                mass_mode="monoisotopic",
                tolerance_value=0.02,
                tolerance_mode="Da",
            )[0]
        psm = PSM(peptidoform="PEPTIDE/2", spectrum_id="scan=1", run="run1", spectrum=spec)
        psm.rescoring_features = {}
        pl = PSMList(psm_list=[psm])
        if mode == "annotate_spectra":
            annotate_spectra(pl)
        return pl

    pl_none = build_psm_list(None)
    pl_empty = build_psm_list(Precursor())
    pl_explicit = build_psm_list(Precursor(mz=400.19, charge=2, rt=51.2))

    original_spec = pl_none[0].spectrum

    fgen = MS2PIPFeatureGenerator(model="HCD", processes=1)
    fgen.add_features(pl_none)
    fgen.add_features(pl_empty)
    fgen.add_features(pl_explicit)

    # All MS2PIP features are populated and match controls identically
    feats_none = pl_none[0].rescoring_features
    feats_empty = pl_empty[0].rescoring_features
    feats_explicit = pl_explicit[0].rescoring_features

    assert set(feats_none.keys()) == set(fgen.feature_names)
    assert feats_none == feats_empty == feats_explicit

    # Original spectrum object is preserved on the PSM after feature generation
    assert pl_none[0].spectrum is original_spec
    if mode == "preannotated_none":
        assert isinstance(pl_none[0].spectrum, AnnotatedMS2Spectrum)
        assert pl_none[0].spectrum.precursor is None
    elif mode == "raw":
        assert isinstance(pl_none[0].spectrum, MS2Spectrum)
        assert pl_none[0].spectrum.precursor is None


def test_ms2pip_feature_generator_mixed_precursors_preserves_explicit(monkeypatch):
    """A PSMList mixing precursor=None and explicit Precursor works without altering explicit ones."""
    _patch_ms2pip_model(monkeypatch)

    mz = [227.1026, 324.1554, 362.1558, 425.2031, 463.2035]
    intensity = [200.0, 500.0, 350.0, 400.0, 600.0]
    spec_none = MS2Spectrum(identifier="scan=1", mz=mz, intensity=intensity)
    spec_explicit = MS2Spectrum(
        identifier="scan=2",
        mz=mz,
        intensity=intensity,
        precursor=Precursor(mz=400.19, charge=2, rt=51.2),
    )

    psm_list = PSMList(
        psm_list=[
            PSM(
                peptidoform="PEPTIDE/2",
                spectrum_id="scan=1",
                run="run1",
                spectrum=spec_none,
                rescoring_features={},
            ),
            PSM(
                peptidoform="PEPTIDE/2",
                spectrum_id="scan=2",
                run="run1",
                spectrum=spec_explicit,
                rescoring_features={},
            ),
        ]
    )

    MS2PIPFeatureGenerator(model="HCD", processes=1).add_features(psm_list)

    assert psm_list[0].rescoring_features == psm_list[1].rescoring_features
    assert psm_list[0].spectrum is spec_none
    assert psm_list[0].spectrum.precursor is None
    assert psm_list[1].spectrum is spec_explicit
    assert psm_list[1].spectrum.precursor.mz == pytest.approx(400.19)
    assert psm_list[1].spectrum.precursor.charge == 2
    assert psm_list[1].spectrum.precursor.rt == pytest.approx(51.2)


def test_ms2pip_feature_generator_missing_precursor_and_missing_charge(monkeypatch):
    """When both ProForma charge and precursor are missing, PSM is skipped cleanly (no features)."""
    _patch_ms2pip_model(monkeypatch)

    mz = [227.1026, 324.1554, 425.2031]
    intensity = [200.0, 500.0, 400.0]

    pl_none = PSMList(
        psm_list=[
            PSM(
                peptidoform="PEPTIDE",
                spectrum_id="scan=1",
                run="run1",
                spectrum=MS2Spectrum(identifier="scan=1", mz=mz, intensity=intensity),
                rescoring_features={},
            )
        ]
    )
    pl_empty = PSMList(
        psm_list=[
            PSM(
                peptidoform="PEPTIDE",
                spectrum_id="scan=1",
                run="run1",
                spectrum=MS2Spectrum(
                    identifier="scan=1", mz=mz, intensity=intensity, precursor=Precursor()
                ),
                rescoring_features={},
            )
        ]
    )

    fgen = MS2PIPFeatureGenerator(model="HCD", processes=1)
    fgen.add_features(pl_none)
    fgen.add_features(pl_empty)

    assert pl_none[0].rescoring_features == {}
    assert pl_empty[0].rescoring_features == {}
