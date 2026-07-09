from types import SimpleNamespace

import numpy as np
import pytest
from ms2rescore_rs import MS2Spectrum, Precursor
from psm_utils import PSM, PSMList

from ms2rescore.feature_generators.base import FeatureGeneratorException
from ms2rescore.feature_generators.ms2pip import MS2PIPFeatureGenerator


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
