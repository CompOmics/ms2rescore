"""Minimal suite for the IM2Deep ion-mobility (CCS) feature generator.

Uses real IM2Deep prediction (bundled model, no network). Marked ``slow`` because it loads
a model.
"""

import math

import numpy as np
import pytest
from psm_utils import PSM, PSMList

from ms2rescore.feature_generators.im2deep import IM2DeepFeatureGenerator
from ms2rescore.parse_spectra import MSDataType

EXPECTED_FEATURE_NAMES = [
    "ccs_observed_im2deep",
    "ccs_predicted_im2deep",
    "ccs_error_im2deep",
    "abs_ccs_error_im2deep",
    "perc_ccs_error_im2deep",
]

_PEPTIDES = [
    "LGGNEQVTR",
    "YILAGVENSK",
    "TPVISGGPYEYR",
    "VEATFGVDESNAK",
    "DAVTPADFSEWSK",
    "GAGSSEPVTGLDAK",
    "LVNELTEFAK",
    "GISNEGQNASIK",
    "YICDNQDTISSK",
    "SAMPLERPEPTIK",
]


def _make_psm_list(is_decoy=False) -> PSMList:
    psms = []
    for i, seq in enumerate(_PEPTIDES):
        psm = PSM(
            peptidoform=f"{seq}/2",
            spectrum_id=f"scan={i + 1}",
            run="run1",
            is_decoy=is_decoy,
            qvalue=0.001,
            score=20.0 - i * 0.1,
            precursor_mz=500.0 + i,
            ion_mobility=0.80 + i * 0.02,
        )
        psm.rescoring_features = {}
        psms.append(psm)
    return PSMList(psm_list=psms)


def test_static_contract():
    """Feature names and required MS data are an external contract; pin them deliberately."""
    assert IM2DeepFeatureGenerator.required_ms_data == {MSDataType.ion_mobility}
    names = IM2DeepFeatureGenerator().feature_names
    assert names == EXPECTED_FEATURE_NAMES
    assert len(names) == len(set(names))


def test_multi_mode_not_implemented():
    with pytest.raises(NotImplementedError, match="Multi-IM"):
        IM2DeepFeatureGenerator(multi=True)


@pytest.mark.slow
def test_prediction_calibrates_onto_observed_scale():
    """Real IM2Deep run: all features populated, ccs_error sign/operands correct, and calibrated
    predictions land on the observed CCS scale (physically plausible residuals)."""
    psm_list = _make_psm_list()
    IM2DeepFeatureGenerator(processes=1).add_features(psm_list)

    for psm in psm_list:
        for name in EXPECTED_FEATURE_NAMES:
            assert name in psm.rescoring_features

    # Per-PSM alignment: ccs_error is predicted - observed (operand order + sign) and must
    # re-derive from each PSM's own values. abs_/perc_ just re-derive from ccs_error, so
    # asserting them would only restate the source; not tested separately.
    for psm in psm_list:
        f = psm.rescoring_features
        assert f["ccs_error_im2deep"] == pytest.approx(
            f["ccs_predicted_im2deep"] - f["ccs_observed_im2deep"]
        )
        assert not math.isnan(f["ccs_predicted_im2deep"])

    # Ion mobility inputs are strictly increasing, so observed CCS must be too. A predict/PSM
    # ordering mismatch would scramble this monotonicity.
    observed = np.array([p.rescoring_features["ccs_observed_im2deep"] for p in psm_list])
    assert np.all(np.diff(observed) > 0)

    perc = np.array([p.rescoring_features["perc_ccs_error_im2deep"] for p in psm_list])
    # Calibration shifts predictions onto the observed CCS scale; residuals stay physically small
    # (empirically median ~8%). Gross breakage (wrong units, model not loaded, calibration
    # skipped) pushes this far higher.
    assert np.median(perc) < 15


@pytest.mark.slow
def test_no_target_psms_raises_on_calibration():
    """A run with only decoys has no PSMs to calibrate against."""
    psm_list = _make_psm_list(is_decoy=True)
    with pytest.raises(ValueError, match="no target PSMs"):
        IM2DeepFeatureGenerator(processes=1).add_features(psm_list)


def test_get_im_calibration_data_excludes_mumble_psms():
    """Mumble-generated candidate PSMs (original_psm=False) are unconfirmed and must never enter
    the calibration set, even if they carry a high-confidence qvalue copied from the original
    hit."""
    import pandas as pd

    run_df = pd.DataFrame(
        {
            "peptidoform": ["A", "B", "C"],
            "ccs_observed_im2deep": [500.0, 999.0, 520.0],
            "qvalue": [0.001, 0.001, 0.005],
            "is_decoy": [False, False, False],
            # "B" is a Mumble mass-shift candidate sharing the original's qvalue.
            "original_psm": [True, False, True],
        }
    )

    gen = IM2DeepFeatureGenerator()
    gen.calibration_set_size = None
    calibration_df = gen._get_im_calibration_data(run_df)

    assert list(calibration_df["peptidoform"]) == ["A", "C"]
    assert list(calibration_df["CCS"]) == [500.0, 520.0]
