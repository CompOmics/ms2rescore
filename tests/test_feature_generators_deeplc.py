"""Minimal suite for the DeepLC retention-time feature generator.

Uses real DeepLC prediction (bundled model, no network) with retraining disabled to keep
the run fast and deterministic in structure. Marked ``slow`` because it loads a model.
"""

import numpy as np
import pytest
from psm_utils import PSM, PSMList

from ms2rescore.feature_generators.deeplc import DeepLCFeatureGenerator
from ms2rescore.parse_spectra import MSDataType

EXPECTED_FEATURE_NAMES = [
    "observed_retention_time",
    "predicted_retention_time",
    "rt_diff",
    "observed_retention_time_best",
    "predicted_retention_time_best",
    "rt_diff_best",
]

# A handful of tryptic peptides with monotonic dummy retention times, single run.
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


def _make_psm_list(run="run1", is_decoy=False) -> PSMList:
    psms = []
    for i, seq in enumerate(_PEPTIDES):
        psm = PSM(
            peptidoform=f"{seq}/2",
            spectrum_id=f"scan={i + 1}",
            run=run,
            is_decoy=is_decoy,
            qvalue=0.001,
            score=20.0 - i * 0.1,
            retention_time=100.0 + i * 30.0,
        )
        psm.rescoring_features = {}
        psms.append(psm)
    return PSMList(psm_list=psms)


def test_static_contract():
    """Feature names and required MS data are an external contract; pin them deliberately."""
    assert DeepLCFeatureGenerator.required_ms_data == {MSDataType.retention_time}
    names = DeepLCFeatureGenerator(deeplc_retrain=False).feature_names
    assert names == EXPECTED_FEATURE_NAMES
    assert len(names) == len(set(names))


@pytest.mark.slow
def test_prediction_calibrates_onto_observed_scale():
    """Real DeepLC run: all features populated, rt_diff arithmetic correct, and calibration
    actually maps predictions onto the observed RT scale (small residuals, not collapsed)."""
    psm_list = _make_psm_list()
    DeepLCFeatureGenerator(deeplc_retrain=False, processes=1).add_features(psm_list)

    for psm in psm_list:
        for name in EXPECTED_FEATURE_NAMES:
            assert name in psm.rescoring_features

    # Per-PSM alignment: each PSM must keep its own observed RT and rt_diff must re-derive from
    # that same PSM's values. A predict/PSM ordering mismatch would break the passthrough here.
    for i, psm in enumerate(psm_list):
        f = psm.rescoring_features
        assert f["observed_retention_time"] == 100.0 + i * 30.0
        assert f["rt_diff"] == pytest.approx(
            abs(f["observed_retention_time"] - f["predicted_retention_time"])
        )

    predicted = np.array([p.rescoring_features["predicted_retention_time"] for p in psm_list])
    rt_diffs = np.array([p.rescoring_features["rt_diff"] for p in psm_list])
    # Calibration fits predicted -> observed on these high-confidence PSMs, so residuals stay
    # small (observed span 270; empirically median ~6). A broken/no-op calibration would leave
    # predictions in raw model units and blow this up.
    assert np.median(rt_diffs) < 60
    # Calibration must not collapse every prediction onto a single value.
    assert predicted.std() > 1.0


@pytest.mark.slow
def test_no_target_psms_raises_on_calibration():
    """A run with only decoys has no PSMs to calibrate against."""
    psm_list = _make_psm_list(is_decoy=True)
    with pytest.raises(ValueError, match="no target PSMs"):
        DeepLCFeatureGenerator(deeplc_retrain=False, processes=1).add_features(psm_list)


def test_best_run_by_shared_proteoforms_picks_max_overlap():
    runs = ["A", "A", "B", "B", "C"]
    proteoforms = ["p1", "p2", "p1", "p2", "p3"]
    # A and B share p1 and p2; C shares nothing -> A (first of the tied max) wins
    best = DeepLCFeatureGenerator._best_run_by_shared_proteoforms(runs, proteoforms)
    assert best == "A"


def test_best_run_by_shared_proteoforms_no_overlap_returns_first():
    runs = ["A", "B", "C"]
    proteoforms = ["p1", "p2", "p3"]
    best = DeepLCFeatureGenerator._best_run_by_shared_proteoforms(runs, proteoforms)
    assert best == "A"
