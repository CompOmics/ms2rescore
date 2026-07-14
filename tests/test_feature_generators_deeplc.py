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


def _fake_predict_factory(matrix):
    """Return a stand-in for ``deeplc.core.predict`` that yields a fixed multitask matrix.

    Lets head-selection tests run deterministically without loading the real model. The generator
    always calls ``predict(..., return_matrix=True)``; the ``return_matrix=False`` branch mirrors
    DeepLC's own head-0 default so a regression to the old behaviour is exercised too.
    """
    arr = np.asarray(matrix, dtype=float)

    def _fake_predict(psm_list, model=None, predict_kwargs=None, return_matrix=False):
        return arr if return_matrix else arr[:, 0]

    return _fake_predict


def _make_multirun_psm_list(runs_obs) -> PSMList:
    """Build a PSMList from ``{run: [observed_rt, ...]}`` in the given per-run order."""
    psms = []
    scan = 0
    for run, observed in runs_obs.items():
        for rt in observed:
            scan += 1
            psm = PSM(
                peptidoform=f"{_PEPTIDES[(scan - 1) % len(_PEPTIDES)]}/2",
                spectrum_id=f"scan={scan}",
                run=run,
                is_decoy=False,
                qvalue=0.001,
                score=20.0,
                retention_time=float(rt),
            )
            psm.rescoring_features = {}
            psms.append(psm)
    return PSMList(psm_list=psms)


def test_selects_best_correlating_head(monkeypatch):
    """The best-correlating head must be chosen, not head 0.

    Head 0 is near-uncorrelated noise, head 1 tracks observed RT, head 2 is constant (skipped).
    Only head-1 selection yields small residuals; a revert to head-0-only would collapse
    calibration and blow up rt_diff.
    """
    psm_list = _make_psm_list()
    obs = np.array([100.0 + i * 30.0 for i in range(len(_PEPTIDES))])
    n = len(obs)
    mat = np.zeros((n, 3))
    mat[:, 0] = np.tile([0.0, 1.0], n)[:n]  # non-zero std, ~zero correlation with obs
    mat[:, 1] = obs + 0.5  # best head
    mat[:, 2] = 5.0  # constant -> zero std -> skipped by head selection
    monkeypatch.setattr(
        "ms2rescore.feature_generators.deeplc.predict", _fake_predict_factory(mat)
    )

    DeepLCFeatureGenerator(deeplc_retrain=False, processes=1).add_features(psm_list)

    preds = np.array([p.rescoring_features["predicted_retention_time"] for p in psm_list])
    rt_diffs = np.array([p.rescoring_features["rt_diff"] for p in psm_list])
    assert np.corrcoef(preds, obs)[0, 1] > 0.99
    assert np.median(rt_diffs) < 5.0


def test_per_run_head_selection(monkeypatch):
    """Different runs may select different heads; each must calibrate onto its own scale.

    Head 1 tracks observed RT for run1 only (constant for run2); head 2 tracks run2 only. Correct
    per-run selection gives small residuals in both runs.
    """
    obs1 = [100.0 + i * 30.0 for i in range(10)]
    obs2 = [50.0 + i * 20.0 for i in range(10)]
    psm_list = _make_multirun_psm_list({"run1": obs1, "run2": obs2})

    n = 20
    mat = np.zeros((n, 3))
    mat[:, 0] = np.tile([0.0, 1.0], n)[:n]  # noise everywhere
    mat[0:10, 1] = np.array(obs1) + 0.5  # run1 best head
    mat[10:20, 1] = 7.0  # constant for run2 -> skipped
    mat[0:10, 2] = 7.0  # constant for run1 -> skipped
    mat[10:20, 2] = np.array(obs2) + 0.5  # run2 best head
    monkeypatch.setattr(
        "ms2rescore.feature_generators.deeplc.predict", _fake_predict_factory(mat)
    )

    DeepLCFeatureGenerator(deeplc_retrain=False, processes=1).add_features(psm_list)

    for run, observed in (("run1", obs1), ("run2", obs2)):
        preds = np.array(
            [p.rescoring_features["predicted_retention_time"] for p in psm_list if p.run == run]
        )
        assert np.corrcoef(preds, np.array(observed))[0, 1] > 0.99


def test_row_alignment_preserved_across_runs(monkeypatch):
    """Each PSM must receive its own prediction after the internal qvalue sort and per-run
    matrix gather. Runs are interleaved and observed RTs are all distinct, so any cross-row
    leakage in the ``pred_matrix[run_idx, head]`` indexing would break the passthrough."""
    # Interleave two runs with distinct observed RTs; input order is not grouped by run.
    # Enough points per run for the calibration's 10% linear trail fit (needs >= 1 sample).
    interleaved = []
    for i in range(12):
        interleaved.append(("A", 10.0 + i * 5.0))
        interleaved.append(("B", 200.0 + i * 7.0))
    psms = []
    for i, (run, rt) in enumerate(interleaved):
        psm = PSM(
            peptidoform=f"{_PEPTIDES[i % len(_PEPTIDES)]}/2",
            spectrum_id=f"scan={i + 1}",
            run=run,
            is_decoy=False,
            qvalue=0.001,
            score=20.0,
            retention_time=rt,
        )
        psm.rescoring_features = {}
        psms.append(psm)
    psm_list = PSMList(psm_list=psms)

    obs = np.array([rt for _, rt in interleaved])
    mat = np.column_stack([obs, obs])  # identity: predicted head == observed
    monkeypatch.setattr(
        "ms2rescore.feature_generators.deeplc.predict", _fake_predict_factory(mat)
    )

    DeepLCFeatureGenerator(deeplc_retrain=False, processes=1).add_features(psm_list)

    for psm, rt in zip(psm_list, obs):
        assert psm.rescoring_features["observed_retention_time"] == rt
        # Identity calibration on each run -> predicted stays on this PSM's own observed RT.
        assert psm.rescoring_features["rt_diff"] < 5.0


def test_get_calibration_data_returns_indices_and_observed():
    """Contract: returns (row indices into pred_matrix, observed RTs) for the best target PSMs,
    excluding decoys, using the qvalue<=0.01 count when calibration_set_size is None."""
    import pandas as pd

    gen = DeepLCFeatureGenerator(deeplc_retrain=False)
    gen.calibration_set_size = None
    # Index labels stand in for pred_matrix row positions.
    run_df = pd.DataFrame(
        {
            "retention_time": [10.0, 20.0, 30.0, 40.0],
            "qvalue": [0.001, 0.005, 0.02, 0.001],
            "is_decoy": [False, False, False, True],
        },
        index=[5, 6, 7, 8],
    ).sort_values("qvalue")  # add_features pre-sorts by qvalue ascending

    idx, observed = gen._get_calibration_data(run_df)

    # Decoy (index 8) excluded; targets with qvalue<=0.01 are indices 5 (0.001) and 6 (0.005).
    assert list(idx) == [5, 6]
    assert list(observed) == [10.0, 20.0]


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
