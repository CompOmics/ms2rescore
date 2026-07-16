import json
from pathlib import Path

import numpy as np
import pytest
from psm_utils import PSM, PSMList

from ms2rescore.exceptions import MS2RescoreConfigurationError
from ms2rescore.parse_psms import infer_score_direction, parse_psms

_RESIDUES = "ACDEFGHIKLMNPQRSTVWY"


def _peptide_seq(i: int) -> str:
    """Build a unique, valid amino-acid sequence for index `i`."""
    a = _RESIDUES[i % len(_RESIDUES)]
    b = _RESIDUES[(i // len(_RESIDUES)) % len(_RESIDUES)]
    return f"PEPT{a}{b}IDEK"


def _make_direction_psm_list(n_spectra=30, seed=0, lower_is_better=False, run="run1"):
    """Separable synthetic PSMList for score-direction inference tests."""
    rng = np.random.default_rng(seed)
    psms = []
    for i in range(n_spectra):
        is_decoy = rng.random() < 0.3
        shift = 0.0 if is_decoy else 8.0
        score = float(rng.normal(0, 1) + shift)
        if lower_is_better:
            score = -score
        psms.append(
            PSM(
                peptidoform=f"{_peptide_seq(i)}/2",
                spectrum_id=str(i),
                run=run,
                is_decoy=is_decoy,
                score=score,
            )
        )
    return PSMList(psm_list=psms)


@pytest.fixture(scope="module")
def default_config():
    cfg_path = Path(__file__).parents[1] / "ms2rescore" / "package_data" / "config_default.json"
    cfg = json.loads(cfg_path.read_text())["ms2rescore"]
    return cfg


@pytest.fixture
def psm_list_factory():
    def _factory(ids):
        return PSMList(
            psm_list=[
                PSM(
                    peptidoform="PEPTIDE/2",
                    run="run1",
                    spectrum_id=sid,
                    retention_time=None,
                    ion_mobility=None,
                    precursor_mz=None,
                )
                for sid in ids
            ]
        )

    return _factory


def test_psm_id_pattern_success(default_config, psm_list_factory):
    psm_list = psm_list_factory(["scan:1:fileA", "scan:2:fileA"])
    # Ensure at least one decoy is present so parse_psms does not raise
    psm_list[0].is_decoy = True
    config = dict(default_config)
    config.update(
        {
            "psm_id_pattern": r"scan:(\d+):.*",
            "psm_file": [],
            "psm_reader_kwargs": {},
            "id_decoy_pattern": None,
        }
    )

    result = parse_psms(config, psm_list)
    assert list(result["spectrum_id"]) == ["1", "2"]


def test_psm_id_pattern_collapses_unique_ids(default_config, psm_list_factory):
    psm_list = psm_list_factory(["scan:1:fileA", "scan:1:fileB"])
    config = dict(default_config)
    config.update(
        {
            "psm_id_pattern": r"scan:(\d+):.*",
            "psm_file": [],
            "psm_reader_kwargs": {},
            "id_decoy_pattern": None,
        }
    )

    with pytest.raises(MS2RescoreConfigurationError):
        parse_psms(config, psm_list)


def test_infer_score_direction_detects_higher_is_better():
    psm_list = _make_direction_psm_list(n_spectra=30, seed=1, lower_is_better=False)
    assert infer_score_direction(psm_list, train_fdr=0.1) is False


def test_infer_score_direction_detects_lower_is_better():
    psm_list = _make_direction_psm_list(n_spectra=30, seed=2, lower_is_better=True)
    assert infer_score_direction(psm_list, train_fdr=0.1) is True


def test_infer_score_direction_defaults_to_higher_when_too_few_psms():
    psm_list = PSMList(
        psm_list=[
            PSM(peptidoform="PEPTIDEK/2", spectrum_id="1", run="run1", is_decoy=False, score=1.0)
        ]
    )
    assert infer_score_direction(psm_list) is False


def test_infer_score_direction_defaults_to_higher_when_no_finite_scores():
    psms = [
        PSM(
            peptidoform="PEPTIDEK/2",
            spectrum_id="1",
            run="run1",
            is_decoy=False,
            score=float("nan"),
        ),
        PSM(
            peptidoform="PEPTIDER/2",
            spectrum_id="2",
            run="run1",
            is_decoy=True,
            score=float("nan"),
        ),
    ]
    assert infer_score_direction(PSMList(psm_list=psms)) is False


def test_infer_score_direction_groups_by_run(monkeypatch):
    """Multiple input files can share native spectrum/scan IDs (e.g. per-run MaxQuant/MSGFPlus
    output both starting at scan 1) -- ristretto's competition must be grouped by run_col, or
    unrelated PSMs from different runs that happen to share a spectrum_id would be treated as
    competing for the same spectrum.
    """
    import ristretto

    from ms2rescore import parse_psms as parse_psms_module

    psm_list_a = _make_direction_psm_list(n_spectra=10, seed=5, run="runA")
    psm_list_b = _make_direction_psm_list(n_spectra=10, seed=6, run="runB")
    psm_list = PSMList(psm_list=list(psm_list_a) + list(psm_list_b))

    calls = []
    original_evaluate = ristretto.evaluate

    def _spy(*args, **kwargs):
        calls.append(kwargs.get("run_col"))
        return original_evaluate(*args, **kwargs)

    monkeypatch.setattr(parse_psms_module.ristretto, "evaluate", _spy)
    infer_score_direction(psm_list, train_fdr=0.1)

    assert calls
    assert all(run_col == "run" for run_col in calls)
