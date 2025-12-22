import json
from pathlib import Path

import pytest
from psm_utils import PSM, PSMList

from ms2rescore.exceptions import MS2RescoreConfigurationError
from ms2rescore.parse_psms import parse_psms


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
            "lower_score_is_better": True,
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
            "lower_score_is_better": True,
            "psm_file": [],
            "psm_reader_kwargs": {},
            "id_decoy_pattern": None,
        }
    )

    with pytest.raises(MS2RescoreConfigurationError):
        parse_psms(config, psm_list)
