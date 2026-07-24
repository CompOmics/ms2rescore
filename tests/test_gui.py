"""Tests for ms2rescore.gui.app."""

import logging
from types import SimpleNamespace

import pytest

pytest.importorskip("tkinter")

from ms2rescore.gui.app import DeepLCConfiguration, MS2PIPConfiguration, _setup_logging  # noqa: E402
from ms2rescore.gui.function2ctk import _apply_selected_log_level  # noqa: E402


def test_setup_logging_writes_txt_and_html_log(tmp_path):
    """_setup_logging must produce both the plain-text log and an HTML log via the returned
    console -- the GUI run used to only ever write the text log, silently skipping the HTML
    log that the CLI and documentation both promise.
    """
    root_logger = logging.getLogger()
    original_handlers = root_logger.handlers[:]
    try:
        log_txt = str(tmp_path / "test.log.txt")
        console = _setup_logging("info", log_txt)
        logging.getLogger("test_gui").info("hello from gui logging test")

        log_html = str(tmp_path / "test.log.html")
        console.save_html(log_html)

        assert (tmp_path / "test.log.txt").is_file()
        assert (tmp_path / "test.log.html").is_file()
        assert "hello from gui logging test" in (tmp_path / "test.log.html").read_text(
            encoding="utf-8"
        )
    finally:
        for handler in root_logger.handlers[:]:
            if handler not in original_handlers:
                root_logger.removeHandler(handler)


def test_apply_selected_log_level_overrides_config_copy():
    config = {"ms2rescore": {"log_level": "info", "psm_file": ["input.tsv"]}}

    updated = _apply_selected_log_level(config, "debug")

    assert updated["ms2rescore"]["log_level"] == "debug"
    assert config["ms2rescore"]["log_level"] == "info"


def test_ms2pip_configuration_get_returns_annotation_settings():
    config = MS2PIPConfiguration.__new__(MS2PIPConfiguration)
    config.enabled = SimpleNamespace(get=lambda: True)
    config.model = SimpleNamespace(get=lambda: "HCD2021")
    config.fragmentation_model = SimpleNamespace(get=lambda: "etd")
    config.ms2_tolerance = SimpleNamespace(get=lambda: 0.05)
    config.tolerance_mode = SimpleNamespace(get=lambda: "ppm")

    enabled, feature_config, annotation_config = MS2PIPConfiguration.get(config)

    assert enabled is True
    assert feature_config == {"model": "HCD2021"}
    assert annotation_config == {
        "fragmentation_model": "etd",
        "tolerance_value": 0.05,
        "tolerance_mode": "ppm",
    }


def test_deeplc_configuration_get_uses_epochs_key():
    config = DeepLCConfiguration.__new__(DeepLCConfiguration)
    config.enabled = SimpleNamespace(get=lambda: True)
    config.transfer_learning = SimpleNamespace(get=lambda: True)
    config.num_epochs = SimpleNamespace(get=lambda: 30)
    config.calibration_set_size = SimpleNamespace(get=lambda: "0.25")

    enabled, deeplc_config = DeepLCConfiguration.get(config)

    assert enabled is True
    assert deeplc_config["finetune"] is True
    assert deeplc_config["epochs"] == 30
    assert "n_epochs" not in deeplc_config
