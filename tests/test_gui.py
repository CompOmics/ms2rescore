"""Tests for ms2rescore.gui.app."""

import logging

from ms2rescore.gui.app import _setup_logging


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
