import importlib
import json
import logging

from packaging.version import Version

import ms2rescore._version


def test_version_matches_pyproject_or_is_nonempty():
    ms2rescore = importlib.import_module("ms2rescore")
    pkg_ver = getattr(ms2rescore, "__version__", None)
    assert isinstance(pkg_ver, str) and pkg_ver, "__version__ must be a non-empty string"

    # The package version should match the resolver used by the package itself
    # (installed metadata or pyproject). Use get_version() to get the authoritative
    # value in the current environment.
    assert pkg_ver == ms2rescore._version.get_version()


class _FakeResp:
    def __init__(self, data: bytes):
        self._data = data

    def read(self) -> bytes:
        return self._data

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        return False


def test_get_version_nonempty():
    ver = ms2rescore._version.get_version()
    assert isinstance(ver, str) and ver, "get_version() should return a non-empty string"


def test_check_for_update_network_error(monkeypatch):
    # Simulate a network error by making urlopen raise URLError
    from urllib.error import URLError

    def _raise(*args, **kwargs):
        raise URLError("no network")

    monkeypatch.setattr(ms2rescore._version, "urlopen", _raise)

    info = ms2rescore._version.check_for_update(timeout_seconds=0.01)
    assert isinstance(info, dict)
    # New API: check_for_update returns minimal keys. Update should not be reported
    assert info.get("update_available") is False
    assert "current_version" in info
    assert info.get("latest_version") is None


def test_check_for_update_detects_update(monkeypatch):
    # Mock the release API to return a very large version tag
    payload = {"tag_name": "v999.0.0", "html_url": "https://example.org/release/999"}
    raw = json.dumps(payload).encode("utf-8")

    def fake_urlopen(req, timeout=None):
        # Ensure the User-Agent header contains the runtime-resolved version
        ua = None
        try:
            ua = req.get_header("User-agent") or req.get_header("User-Agent")
        except Exception:
            ua = getattr(req, "headers", {}).get("User-Agent")
        assert ua is not None and "ms2rescore/" in ua
        return _FakeResp(raw)

    # Make the module helpers report an old version so the check detects an update
    monkeypatch.setattr(ms2rescore._version, "urlopen", fake_urlopen)
    monkeypatch.setattr(ms2rescore._version, "_version_from_metadata", lambda: Version("0.0.1"))
    monkeypatch.setattr(ms2rescore._version, "_version_from_pyproject", lambda: None)

    info = ms2rescore._version.check_for_update(timeout_seconds=1)
    assert isinstance(info, dict)
    assert info.get("update_available") is True
    assert info.get("latest_version") == "999.0.0"
    assert info.get("html_url") == "https://example.org/release/999"


def test_check_for_update_no_current_version(monkeypatch, caplog):
    # Simulate missing metadata and pyproject
    monkeypatch.setattr(ms2rescore._version, "_version_from_metadata", lambda: None)
    monkeypatch.setattr(ms2rescore._version, "_version_from_pyproject", lambda: None)
    with caplog.at_level(logging.WARNING):
        info = ms2rescore._version.check_for_update(timeout_seconds=0.01)
    assert isinstance(info, dict)
    assert info.get("current_version") is None
    assert info.get("update_available") is False
    assert "Update check failed" in caplog.text


def test_check_for_update_missing_tag(monkeypatch):
    # API returns JSON without tag_name/name
    payload = {}
    raw = json.dumps(payload).encode("utf-8")

    def fake_urlopen(req, timeout=None):
        return _FakeResp(raw)

    monkeypatch.setattr(ms2rescore._version, "urlopen", fake_urlopen)
    monkeypatch.setattr(ms2rescore._version, "_version_from_metadata", lambda: Version("1.0.0"))
    monkeypatch.setattr(ms2rescore._version, "_version_from_pyproject", lambda: None)

    info = ms2rescore._version.check_for_update(timeout_seconds=1)
    assert isinstance(info, dict)
    assert info.get("update_available") is False
    assert info.get("latest_version") is None


def test_check_for_update_handles_updatecheckerror(monkeypatch, caplog):
    # Simulate _get_latest_version raising UpdateCheckError
    def raise_err(timeout):
        raise ms2rescore._version.UpdateCheckError("HTTP 500")

    monkeypatch.setattr(ms2rescore._version, "_get_latest_version", raise_err)
    monkeypatch.setattr(ms2rescore._version, "_version_from_metadata", lambda: Version("1.0.0"))
    monkeypatch.setattr(ms2rescore._version, "_version_from_pyproject", lambda: None)

    with caplog.at_level(logging.WARNING):
        info = ms2rescore._version.check_for_update(timeout_seconds=0.01)
    assert isinstance(info, dict)
    assert info.get("update_available") is False
    assert "Update check failed" in caplog.text
