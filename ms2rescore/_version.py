"""
Resolve the package version.

Prefer installed distribution metadata (importlib.metadata); if the package is not installed, fall
back to reading ``pyproject.toml`` from the repository root (PEP 621: ``[project].version``).
"""

import importlib.metadata
import json
import logging
from pathlib import Path
from typing import Dict, Optional, Tuple, Union
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

from packaging.version import Version

try:
    import tomllib as toml  # type: ignore
except ImportError:
    import tomli as toml  # type: ignore

from ms2rescore.exceptions import MS2RescoreError

LOGGER = logging.getLogger(__name__)

_GITHUB_REPO = "CompOmics/ms2rescore"
_GITHUB_TIMEOUT_SECONDS = 2.5


class UpdateCheckError(MS2RescoreError):
    """An error occurred while checking for software updates."""

    pass


def _version_from_metadata() -> Optional[Version]:
    try:
        return Version(importlib.metadata.version("ms2rescore"))
    except importlib.metadata.PackageNotFoundError:
        return None


def _version_from_pyproject() -> Optional[Version]:
    pyproject = Path(__file__).resolve().parents[1] / "pyproject.toml"
    if not pyproject.is_file():
        return None

    try:
        data = toml.loads(pyproject.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None

    project = data.get("project") if isinstance(data, dict) else None
    if isinstance(project, dict):
        ver = project.get("version")
        if isinstance(ver, str):
            return Version(ver)
    return None


def _get_latest_version(timeout_seconds: float) -> Tuple[Version, Optional[str]]:
    """Check GitHub latest release and return the version string."""
    # Prepare GitHub API request
    url = f"https://api.github.com/repos/{_GITHUB_REPO}/releases/latest"
    user_agent = "ms2rescore/ (+https://github.com/compomics/ms2rescore)"
    req = Request(url, headers={"Accept": "application/vnd.github+json", "User-Agent": user_agent})

    # Fetch and parse release info
    try:
        with urlopen(req, timeout=timeout_seconds) as resp:
            raw = resp.read()
        data = json.loads(raw.decode("utf-8", errors="replace")) if raw else {}
    except HTTPError as e:
        raise UpdateCheckError(f"HTTP {e.code}") from e
    except URLError as e:
        raise UpdateCheckError("Network unavailable or timed out") from e
    except Exception as e:
        raise UpdateCheckError("Failed to fetch release info") from e

    tag = data.get("tag_name") or data.get("name") or ""
    if not tag:
        raise UpdateCheckError("Latest release tag not found in API response")

    # Validate version string
    try:
        latest_version = Version(tag)
    except Exception as e:
        raise UpdateCheckError("Latest release tag is not a valid version") from e

    return latest_version, data.get("html_url") or data.get("url")


def get_version() -> str:
    """Return the version of this MS²Rescore installation."""
    version = _version_from_metadata() or _version_from_pyproject() or "0+unknown"
    return str(version)


def check_for_update(
    timeout_seconds: Optional[float] = None,
) -> Dict[str, Optional[Union[str, bool]]]:
    """Check GitHub latest release and report whether an update exists."""
    timeout_seconds = timeout_seconds or _GITHUB_TIMEOUT_SECONDS

    # Initialize result dictionary
    result: Dict[str, Optional[Union[str, bool]]] = {
        "update_available": False,
        "current_version": None,
        "latest_version": None,
        "html_url": None,
    }

    # Get current version
    current_version = _version_from_metadata() or _version_from_pyproject()
    if current_version is None:
        LOGGER.warning("Update check failed: Cannot determine current version")
        return result
    result["current_version"] = str(current_version)

    # Get latest version from GitHub
    try:
        latest_version, html_url = _get_latest_version(timeout_seconds)
    except UpdateCheckError as e:
        LOGGER.warning("Update check failed: %s", e)
        return result
    result["latest_version"] = str(latest_version)
    result["html_url"] = html_url

    try:
        result["update_available"] = latest_version > current_version
    except Exception:
        # If current_version can't be parsed, don't treat as updateable
        result["update_available"] = False
    return result
