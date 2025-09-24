import logging
import os
from glob import glob
from pathlib import Path
from typing import Optional, Union
import json
import re
import socket
from typing import Dict
from urllib.request import Request, urlopen
from urllib.error import URLError, HTTPError

from ms2rescore.exceptions import MS2RescoreConfigurationError
from ms2rescore_rs import is_supported_file_type

logger = logging.getLogger(__name__)


def infer_spectrum_path(
    configured_path: Union[str, Path, None],
    run_name: Optional[str] = None,
) -> Union[str, Path]:
    """
    Infer spectrum path from passed path and expected filename (e.g. from PSM file).

    Parameters
    ----------
    configured_path: str, Path, None
        User-defined path to spectrum file or directory containing spectrum file
    run_name : str, optional
        MS run name (stem of spectrum filename), e.g., as expected from PSM file.

    """
    # If no spectrum path configured, use expected run_name in default dir
    if not configured_path:
        if run_name:
            resolved_path = os.path.join(".", run_name)
        else:
            raise MS2RescoreConfigurationError(
                "Could not resolve spectrum file name: No spectrum path configured "
                "and no run name in PSM file found."
            )

    else:
        is_bruker_dir = configured_path.endswith(".d") or _is_minitdf(configured_path)

        # If passed path is directory (that is not Bruker raw), join with run name
        if os.path.isdir(configured_path) and not is_bruker_dir:
            if run_name:
                resolved_path = os.path.join(configured_path, run_name)
            else:
                raise MS2RescoreConfigurationError(
                    "Could not resolve spectrum file name: Spectrum path is directory "
                    "but no run name in PSM file found."
                )

        # If passed path is file, use that, but warn if basename doesn't match expected
        elif os.path.isfile(configured_path) or (os.path.isdir(configured_path) and is_bruker_dir):
            if run_name and Path(configured_path).stem != Path(run_name).stem:
                logger.warning(
                    "Passed spectrum path (`%s`) does not match run name found in PSM "
                    "file (`%s`). Continuing with passed spectrum path.",
                    configured_path,
                    run_name,
                )
            resolved_path = configured_path
        else:
            raise MS2RescoreConfigurationError(
                "Configured `spectrum_path` must be `None` or a path to an existing file "
                "or directory. If `None` or path to directory, spectrum run information "
                "should be present in the PSM file."
            )

    # Match with file extension if not in resolved_path yet
    if not is_supported_file_type(resolved_path) or not os.path.exists(resolved_path):
        for filename in glob(resolved_path + "*"):
            if is_supported_file_type(filename):
                resolved_path = filename
                break
        else:
            raise MS2RescoreConfigurationError(
                f"Resolved spectrum filename ('{resolved_path}') does not contain a supported "
                "file extension (mzML, MGF, or .d) and could not find any matching existing "
                "files."
            )

    return Path(resolved_path)


def _is_minitdf(spectrum_file: str) -> bool:
    """
    Check if the spectrum file is a Bruker miniTDF folder.

    A Bruker miniTDF folder has no fixed name, but contains files matching the patterns
    ``*ms2spectrum.bin`` and ``*ms2spectrum.parquet``.
    """
    files = set(Path(spectrum_file).glob("*ms2spectrum.bin"))
    files.update(Path(spectrum_file).glob("*ms2spectrum.parquet"))
    return len(files) >= 2


def _strip_v(s: str) -> str:
    return s.lstrip("vV").strip() if isinstance(s, str) else s


def _normalize_for_tuple(v: str):
    """
    Very small-footprint normalizer that yields a tuple suitable for comparison
    when 'packaging' is unavailable. It prefers numeric precedence and
    compares leftover text segments lexically. Pre-releases are not modeled,
    but common tags like 'v1.2.3' vs '1.10.0' work as expected.
    """
    v = _strip_v(v)
    # Split into alternating number / text chunks: "1.2.3-rc1" -> ['1','.', '2','.', '3','-','rc','1']
    chunks = re.findall(r"\d+|[A-Za-z]+|[^A-Za-z0-9]+", v)
    norm = []
    for c in chunks:
        if c.isdigit():
            norm.append((0, int(c)))
        else:
            # collapse punctuation to a single separator to keep ordering stable
            txt = c if re.search(r"[A-Za-z]", c) else "."
            norm.append((1, txt))
    return tuple(norm)


def _compare_versions(a: str, b: str) -> int:
    """
    Returns -1 if a<b, 0 if a==b, 1 if a>b.
    Tries packaging.version.Version if available; otherwise falls back to tuple comparison.
    """
    a = _strip_v(a)
    b = _strip_v(b)
    try:
        from packaging.version import Version  # optional

        va, vb = Version(a), Version(b)
        return (va > vb) - (va < vb)
    except Exception:
        ta, tb = _normalize_for_tuple(a), _normalize_for_tuple(b)
        return (ta > tb) - (ta < tb)


def check_for_update(
    current_version: str,
    repo: str = "CompOmics/ms2rescore",
    timeout_seconds: float = 2.5,
    user_agent: str = "ms2rescore-update-checker",
) -> Dict[str, Optional[str]]:
    """
    Query GitHub for the latest release of `repo` and compare to `current_version`.

    Returns a dict:
      {
        "ok": bool,                    # True if HTTP+parse succeeded
        "is_update": bool,             # True if a newer version exists
        "current_version": str,        # Echo of your input
        "latest_version": Optional[str],# Latest tag (without leading 'v')
        "html_url": Optional[str],     # Release page
        "error": Optional[str],        # Non-fatal reason in case of failure
      }

    Notes
    -----
    - Never raises: all exceptions are caught.
    - Safe offline: on network issues, returns ok=False, is_update=False.
    - Works in PyInstaller (no extra deps). If 'packaging' is present, it will use it.
    """
    result = {
        "ok": False,
        "is_update": False,
        "current_version": current_version,
        "latest_version": None,
        "html_url": None,
        "error": None,
    }

    url = f"https://api.github.com/repos/{repo}/releases/latest"
    req = Request(
        url,
        headers={
            "Accept": "application/vnd.github+json",
            "User-Agent": user_agent,
        },
        method="GET",
    )

    try:
        with urlopen(req, timeout=timeout_seconds) as resp:
            raw = resp.read()
        data = json.loads(raw.decode("utf-8", errors="replace")) if raw else {}

        tag = data.get("tag_name") or data.get("name") or ""
        latest = _strip_v(tag)
        html_url = data.get("html_url") or data.get("url")

        if not latest:
            result["error"] = "Latest release tag not found in API response."
            return result

        result["latest_version"] = latest
        result["html_url"] = html_url
        result["ok"] = True

        cmp = _compare_versions(latest, current_version)
        result["is_update"] = cmp == 1
        return result

    except HTTPError as e:
        result["error"] = f"HTTPError {e.code}"
        return result
    except (URLError, socket.timeout):
        result["error"] = "Network unavailable or timed out"
        return result
    except (json.JSONDecodeError, ValueError):
        result["error"] = "Failed to parse API response"
        return result
    except Exception:
        # Final safety net; do not raise
        result["error"] = "Unexpected error"
        return result
