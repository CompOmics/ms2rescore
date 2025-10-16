"""Modular and user-friendly platform for AI-assisted rescoring of peptide identifications ."""

# __version__ is provided by the small helper module _version which prefers
# installed distribution metadata and falls back to reading pyproject.toml.
# The single source of truth for the version is pyproject.toml.

__all__ = [
    "parse_configurations",
    "rescore",
]

from warnings import filterwarnings

# mzmlb is not used, so hdf5plugin is not needed
filterwarnings(
    "ignore",
    message="hdf5plugin is missing",
    category=UserWarning,
    module="psims.mzmlb",
)

# Attempt to suppress pyopenms warning about OPENMS_DATA_PATH already set
filterwarnings(
    "ignore",
    message="Warning: OPENMS_DATA_PATH environment variable already exists.*",
    category=UserWarning,
    module="pyopenms",
)

from ms2rescore._version import get_version  # noqa: E402
from ms2rescore.config_parser import parse_configurations  # noqa: E402
from ms2rescore.core import rescore  # noqa: E402

__version__ = get_version()
