"""Modular and user-friendly platform for AI-assisted rescoring of peptide identifications ."""

__version__ = "3.2.0.dev2"
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

from ms2rescore.config_parser import parse_configurations  # noqa: E402
from ms2rescore.core import rescore  # noqa: E402
