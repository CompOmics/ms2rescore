"""MS²Rescore: Sensitive PSM rescoring with predicted MS² peak intensities and RTs."""

import argparse
import cProfile
import importlib.resources
import json
import logging
import sys
from pathlib import Path
from typing import Union

from rich.console import Console
from rich.logging import RichHandler
from rich.text import Text

from ms2rescore import __version__, package_data
from ms2rescore.config_parser import parse_configurations
from ms2rescore.core import rescore
from ms2rescore.exceptions import MS2RescoreConfigurationError
from ms2rescore.utils import check_for_update

try:
    import matplotlib.pyplot as plt

    plt.set_loglevel("warning")
except ImportError:
    pass

LOG_MAPPING = {
    "critical": logging.CRITICAL,
    "error": logging.ERROR,
    "warning": logging.WARNING,
    "info": logging.INFO,
    "debug": logging.DEBUG,
}
LOGGER = logging.getLogger(__name__)
CONSOLE = Console(record=True)


def _print_credits(tims=False):
    """Print software credits to terminal."""
    text = Text()
    text.append("\n")
    if tims:
        text.append("TIMS²Rescore", style="bold link https://github.com/compomics/tims2rescore")
    else:
        text.append("MS²Rescore", style="bold link https://github.com/compomics/ms2rescore")
    text.append(f" (v{__version__})\n", style="bold")
    if tims:
        text.append("MS²Rescore tuned for timsTOF DDA-PASEF data.\n", style="italic")
    text.append("Developed at CompOmics, VIB / Ghent University, Belgium.\n")
    text.append("Please cite: ")
    if tims:
        text.append(
            "Declercq & Devreese et al. bioRxiv (2024)",
            style="link https://doi.org/10.1101/2024.05.29.596400",
        )
    else:
        text.append(
            "Buur & Declercq et al. JPR (2024)",
            style="link https://doi.org/10.1021/acs.jproteome.3c00785",
        )
    text.append("\n")
    if tims:
        text.stylize("#006cb5")
    CONSOLE.print(text)


def _argument_parser() -> argparse.ArgumentParser:
    """Parse CLI arguments."""
    parser = argparse.ArgumentParser(
        description="MS²Rescore: Sensitive PSM rescoring with predicted features.",
        formatter_class=lambda prog: argparse.HelpFormatter(prog, max_help_position=42),
    )
    parser.add_argument("-v", "--version", action="version", version=__version__)
    parser.add_argument(
        "-p",
        "--psm-file",
        metavar="FILE",
        action="store",
        type=str,
        nargs="*",
        dest="psm_file",
        help="path to PSM file (PIN, mzIdentML, MaxQuant msms, X!Tandem XML...)",
    )
    parser.add_argument(
        "-t",
        "--psm-file-type",
        metavar="STR",
        action="store",
        type=str,
        dest="psm_file_type",
        help="PSM file type (default: 'infer')",
    )
    parser.add_argument(
        "-s",
        "--spectrum-path",
        metavar="FILE/DIR",
        action="store",
        type=str,
        dest="spectrum_path",
        help="path to MGF/mzML spectrum file or directory with spectrum files (default: derived\
            from identification file)",
    )
    parser.add_argument(
        "-c",
        "--config-file",
        metavar="FILE",
        action="store",
        type=str,
        dest="config_file",
        help="path to MS²Rescore configuration file (see README.md)",
    )
    parser.add_argument(
        "-o",
        "--output-path",
        metavar="FILE",
        action="store",
        type=str,
        dest="output_path",
        help="Path and stem for output file names (default: derive from identification file)",
    )
    parser.add_argument(
        "-l",
        "--log-level",
        metavar="STR",
        action="store",
        type=str,
        dest="log_level",
        help="logging level (default: `info`)",
    )
    parser.add_argument(
        "-n",
        "--processes",
        metavar="INT",
        action="store",
        type=int,
        dest="processes",
        help="number of parallel processes available to MS²Rescore",
    )
    parser.add_argument(
        "-f",
        "--fasta-file",
        metavar="FILE",
        action="store",
        type=str,
        dest="fasta_file",
        help="path to FASTA file",
    )
    parser.add_argument(
        "--write-report",
        action="store_true",
        default=None,
        dest="write_report",
        help="boolean whether to write an HTML report (default: True)",
    )
    parser.add_argument(
        "--profile",
        action="store_true",
        default=None,
        dest="profile",
        help="boolean to enable profiling with cProfile",
    )

    return parser


def _setup_logging(passed_level: str, log_file: Union[str, Path]):
    """Setup logging for writing to log file and Rich Console."""
    if passed_level not in LOG_MAPPING:
        raise MS2RescoreConfigurationError(
            f"Invalid log level '{passed_level}'. "
            f"Valid levels are: {', '.join(LOG_MAPPING.keys())}"
        )
    logging.basicConfig(
        format="%(name)s // %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        level=LOG_MAPPING[passed_level],
        handlers=[
            logging.FileHandler(log_file, mode="w", encoding="utf-8"),
            RichHandler(rich_tracebacks=True, console=CONSOLE, show_path=False),
        ],
    )


def profile(fnc, filepath):
    """A decorator that uses cProfile to profile a function"""

    def inner(*args, **kwargs):
        with cProfile.Profile() as profiler:
            return_value = fnc(*args, **kwargs)
        profiler.dump_stats(filepath + ".profile.prof")
        return return_value

    return inner


def main_tims():
    """Run MS²Rescore command-line interface in TIMS²Rescore mode."""
    main(tims=True)


def main(tims=False):
    """Run MS²Rescore command-line interface."""
    _print_credits(tims)

    # Parse CLI arguments and configuration file
    parser = _argument_parser()
    cli_args = parser.parse_args()

    configurations = []
    if tims:
        configurations.append(
            json.load(importlib.resources.open_text(package_data, "config_default_tims.json"))
        )
    if cli_args.config_file:
        configurations.append(cli_args.config_file)
    configurations.append(cli_args)

    try:
        config = parse_configurations(configurations)
    except MS2RescoreConfigurationError as e:
        LOGGER.critical(e)
        sys.exit(1)

    # Setup logging
    _setup_logging(
        config["ms2rescore"]["log_level"], config["ms2rescore"]["output_path"] + ".log.txt"
    )

    # check for software updates
    info = check_for_update(__version__, repo="CompOmics/ms2rescore")
    if info["is_update"]:
        LOGGER.info(
            f"New version of MS²Rescore available: {info['latest_version']} "
            f"(you are using {__version__})"
        )
        LOGGER.info(f"Download the latest version at: {info['html_url']}")
    else:
        LOGGER.debug("You are using the latest version of MS²Rescore.")

    # Run MS²Rescore
    try:
        if cli_args.profile:
            profiled_rescore = profile(rescore, config["ms2rescore"]["output_path"])
            profiled_rescore(configuration=config)
        else:
            rescore(configuration=config)
    except Exception as e:
        LOGGER.exception(e)
        sys.exit(1)
    finally:
        CONSOLE.save_html(config["ms2rescore"]["output_path"] + ".log.html")


if __name__ == "__main__":
    main()
