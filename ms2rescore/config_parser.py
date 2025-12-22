"""Parse configuration from command line arguments and configuration files."""

import importlib.resources
import json
import multiprocessing as mp
import re
from argparse import Namespace
from pathlib import Path
from typing import Dict, List, Union

try:
    import tomllib
except ImportError:
    import tomli as tomllib  # type: ignore

from cascade_config import CascadeConfig

from ms2rescore import package_data
from ms2rescore.exceptions import MS2RescoreConfigurationError


def _parse_output_path(configured_path, psm_file_path):
    """Parse output path and make parent dirs if required."""
    psm_file_stem = Path(psm_file_path).stem + ".ms2rescore"
    if configured_path:
        configured_path = Path(configured_path)
        # If existing dir, add psm_file stem
        if configured_path.is_dir():
            return (configured_path / psm_file_stem).as_posix()
        # If parent is existing dir, use as is (user intended as path + stem)
        elif configured_path.parent.is_dir():
            return configured_path.as_posix()
        # If none-existing dir, create dirs and add psm_file stem
        else:
            configured_path.mkdir(parents=True, exist_ok=True)
            return (configured_path / psm_file_stem).as_posix()
    else:
        # If none, use psm_file path and stem
        return (Path(psm_file_path).parent / psm_file_stem).as_posix()


def _validate_filenames(config: Dict) -> Dict:
    """Validate and infer input/output filenames."""
    # psm_file should be provided
    if not config["ms2rescore"]["psm_file"]:
        raise MS2RescoreConfigurationError("PSM file should be provided.")

    # if psm_file is a string turn into a list else leave as is
    if isinstance(config["ms2rescore"]["psm_file"], str):
        config["ms2rescore"]["psm_file"] = [config["ms2rescore"]["psm_file"]]

    # all provided psm_file(s) should exist
    psm_files = []
    for psm_file in config["ms2rescore"]["psm_file"]:
        id_file = Path(psm_file)
        if not id_file.is_file():
            raise FileNotFoundError(id_file)
        psm_files.append(id_file.as_posix())
    config["ms2rescore"]["psm_file"] = psm_files

    # spectrum_path should either be None, or existing path to file or dir
    if config["ms2rescore"]["spectrum_path"]:
        spectrum_path = Path(config["ms2rescore"]["spectrum_path"])
        if not spectrum_path.exists():
            raise FileNotFoundError(spectrum_path)
        config["ms2rescore"]["spectrum_path"] = spectrum_path.as_posix()

    # Parse output_path
    config["ms2rescore"]["output_path"] = _parse_output_path(
        config["ms2rescore"]["output_path"], config["ms2rescore"]["psm_file"][0]
    )

    # Parse config_file as posix path to avoid combination of forward and backward slashes
    if config["ms2rescore"]["config_file"]:
        config["ms2rescore"]["config_file"] = Path(config["ms2rescore"]["config_file"]).as_posix()

    return config


def _validate_processes(config: Dict) -> Dict:
    """Validate requested processes with available cpu count."""
    n_available = mp.cpu_count()
    if (config["ms2rescore"]["processes"] == -1) or (
        config["ms2rescore"]["processes"] > n_available
    ):
        config["ms2rescore"]["processes"] = n_available
    return config


def _validate_regular_expressions(config: Dict) -> Dict:
    """Validate regular expressions in configuration."""
    for field in [
        "psm_id_pattern",
        "spectrum_id_pattern",
        "psm_id_rt_pattern",
        "psm_id_im_pattern",
    ]:
        if config["ms2rescore"][field]:

            # Check if valid regex
            try:
                pattern = re.compile(config["ms2rescore"][field])
            except re.error as e:
                raise MS2RescoreConfigurationError(
                    f"Invalid regular expression provided for '{field}': {e}"
                ) from e

            # Check if regex has exactly one capturing group
            if pattern.groups != 1:
                raise MS2RescoreConfigurationError(
                    f"Regular expression for '{field}' should contain exactly one "
                    "capturing group. Please check and try again. "
                    "See https://ms2rescore.readthedocs.io/en/stable/userguide/configuration/#mapping-psms-to-spectra "
                    "for more information."
                )

    return config


def parse_configurations(configurations: List[Union[dict, str, Path, Namespace]]) -> Dict:
    """
    Parse and validate MS²Rescore configuration files and CLI arguments.

    Default configuration, user configuration files, and CLI/class arguments are parsed
    in cascading order, with each successive configuration taking priority over the
    previous.

    Parameters
    ----------
    configurations: Dict, str, Path, Namespace, List[Dict, str, Path, Namespace]
        configuration dictionary, path to configuration files, argparse Namespace, or a list of the
        above.
    """
    if not isinstance(configurations, list):
        configurations = [configurations]

    # Initialize CascadeConfig with validation schema and defaults
    config_schema = importlib.resources.open_text(package_data, "config_schema.json")
    config_default = importlib.resources.open_text(package_data, "config_default.json")
    cascade_conf = CascadeConfig(
        validation_schema=json.load(config_schema),
        none_overrides_value=False,
        max_recursion_depth=1,
    )
    cascade_conf.add_dict(json.load(config_default))

    # Add configurations
    for config in configurations:
        if not config:
            continue
        if isinstance(config, dict):
            cascade_conf.add_dict(config)
        elif isinstance(config, str) or isinstance(config, Path):
            if Path(config).suffix.lower() == ".json":
                cascade_conf.add_json(config)
            elif Path(config).suffix.lower() == ".toml":
                cascade_conf.add_dict(dict(tomllib.load(Path(config).open("rb"))))
            else:
                raise MS2RescoreConfigurationError(
                    "Unknown file extension for configuration file. Should be `json` or " "`toml`."
                )
        elif isinstance(config, Namespace):
            cascade_conf.add_namespace(config, subkey="ms2rescore")
        else:
            raise ValueError(
                "Configuration should be a dictionary, argparse Namespace, or path to a "
                "configuration file."
            )

    # Parse configurations
    config = cascade_conf.parse()

    # Validate and infer filenames and number of parallel processes
    config = _validate_filenames(config)
    config = _validate_processes(config)
    config = _validate_regular_expressions(config)

    # Convert feature_generators and rescoring_engine names to lowercase
    config["ms2rescore"]["feature_generators"] = {
        k.lower(): v for k, v in config["ms2rescore"]["feature_generators"].items()
    }
    config["ms2rescore"]["rescoring_engine"] = {
        k.lower(): v for k, v in config["ms2rescore"]["rescoring_engine"].items()
    }

    return config
