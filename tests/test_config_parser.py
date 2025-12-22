from pathlib import Path

import pytest

from ms2rescore.config_parser import _parse_output_path, _validate_regular_expressions
from ms2rescore.exceptions import MS2RescoreConfigurationError


def test__parse_output_path():
    # Ensure that test dir exists
    Path("examples/id").mkdir(parents=True, exist_ok=True)
    test_psm_file = "some/dir/psm_file.mzid"

    test_cases = [
        ("examples/id", "examples/id/psm_file.ms2rescore"),  # Existing dir
        ("examples/id/custom_stem", "examples/id/custom_stem"),  # Parent is existing dir
        ("some/other_dir", "some/other_dir/psm_file.ms2rescore"),  # None-existing dir
        (
            "some/other_dir/",
            "some/other_dir/psm_file.ms2rescore",
        ),  # None-existing dir, with trailing slash
        (None, "some/dir/psm_file.ms2rescore"),
    ]

    for output_path, expected in test_cases:
        assert _parse_output_path(output_path, test_psm_file) == expected


def test__validate_regular_expressions_valid():
    """Test _validate_regular_expressions with valid regex patterns."""
    # Valid pattern with one capturing group
    config = {
        "ms2rescore": {
            "psm_id_pattern": r"scan:(\d+):.*",
            "spectrum_id_pattern": r"spectrum_(\d+)",
            "psm_id_rt_pattern": None,
            "psm_id_im_pattern": None,
        }
    }
    result = _validate_regular_expressions(config)
    assert result == config


def test__validate_regular_expressions_none():
    """Test _validate_regular_expressions with None patterns."""
    config = {
        "ms2rescore": {
            "psm_id_pattern": None,
            "spectrum_id_pattern": None,
            "psm_id_rt_pattern": None,
            "psm_id_im_pattern": None,
        }
    }
    result = _validate_regular_expressions(config)
    assert result == config


def test__validate_regular_expressions_invalid_regex():
    """Test _validate_regular_expressions with invalid regex syntax."""
    config = {
        "ms2rescore": {
            "psm_id_pattern": r"scan:(\d+",  # Missing closing parenthesis
            "spectrum_id_pattern": None,
            "psm_id_rt_pattern": None,
            "psm_id_im_pattern": None,
        }
    }
    with pytest.raises(MS2RescoreConfigurationError, match="Invalid regular expression"):
        _validate_regular_expressions(config)


def test__validate_regular_expressions_no_capturing_group():
    """Test _validate_regular_expressions with no capturing groups."""
    config = {
        "ms2rescore": {
            "psm_id_pattern": r"scan:\d+:.*",  # No capturing group
            "spectrum_id_pattern": None,
            "psm_id_rt_pattern": None,
            "psm_id_im_pattern": None,
        }
    }
    with pytest.raises(
        MS2RescoreConfigurationError, match="should contain exactly one capturing group"
    ):
        _validate_regular_expressions(config)


def test__validate_regular_expressions_multiple_capturing_groups():
    """Test _validate_regular_expressions with multiple capturing groups."""
    config = {
        "ms2rescore": {
            "psm_id_pattern": r"scan:(\d+):(.*)",  # Two capturing groups
            "spectrum_id_pattern": None,
            "psm_id_rt_pattern": None,
            "psm_id_im_pattern": None,
        }
    }
    with pytest.raises(
        MS2RescoreConfigurationError, match="should contain exactly one capturing group"
    ):
        _validate_regular_expressions(config)


def test__validate_regular_expressions_spectrum_id_pattern_invalid():
    """Test _validate_regular_expressions with invalid spectrum_id_pattern."""
    config = {
        "ms2rescore": {
            "psm_id_pattern": None,
            "spectrum_id_pattern": r"spectrum_\d+",  # No capturing group
            "psm_id_rt_pattern": None,
            "psm_id_im_pattern": None,
        }
    }
    with pytest.raises(
        MS2RescoreConfigurationError, match="should contain exactly one capturing group"
    ):
        _validate_regular_expressions(config)
