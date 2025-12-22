from unittest.mock import MagicMock, patch

import numpy as np
import pytest
from psm_utils import PSM, PSMList

from ms2rescore.exceptions import MS2RescoreConfigurationError
from ms2rescore.parse_spectra import (
    MSDataType,
    SpectrumParsingError,
    _get_precursor_values,
    add_precursor_values,
)


def psm_list_factory(ids):
    return PSMList(
        psm_list=[
            PSM(
                peptidoform="PEPTIDE/2",
                run="run1",
                spectrum_id=sid,
                retention_time=None,
                ion_mobility=None,
                precursor_mz=None,
            )
            for sid in ids
        ]
    )


@pytest.fixture
def mock_psm_list():
    return PSMList(
        psm_list=[
            PSM(
                peptidoform="PEPTIDE/2",
                run="run1",
                spectrum_id="spectrum1",
                retention_time=None,
                ion_mobility=None,
                precursor_mz=None,
            ),
            PSM(
                peptidoform="PEPTIDE/2",
                run="run1",
                spectrum_id="spectrum2",
                retention_time=None,
                ion_mobility=None,
                precursor_mz=None,
            ),
        ]
    )


@pytest.fixture
def mock_precursor_info():
    return {
        "spectrum1": MagicMock(mz=529.7935187324, rt=10.5, im=1.0),
        "spectrum2": MagicMock(mz=651.83, rt=12.3, im=1.2),
    }


@pytest.fixture
def mock_precursor_info_missing_im():
    return {
        "spectrum1": MagicMock(mz=529.7935187324, rt=10.5, im=0.0),
        "spectrum2": MagicMock(mz=651.83, rt=12.3, im=0.0),
    }


@pytest.fixture
def mock_precursor_info_incomplete():
    return {
        "spectrum1": MagicMock(mz=529.7935187324, rt=10.5, im=1.0),
        # spectrum2 intentionally missing
    }


def test_spectrum_id_pattern_nonmatching(monkeypatch):
    """If the provided spectrum_id_pattern doesn't match any spectrum-file IDs, raise."""
    psm_list = psm_list_factory(["scan:1:fileA"])  # PSM id doesn't matter for this test

    # Fake precursor info returns IDs that do not match the regex below
    def fake_get_precursor_info(path):
        class FakePrecursor:
            def __init__(self, mz, rt, im):
                self.mz = mz
                self.rt = rt
                self.im = im

        return {"specA": FakePrecursor(100.0, 1.0, 0.1), "specB": FakePrecursor(200.0, 2.0, 0.2)}

    monkeypatch.setattr("ms2rescore.parse_spectra.get_precursor_info", fake_get_precursor_info)
    # Avoid filesystem validation in infer_spectrum_path by returning a dummy path
    monkeypatch.setattr("ms2rescore.parse_spectra.infer_spectrum_path", lambda cfg, rn: "dummy")

    # Request a data type that is missing (retention_time) so the function will
    # parse spectrum files and therefore apply the spectrum_id_pattern.
    with pytest.raises(MS2RescoreConfigurationError):
        add_precursor_values(
            psm_list,
            {MSDataType.retention_time},
            spectrum_path="/not/used",
            spectrum_id_pattern=r"scan:(\d+):.*",
        )


@patch("ms2rescore.parse_spectra.get_precursor_info")
@patch("ms2rescore.parse_spectra.infer_spectrum_path")
def test_add_precursor_values(
    mock_infer_spectrum_path, mock_get_precursor_info, mock_psm_list, mock_precursor_info
):
    mock_infer_spectrum_path.return_value = "test_data/test_spectrum_file.mgf"
    mock_get_precursor_info.return_value = mock_precursor_info

    required_data_types = {
        MSDataType.retention_time,
        MSDataType.ion_mobility,
        MSDataType.precursor_mz,
    }
    available_ms_data = add_precursor_values(
        mock_psm_list, required_data_types, spectrum_path="test_data"
    )

    assert MSDataType.retention_time in available_ms_data
    assert MSDataType.ion_mobility in available_ms_data
    assert MSDataType.precursor_mz in available_ms_data

    for psm in mock_psm_list:
        assert psm.retention_time is not None
        assert psm.ion_mobility is not None
        assert psm.precursor_mz is not None


@patch("ms2rescore.parse_spectra.get_precursor_info")
@patch("ms2rescore.parse_spectra.infer_spectrum_path")
def test_add_precursor_values_missing_im(
    mock_infer_spectrum_path,
    mock_get_precursor_info,
    mock_psm_list,
    mock_precursor_info_missing_im,
):
    mock_infer_spectrum_path.return_value = "test_data/test_spectrum_file.mgf"
    mock_get_precursor_info.return_value = mock_precursor_info_missing_im

    required_data_types = {
        MSDataType.retention_time,
        MSDataType.ion_mobility,
        MSDataType.precursor_mz,
    }
    available_ms_data = add_precursor_values(
        mock_psm_list, required_data_types, spectrum_path="test_data"
    )

    assert MSDataType.retention_time in available_ms_data
    assert MSDataType.ion_mobility not in available_ms_data
    assert MSDataType.precursor_mz in available_ms_data

    for psm in mock_psm_list:
        assert psm.retention_time is not None
        assert psm.ion_mobility is None
        assert psm.precursor_mz is not None


@patch("ms2rescore.parse_spectra.get_precursor_info")
@patch("ms2rescore.parse_spectra.infer_spectrum_path")
def test_get_precursor_values(
    mock_infer_spectrum_path, mock_get_precursor_info, mock_psm_list, mock_precursor_info
):
    mock_infer_spectrum_path.return_value = "test_data/test_spectrum_file.mgf"
    mock_get_precursor_info.return_value = mock_precursor_info

    mzs, rts, ims = _get_precursor_values(mock_psm_list, "test_data", None)

    expected_mzs = np.array([529.7935187324, 651.83])
    expected_rts = np.array([10.5, 12.3])
    expected_ims = np.array([1.0, 1.2])

    np.testing.assert_array_equal(mzs, expected_mzs)
    np.testing.assert_array_equal(rts, expected_rts)
    np.testing.assert_array_equal(ims, expected_ims)


@patch("ms2rescore.parse_spectra.get_precursor_info")
@patch("ms2rescore.parse_spectra.infer_spectrum_path")
def test_get_precursor_values_missing_spectrum_id(
    mock_infer_spectrum_path,
    mock_get_precursor_info,
    mock_psm_list,
    mock_precursor_info_incomplete,
):
    mock_infer_spectrum_path.return_value = "test_data/test_spectrum_file.mgf"
    mock_get_precursor_info.return_value = mock_precursor_info_incomplete

    with pytest.raises(MS2RescoreConfigurationError):
        _get_precursor_values(mock_psm_list, "test_data", None)


def test_add_precursor_values_no_missing_data():
    """Test early return when all required data is already available."""
    psm_list = PSMList(
        psm_list=[
            PSM(
                peptidoform="PEPTIDE/2",
                run="run1",
                spectrum_id="spectrum1",
                retention_time=10.5,
                ion_mobility=1.0,
                precursor_mz=529.79,
            ),
        ]
    )
    required_data_types = {MSDataType.retention_time, MSDataType.ion_mobility}

    available_ms_data = add_precursor_values(psm_list, required_data_types)

    assert MSDataType.retention_time in available_ms_data
    assert MSDataType.ion_mobility in available_ms_data
    assert MSDataType.precursor_mz in available_ms_data  # Available but not required


def test_add_precursor_values_no_spectrum_path_error():
    """Test that error is raised when spectrum path is needed but not provided."""
    psm_list = PSMList(
        psm_list=[
            PSM(
                peptidoform="PEPTIDE/2",
                run="run1",
                spectrum_id="spectrum1",
                retention_time=None,
                ion_mobility=None,
                precursor_mz=None,
            ),
        ]
    )
    required_data_types = {MSDataType.retention_time}

    with pytest.raises(SpectrumParsingError, match="Spectrum path must be provided"):
        add_precursor_values(psm_list, required_data_types)


def test_add_precursor_values_ms2_spectra_availability():
    """Test that MS2 spectra availability depends on spectrum_path."""
    psm_list = PSMList(
        psm_list=[
            PSM(
                peptidoform="PEPTIDE/2",
                run="run1",
                spectrum_id="spectrum1",
                retention_time=10.5,
                ion_mobility=1.0,
                precursor_mz=529.79,
            ),
        ]
    )

    # Without spectrum path - MS2 spectra not available
    required_data_types = {MSDataType.retention_time}
    available_ms_data = add_precursor_values(psm_list, required_data_types)
    assert MSDataType.ms2_spectra not in available_ms_data


def test_spectrum_parsing_error():
    with pytest.raises(SpectrumParsingError):
        raise SpectrumParsingError("Test error message")
