from unittest.mock import MagicMock, patch

import pytest
from psm_utils import PSM, PSMList

from ms2rescore.exceptions import MS2RescoreConfigurationError
from ms2rescore.parse_spectra import (
    MSDataType,
    SpectrumParsingError,
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


def _make_mock_ms2_spectrum(identifier, mz, rt, im):
    """Create a mock MS2Spectrum with the given precursor values."""
    precursor = MagicMock(mz=mz, rt=rt, im=im)
    spectrum = MagicMock(identifier=identifier, precursor=precursor)
    return spectrum


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
def mock_ms2_spectra():
    return [
        _make_mock_ms2_spectrum("spectrum1", mz=529.7935187324, rt=10.5, im=1.0),
        _make_mock_ms2_spectrum("spectrum2", mz=651.83, rt=12.3, im=1.2),
    ]


@pytest.fixture
def mock_ms2_spectra_missing_im():
    return [
        _make_mock_ms2_spectrum("spectrum1", mz=529.7935187324, rt=10.5, im=0.0),
        _make_mock_ms2_spectrum("spectrum2", mz=651.83, rt=12.3, im=0.0),
    ]


def test_spectrum_id_pattern_nonmatching(monkeypatch):
    """If the provided spectrum_id_pattern doesn't match any spectrum-file IDs, raise."""
    psm_list = psm_list_factory(["scan:1:fileA"])

    # Fake spectra with IDs that do not match the regex below
    fake_spectra = [
        _make_mock_ms2_spectrum("specA", mz=100.0, rt=1.0, im=0.1),
        _make_mock_ms2_spectrum("specB", mz=200.0, rt=2.0, im=0.2),
    ]

    monkeypatch.setattr(
        "ms2rescore.parse_spectra.get_ms2_spectra", lambda path: fake_spectra
    )
    monkeypatch.setattr(
        "ms2rescore.parse_spectra.infer_spectrum_path", lambda cfg, rn: "dummy"
    )

    with pytest.raises(MS2RescoreConfigurationError):
        add_precursor_values(
            psm_list,
            {MSDataType.retention_time},
            spectrum_path="/not/used",
            spectrum_id_pattern=r"scan:(\d+):.*",
        )


@patch("ms2rescore.parse_spectra.get_ms2_spectra")
@patch("ms2rescore.parse_spectra.infer_spectrum_path")
def test_add_precursor_values(
    mock_infer_spectrum_path, mock_get_ms2_spectra, mock_psm_list, mock_ms2_spectra
):
    mock_infer_spectrum_path.return_value = "test_data/test_spectrum_file.mgf"
    mock_get_ms2_spectra.return_value = mock_ms2_spectra

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


@patch("ms2rescore.parse_spectra.get_ms2_spectra")
@patch("ms2rescore.parse_spectra.infer_spectrum_path")
def test_add_precursor_values_missing_im(
    mock_infer_spectrum_path,
    mock_get_ms2_spectra,
    mock_psm_list,
    mock_ms2_spectra_missing_im,
):
    """When IM is all zeros in spectrum files but not required, it should not be available."""
    mock_infer_spectrum_path.return_value = "test_data/test_spectrum_file.mgf"
    mock_get_ms2_spectra.return_value = mock_ms2_spectra_missing_im

    # Only require RT and precursor_mz, not ion_mobility
    required_data_types = {
        MSDataType.retention_time,
        MSDataType.precursor_mz,
    }
    available_ms_data = add_precursor_values(
        mock_psm_list, required_data_types, spectrum_path="test_data"
    )

    assert MSDataType.retention_time in available_ms_data
    assert MSDataType.ion_mobility not in available_ms_data
    assert MSDataType.precursor_mz in available_ms_data


@patch("ms2rescore.parse_spectra.get_ms2_spectra")
@patch("ms2rescore.parse_spectra.infer_spectrum_path")
def test_add_precursor_values_missing_im_required_raises(
    mock_infer_spectrum_path,
    mock_get_ms2_spectra,
    mock_psm_list,
    mock_ms2_spectra_missing_im,
):
    """When IM is required but all zeros in spectrum files, raise SpectrumParsingError."""
    mock_infer_spectrum_path.return_value = "test_data/test_spectrum_file.mgf"
    mock_get_ms2_spectra.return_value = mock_ms2_spectra_missing_im

    required_data_types = {
        MSDataType.retention_time,
        MSDataType.ion_mobility,
        MSDataType.precursor_mz,
    }
    with pytest.raises(SpectrumParsingError, match="Ion mobility values are required"):
        add_precursor_values(
            mock_psm_list, required_data_types, spectrum_path="test_data"
        )


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


def test_add_precursor_values_ms2_spectra_not_available_without_spectrum_path():
    """Test that MS2 spectra are not available when no spectrum files are parsed."""
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

    # All required data present, so no spectrum files are parsed
    required_data_types = {MSDataType.retention_time}
    available_ms_data = add_precursor_values(psm_list, required_data_types)
    # ms2_spectra is only available when spectrum files were actually parsed
    assert MSDataType.ms2_spectra not in available_ms_data


def test_spectrum_parsing_error():
    with pytest.raises(SpectrumParsingError):
        raise SpectrumParsingError("Test error message")
