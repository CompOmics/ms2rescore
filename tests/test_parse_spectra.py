from unittest.mock import MagicMock, patch

import pytest
from ms2pip._spectrum_processing import _to_observed_spectrum
from ms2rescore_rs import MS2Spectrum, Precursor
from psm_utils import PSM, PSMList

from ms2rescore.exceptions import MS2RescoreConfigurationError
from ms2rescore.parse_spectra import (
    MSDataType,
    SpectrumParsingError,
    add_precursor_values,
    annotate_spectra,
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


@patch("ms2rescore.parse_spectra.get_ms2_spectra")
@patch("ms2rescore.parse_spectra.infer_spectrum_path")
def test_add_precursor_values_missing_precursor(
    mock_infer_spectrum_path, mock_get_ms2_spectra, mock_psm_list
):
    """MS2Spectrum objects with precursor=None should be treated like empty Precursor()."""
    mock_infer_spectrum_path.return_value = "test_data/test_spectrum_file.mgf"
    mock_get_ms2_spectra.return_value = [
        MS2Spectrum(identifier="spectrum1", mz=[100.0, 200.0], intensity=[10.0, 20.0]),
        MS2Spectrum(identifier="spectrum2", mz=[150.0, 250.0], intensity=[15.0, 25.0]),
    ]

    available_ms_data = add_precursor_values(
        mock_psm_list, {MSDataType.ms2_spectra}, spectrum_path="test_data"
    )

    assert available_ms_data == {MSDataType.ms2_spectra}


@patch("ms2rescore.parse_spectra.get_ms2_spectra")
@patch("ms2rescore.parse_spectra.infer_spectrum_path")
def test_add_precursor_values_missing_precursor_required_raises(
    mock_infer_spectrum_path, mock_get_ms2_spectra, mock_psm_list
):
    """When precursor=None and a precursor value is required, raise SpectrumParsingError."""
    mock_infer_spectrum_path.return_value = "test_data/test_spectrum_file.mgf"
    mock_get_ms2_spectra.return_value = [
        MS2Spectrum(identifier="spectrum1", mz=[100.0, 200.0], intensity=[10.0, 20.0]),
        MS2Spectrum(identifier="spectrum2", mz=[150.0, 250.0], intensity=[15.0, 25.0]),
    ]

    with pytest.raises(SpectrumParsingError, match="Retention time values are required"):
        add_precursor_values(mock_psm_list, {MSDataType.retention_time}, spectrum_path="test_data")


@pytest.mark.parametrize("extended", [False, True])
def test_annotate_spectra_missing_precursor_matches_empty_and_explicit(extended):
    """annotate_spectra handles precursor=None identically to Precursor() and supports ms2pip."""
    mz = [324.1554, 393.1769, 435.1980, 644.2175, 999.0]
    intensity = [10.0, 20.0, 30.0, 40.0, 100.0]

    raw_none = MS2Spectrum(identifier="scan=1", mz=mz, intensity=intensity)
    raw_empty = MS2Spectrum(identifier="scan=1", mz=mz, intensity=intensity, precursor=Precursor())
    explicit_precursor = Precursor(mz=484.19, charge=2, rt=51.2, im=0.95)
    raw_explicit = MS2Spectrum(
        identifier="scan=1", mz=mz, intensity=intensity, precursor=explicit_precursor
    )

    psm_none = PSMList(
        psm_list=[
            PSM(
                peptidoform="PEPS[Phospho]TIDE/2",
                spectrum_id="scan=1",
                run="run1",
                spectrum=raw_none,
            )
        ]
    )
    psm_empty = PSMList(
        psm_list=[
            PSM(
                peptidoform="PEPS[Phospho]TIDE/2",
                spectrum_id="scan=1",
                run="run1",
                spectrum=raw_empty,
            )
        ]
    )
    psm_explicit = PSMList(
        psm_list=[
            PSM(
                peptidoform="PEPS[Phospho]TIDE/2",
                spectrum_id="scan=1",
                run="run1",
                spectrum=raw_explicit,
            )
        ]
    )

    annotate_spectra(psm_none, "cidhcd", 20.0, "ppm", extended=extended)
    annotate_spectra(psm_empty, "cidhcd", 20.0, "ppm", extended=extended)
    annotate_spectra(psm_explicit, "cidhcd", 20.0, "ppm", extended=extended)

    # Caller's original MS2Spectrum object is not mutated in place
    assert raw_none.precursor is None

    ann_none = psm_none[0].spectrum
    ann_empty = psm_empty[0].spectrum
    ann_explicit = psm_explicit[0].spectrum

    # Annotated spectrum carries a default Precursor() compatible with ms2pip._to_observed_spectrum
    assert ann_none.precursor is not None
    assert (
        ann_none.precursor.mz,
        ann_none.precursor.charge,
        ann_none.precursor.rt,
        ann_none.precursor.im,
        ann_none.precursor.intensity,
    ) == (0.0, 0, 0.0, 0.0, 0.0)

    obs_none = _to_observed_spectrum(ann_none)
    obs_empty = _to_observed_spectrum(ann_empty)
    assert obs_none.precursor_mz == obs_empty.precursor_mz == 0.0
    assert obs_none.precursor_charge == obs_empty.precursor_charge == 0
    assert obs_none.retention_time == obs_empty.retention_time == 0.0

    # Peak annotations match bit-for-bit across precursor=None, Precursor(), and explicit Precursor
    assert len(ann_none.backbone) > 0
    assert repr(ann_none.backbone) == repr(ann_empty.backbone) == repr(ann_explicit.backbone)
    assert repr(ann_none.extended) == repr(ann_empty.extended) == repr(ann_explicit.extended)

    # Explicit precursor metadata on the control spectrum is preserved unchanged
    assert ann_explicit.precursor.mz == pytest.approx(484.19)
    assert ann_explicit.precursor.charge == 2
    assert ann_explicit.precursor.rt == pytest.approx(51.2)
    assert ann_explicit.precursor.im == pytest.approx(0.95)


def test_annotate_spectra_missing_precursor_and_missing_charge():
    """When both ProForma charge and precursor are absent, no charge is fabricated."""
    mz = [227.1026, 324.1554, 425.2031]
    intensity = [200.0, 300.0, 400.0]

    psm_none = PSMList(
        psm_list=[
            PSM(
                peptidoform="PEPTIDE",
                spectrum_id="scan=1",
                run="run1",
                spectrum=MS2Spectrum(identifier="scan=1", mz=mz, intensity=intensity),
            )
        ]
    )
    psm_empty = PSMList(
        psm_list=[
            PSM(
                peptidoform="PEPTIDE",
                spectrum_id="scan=1",
                run="run1",
                spectrum=MS2Spectrum(
                    identifier="scan=1", mz=mz, intensity=intensity, precursor=Precursor()
                ),
            )
        ]
    )
    psm_fallback = PSMList(
        psm_list=[
            PSM(
                peptidoform="PEPTIDE",
                spectrum_id="scan=1",
                run="run1",
                spectrum=MS2Spectrum(
                    identifier="scan=1",
                    mz=mz,
                    intensity=intensity,
                    precursor=Precursor(mz=400.19, charge=2),
                ),
            )
        ]
    )

    annotate_spectra(psm_none)
    annotate_spectra(psm_empty)
    annotate_spectra(psm_fallback)

    assert psm_none[0].spectrum.backbone == []
    assert psm_empty[0].spectrum.backbone == []
    assert len(psm_fallback[0].spectrum.backbone) == 3
