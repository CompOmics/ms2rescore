"""Parse MGF files."""

import logging
import re
from enum import Enum
from typing import Optional, Set

import numpy as np
from ms2pip._spectrum_processing import proforma_to_mass_shift
from ms2rescore_rs import MS2Spectrum, annotate_ms2_spectra, get_ms2_spectra
from psm_utils import PSMList
from rich.progress import track

from ms2rescore.exceptions import MS2RescoreConfigurationError, MS2RescoreError
from ms2rescore.utils import infer_spectrum_path

LOGGER = logging.getLogger(__name__)


class MSDataType(str, Enum):
    """Enum for MS data types required for feature generation."""

    retention_time = "retention time"
    ion_mobility = "ion mobility"
    precursor_mz = "precursor m/z"
    ms2_spectra = "MS2 spectra"

    def __str__(self):
        return self.value


ALL_MS_DATA_TYPES: Set[MSDataType] = {
    MSDataType.retention_time,
    MSDataType.ion_mobility,
    MSDataType.precursor_mz,
    MSDataType.ms2_spectra,
}


def add_precursor_values(
    psm_list: PSMList,
    required_data_types: Set[MSDataType],
    spectrum_path: Optional[str] = None,
    spectrum_id_pattern: Optional[str] = None,
) -> Set[MSDataType]:
    """
    Add precursor m/z, retention time, and ion mobility values to a PSM list.

    Parameters
    ----------
    psm_list
        PSM list to add precursor values to.
    required_data_types
        Set of MS data types required for feature generation. Only the missing precursor values
        will be added to the PSM list.
    spectrum_path
        Path to the spectrum files. Default is None.
    spectrum_id_pattern
        Regular expression pattern to extract spectrum IDs from file names. If provided, the
        pattern must contain a single capturing group that matches the spectrum ID. Default is
        None.

    Returns
    -------
    available_data_types
        Set of available MS data types in the PSM list.

    """
    # Check which data types are missing
    # Missing if: all values are 0, OR any values are None/NaN
    missing_data_types = {MSDataType.ms2_spectra}  # Always missing until spectrum files are parsed

    rt_values = np.asarray(psm_list["retention_time"])
    if np.any(np.isnan(rt_values)) or np.all(rt_values == 0):
        missing_data_types.add(MSDataType.retention_time)

    im_values = np.asarray(psm_list["ion_mobility"])
    if np.any(np.isnan(im_values)) or np.all(im_values == 0):
        missing_data_types.add(MSDataType.ion_mobility)

    mz_values = np.asarray(psm_list["precursor_mz"])
    if np.any(np.isnan(mz_values)) or np.all(mz_values == 0):
        missing_data_types.add(MSDataType.precursor_mz)

    # Find data types that are both missing and required
    data_types_to_parse = missing_data_types & required_data_types

    # If no data types need to be parsed, return available data types
    if not data_types_to_parse:
        LOGGER.debug("All required data types are already available.")
        # Use same logic as final return: available = all - missing + found (found is empty here)
        found_data_types: set[MSDataType] = set()  # No spectrum file processing done
        available_data_types = ALL_MS_DATA_TYPES - missing_data_types | found_data_types
        return available_data_types

    # If no spectrum path is provided, cannot parse missing precursor values
    elif spectrum_path is None:
        raise SpectrumParsingError(
            "Spectrum path must be provided to parse precursor values that are not present in the"
            " PSM list."
        )
    else:
        LOGGER.debug(
            "Missing required data types: %s. Parsing from spectrum files.",
            ", ".join(str(dt) for dt in data_types_to_parse),
        )

    # Get precursor values from spectrum files
    LOGGER.info("Parsing precursor info from spectrum files...")
    _add_precursor_values(psm_list, spectrum_path, spectrum_id_pattern)

    # Extract precursor values from MS2 spectrum objects in a single pass
    precursor_data = [
        (ms2.precursor.rt, ms2.precursor.im, ms2.precursor.mz) for ms2 in psm_list["spectrum"]
    ]
    rts, ims, mzs = map(np.array, zip(*precursor_data))

    # Determine which data types were successfully found in spectrum files
    found_data_types = {MSDataType.ms2_spectra}

    # Add found data types: if missing and all zeros, raise error
    if not np.all(rts == 0.0):
        found_data_types.add(MSDataType.retention_time)
        if MSDataType.retention_time in data_types_to_parse:
            LOGGER.debug(
                "Missing retention time values in PSM list. Updating from spectrum files."
            )
            psm_list["retention_time"] = rts
    elif MSDataType.retention_time in data_types_to_parse:
        raise SpectrumParsingError(
            "Retention time values are required but not available in spectrum files "
            "(all values are zero)."
        )

    if not np.all(ims == 0.0):
        found_data_types.add(MSDataType.ion_mobility)
        if MSDataType.ion_mobility in data_types_to_parse:
            LOGGER.debug("Missing ion mobility values in PSM list. Updating from spectrum files.")
            psm_list["ion_mobility"] = ims
    elif MSDataType.ion_mobility in data_types_to_parse:
        raise SpectrumParsingError(
            "Ion mobility values are required but not available in spectrum files "
            "(all values are zero)."
        )

    if not np.all(mzs == 0.0):
        found_data_types.add(MSDataType.precursor_mz)
        if MSDataType.precursor_mz in data_types_to_parse:
            LOGGER.debug("Missing precursor m/z values in PSM list. Updating from spectrum files.")
            psm_list["precursor_mz"] = mzs
    elif MSDataType.precursor_mz in data_types_to_parse:
        raise SpectrumParsingError(
            "Precursor m/z values are required but not available in spectrum files "
            "(all values are zero)."
        )

    # Return available data types: (all types - missing types) + found types
    available_data_types = ALL_MS_DATA_TYPES - missing_data_types | found_data_types
    return available_data_types


def _acquire_observed_spectra_dict(
    ms2: list[MS2Spectrum], pattern: str, spectrum_ids: list[str]
) -> dict[str, MS2Spectrum]:
    """Apply spectrum ID pattern to precursor IDs."""
    # Map precursor IDs using regex pattern
    compiled_pattern = re.compile(pattern)
    spectrum_ids_set = set(spectrum_ids)  # For faster lookup

    ms2_observed_spectra_mapping = {
        match.group(1): ms2_spectrum
        for ms2_spectrum in ms2
        if (match := compiled_pattern.search(str(ms2_spectrum.identifier))) is not None
        and match.group(1) in spectrum_ids_set
    }

    # Validate that any IDs were matched
    if not ms2_observed_spectra_mapping:
        raise MS2RescoreConfigurationError(
            "'spectrum_id_pattern' did not match any spectrum-file IDs. Please check and try "
            "again. See "
            "https://ms2rescore.readthedocs.io/en/stable/userguide/configuration/#mapping-psms-to-spectra "
            "for more information."
        )

    return ms2_observed_spectra_mapping


def _add_precursor_values(
    psm_list: PSMList, spectrum_path: str, spectrum_id_pattern: Optional[str] = None
) -> None:
    """Get precursor m/z, RT, and IM from spectrum files."""
    # Iterate over different runs in PSM list
    if spectrum_id_pattern is None:
        spectrum_id_pattern = r"^(.*)$"  # Match entire identifier if no pattern provided

    for run_name in track(set(psm_list["run"])):
        run_mask = psm_list["run"] == run_name
        psm_list_run = psm_list[run_mask]
        spectrum_file = infer_spectrum_path(spectrum_path, run_name)

        LOGGER.debug("Reading spectrum file: '%s'", spectrum_file)
        ms2_spectra: list[MS2Spectrum] = get_ms2_spectra(str(spectrum_file))

        # Parse spectrum IDs with regex pattern if provided
        ms2_spectra_dict = _acquire_observed_spectra_dict(
            ms2_spectra, spectrum_id_pattern, psm_list_run["spectrum_id"]
        )

        try:
            psm_list_run["spectrum"] = [
                ms2_spectra_dict[spec_id] for spec_id in psm_list_run["spectrum_id"]
            ]
        except KeyError as e:
            raise SpectrumParsingError(
                f"Could not find spectrum {e} in spectrum file '{spectrum_file}'. Please "
                "check the 'spectrum_id_pattern' and 'psm_id_pattern' configuration options. See "
                "https://ms2rescore.readthedocs.io/en/stable/userguide/configuration/#mapping-psms-to-spectra "
                "for more information."
            ) from e


def annotate_spectra(
    psm_list: PSMList,
    fragmentation_model: str = "cidhcd",
    ms2_tolerance: float = 0.02,
    ms2_tolerance_mode: str = "Da",
) -> None:
    """
    Annotate MS2 spectra with fragment ion matches, in place.

    Replaces raw ``MS2Spectrum`` objects in ``psm_list["spectrum"]`` with
    ``AnnotatedMS2Spectrum`` objects that include peak annotations. These annotated
    spectra can then be consumed by any feature generator that needs fragment ion
    information (e.g., MS2, MS2PIP). In particular, the unified
    ``ms2pip.correlate()`` API can reuse these annotations directly when they are
    already attached to ``psm.spectrum``.

    Parameters
    ----------
    psm_list
        PSM list with loaded MS2 spectra in the ``spectrum`` field.
    fragmentation_model
        Fragmentation model: ``cidhcd``, ``etd``, ``ethcd``, or ``all``.
    ms2_tolerance
        Fragment mass tolerance value.
    ms2_tolerance_mode
        Fragment mass tolerance mode: ``ppm`` or ``Da``.

    """
    LOGGER.info("Annotating MS2 spectra based on search engine identifications...")

    spectra = list(psm_list["spectrum"])
    # Convert modification labels to numeric mass shifts so rustyms can parse them
    # regardless of modification name convention (Unimod name, accession, formula, etc.).
    proformas = [proforma_to_mass_shift(psm.peptidoform) for psm in psm_list]

    annotated = annotate_ms2_spectra(
        spectra=spectra,
        proformas=proformas,
        fragmentation_model=fragmentation_model,
        mass_mode="monoisotopic",
        tolerance_value=ms2_tolerance,
        tolerance_mode=ms2_tolerance_mode,
    )

    psm_list["spectrum"] = annotated


class SpectrumParsingError(MS2RescoreError):
    """Error while parsing spectrum file."""

    pass
