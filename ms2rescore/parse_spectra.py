"""Parse MGF files."""

import logging
import re
from enum import Enum
from itertools import chain
from typing import Optional, Set, Tuple

import numpy as np
from ms2rescore_rs import get_precursor_info
from psm_utils import PSMList

from ms2rescore.exceptions import MS2RescoreError
from ms2rescore.utils import infer_spectrum_path

LOGGER = logging.getLogger(__name__)


class MSDataType(str, Enum):
    """Enum for MS data types required for feature generation."""

    retention_time = "retention time"
    ion_mobility = "ion mobility"
    precursor_mz = "precursor m/z"
    ms2_spectra = "MS2 spectra"

    # Mimic behavior of StrEnum (Python >=3.11)
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
    spectrum_path: str | None = None,
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
    missing_data_types = set()
    if spectrum_path is None:
        missing_data_types.add(MSDataType.ms2_spectra)

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

    # Get precursor values from spectrum files
    LOGGER.info("Parsing precursor info from spectrum files...")
    mz, rt, im = _get_precursor_values(psm_list, spectrum_path, spectrum_id_pattern)

    # Determine which data types were successfully found in spectrum files
    # ms2rescore_rs always returns 0.0 for missing values
    found_data_types = {MSDataType.ms2_spectra}  # MS2 spectra available when processing files
    if np.all(rt != 0.0):
        found_data_types.add(MSDataType.retention_time)
    if np.all(im != 0.0):
        found_data_types.add(MSDataType.ion_mobility)
    if np.all(mz != 0.0):
        found_data_types.add(MSDataType.precursor_mz)

    # Update PSM list with missing precursor values that were found
    update_types = data_types_to_parse & found_data_types

    if MSDataType.retention_time in update_types:
        LOGGER.debug("Missing retention time values in PSM list. Updating from spectrum files.")
        psm_list["retention_time"] = rt
    if MSDataType.ion_mobility in update_types:
        LOGGER.debug("Missing ion mobility values in PSM list. Updating from spectrum files.")
        psm_list["ion_mobility"] = im
    if MSDataType.precursor_mz in update_types:
        LOGGER.debug("Missing precursor m/z values in PSM list. Updating from spectrum files.")
        psm_list["precursor_mz"] = mz
    elif (
        MSDataType.precursor_mz not in missing_data_types
        and MSDataType.precursor_mz in found_data_types
    ):
        # Check if precursor m/z values are consistent between PSMs and spectrum files
        mz_diff = np.abs(psm_list["precursor_mz"] - mz)
        if np.mean(mz_diff) > 1e-2:
            LOGGER.warning(
                "Mismatch between precursor m/z values in PSM list and spectrum files (mean "
                "difference exceeds 0.01 Da). Please ensure that the correct spectrum files are "
                "provided and that the `spectrum_id_pattern` and `psm_id_pattern` options are "
                "configured correctly. See "
                "https://ms2rescore.readthedocs.io/en/stable/userguide/configuration/#mapping-psms-to-spectra "
                "for more information."
            )

    # Return available data types: (all types - missing types) + found types
    available_data_types = ALL_MS_DATA_TYPES - missing_data_types | found_data_types
    return available_data_types


def _get_precursor_values(
    psm_list: PSMList, spectrum_path: str, spectrum_id_pattern: str | None = None
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Get precursor m/z, RT, and IM from spectrum files."""
    # Iterate over different runs in PSM list
    precursor_dict = dict()
    psm_dict = psm_list.get_psm_dict()
    for runs in psm_dict.values():
        for run_name, psms in runs.items():
            psm_list_run = PSMList(psm_list=list(chain.from_iterable(psms.values())))
            spectrum_file = infer_spectrum_path(spectrum_path, run_name)

            LOGGER.debug("Reading spectrum file: '%s'", spectrum_file)
            precursors = get_precursor_info(str(spectrum_file))

            # Parse spectrum IDs with regex pattern if provided
            if spectrum_id_pattern:
                compiled_pattern = re.compile(spectrum_id_pattern)
                precursors = {
                    match.group(1): precursor
                    for spectrum_id, precursor in precursors.items()
                    if (match := compiled_pattern.search(spectrum_id)) is not None
                }

            # Ensure all PSMs have a precursor values
            for psm in psm_list_run:
                if psm.spectrum_id not in precursors:
                    raise SpectrumParsingError(
                        "Mismatch between PSM and spectrum file IDs. Could find precursor values "
                        f"for PSM with ID {psm.spectrum_id} in run {run_name}.\n"
                        "Please check that the `spectrum_id_pattern` and `psm_id_pattern` options "
                        "are configured correctly. See "
                        "https://ms2rescore.readthedocs.io/en/stable/userguide/configuration/#mapping-psms-to-spectra"
                        " for more information.\n"
                        f"Example ID from PSM file: {psm.spectrum_id}\n"
                        f"Example ID from spectrum file: {list(precursors.keys())[0]}"
                    )

            # Store precursor values in dictionary
            precursor_dict[run_name] = precursors

    # Reshape precursor values into arrays matching PSM list
    mzs = np.fromiter((precursor_dict[psm.run][psm.spectrum_id].mz for psm in psm_list), float)
    rts = np.fromiter((precursor_dict[psm.run][psm.spectrum_id].rt for psm in psm_list), float)
    ims = np.fromiter((precursor_dict[psm.run][psm.spectrum_id].im for psm in psm_list), float)

    return mzs, rts, ims


class SpectrumParsingError(MS2RescoreError):
    """Error parsing retention time from spectrum file."""

    pass
