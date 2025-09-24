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


def add_precursor_values(
    psm_list: PSMList,
    spectrum_path: str,
    spectrum_id_pattern: Optional[str] = None,
    required_ms_data: Optional[Set[MSDataType]] = None,
) -> Set[MSDataType]:
    """
    Add precursor m/z, retention time, and ion mobility values to a PSM list.

    Parameters
    ----------
    psm_list
        PSM list to add precursor values to.
    spectrum_path
        Path to the spectrum files.
    spectrum_id_pattern
        Regular expression pattern to extract spectrum IDs from file names. If provided, the
        pattern must contain a single capturing group that matches the spectrum ID. Default is
        None.
    required_ms_data
        Set of MS data types required for feature generation. If provided, only the missing
        precursor values will be added to the PSM list. Default is None.

    Returns
    -------
    available_ms_data
        Set of available MS data types in the PSM list.

    """
    # Check if precursor values are missing in PSM list
    missing_ms_data = {}

    rt_missing = any(v is None or v == 0 or np.isnan(v) for v in psm_list["retention_time"])
    im_missing = any(v is None or v == 0 or np.isnan(v) for v in psm_list["ion_mobility"])
    mz_missing = any(v is None or v == 0 or np.isnan(v) for v in psm_list["precursor_mz"])

    if required_ms_data:
        if MSDataType.retention_time in required_ms_data:
            missing_ms_data[MSDataType.retention_time] = rt_missing
        if MSDataType.ion_mobility in required_ms_data:
            missing_ms_data[MSDataType.ion_mobility] = im_missing
        if MSDataType.precursor_mz in required_ms_data:
            missing_ms_data[MSDataType.precursor_mz] = mz_missing
    else:
        missing_ms_data = {
            MSDataType.retention_time: rt_missing,
            MSDataType.ion_mobility: im_missing,
            MSDataType.precursor_mz: mz_missing,
        }
    if not any(missing_ms_data.values()):
        # No missing precursor values
        LOGGER.debug("No missing precursor values in PSM list.")
        return {
            MSDataType.ms2_spectra,  # Assume MS2 spectra are always present
            MSDataType.retention_time,
            MSDataType.ion_mobility,
            MSDataType.precursor_mz,
        }
    else:
        # Get precursor values from spectrum files
        LOGGER.info("Parsing precursor info from spectrum files...")
        mz, rt, im = _get_precursor_values(psm_list, spectrum_path, spectrum_id_pattern)
        mz_found, rt_found, im_found = np.all(mz != 0.0), np.all(rt != 0.0), np.all(im != 0.0)
        # ms2rescore_rs always returns 0.0 for missing values

        # Update PSM list with missing precursor values
        if rt_missing and rt_found:
            LOGGER.debug(
                "Missing retention time values in PSM list. Updating from spectrum files."
            )
            psm_list["retention_time"] = rt
        if im_missing and im_found:
            LOGGER.debug("Missing ion mobility values in PSM list. Updating from spectrum files.")
            psm_list["ion_mobility"] = im
        if mz_missing and mz_found:
            LOGGER.debug("Missing precursor m/z values in PSM list. Updating from spectrum files.")
            psm_list["precursor_mz"] = mz
        else:
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

    # Return available MS data types
    available_ms_data = {
        MSDataType.ms2_spectra,  # Assume MS2 spectra are always present
        MSDataType.retention_time if not rt_missing or rt_found else None,
        MSDataType.ion_mobility if not im_missing or im_found else None,
        MSDataType.precursor_mz if not mz_missing or mz_found else None,
    }
    available_ms_data.discard(None)

    return available_ms_data


def _get_precursor_values(
    psm_list: PSMList, spectrum_path: str, spectrum_id_pattern: str
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
                    compiled_pattern.search(spectrum_id).group(1): precursor
                    for spectrum_id, precursor in precursors.items()
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
