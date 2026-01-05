"""
IM2Deep ion mobility-based feature generator.

IM2Deep is a fully modification-aware peptide ion mobility predictor. It uses a deep convolutional
neural network to predict retention times based on the atomic composition of the (modified) amino
acid residues in the peptide. See
`github.com/compomics/IM2Deep <https://github.com/compomics/IM2Deep>`_ for more information.

"""

import logging
import os
from inspect import getfullargspec
from typing import List

import numpy as np
import pandas as pd
from im2deep.utils import im2ccs
from im2deep.im2deep import predict_ccs, REFERENCE_DATASET_PATH
from im2deep.calibrate import calculate_ccs_shift
from psm_utils import PSMList, Peptidoform

from ms2rescore.feature_generators.base import FeatureGeneratorBase
from ms2rescore.parse_spectra import MSDataType

os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"
logger = logging.getLogger(__name__)


class IM2DeepFeatureGenerator(FeatureGeneratorBase):
    """IM2Deep collision cross section feature generator."""

    required_ms_data = {MSDataType.ion_mobility}

    def __init__(
        self,
        *args,
        processes: int = 1,
        **kwargs,
    ):
        """
        Initialize the IM2DeepFeatureGenerator.

        Parameters
        ----------
        processes : int, optional
            Number of parallel processes to use for IM2Deep predictions. Default is 1.
        **kwargs : dict, optional
            Additional keyword arguments to `im2deep.predict_ccs`.

        """
        super().__init__(*args, **kwargs)

        self._verbose = logger.getEffectiveLevel() <= logging.DEBUG

        # Remove any kwargs that are not IM2Deep arguments
        self.im2deep_kwargs = kwargs or {}
        self.im2deep_kwargs = {
            k: v for k, v in self.im2deep_kwargs.items() if k in getfullargspec(predict_ccs).args
        }
        self.im2deep_kwargs["n_jobs"] = processes

    @property
    def feature_names(self) -> List[str]:
        return [
            "ccs_observed_im2deep",
            "ccs_predicted_im2deep",
            "ccs_error_im2deep",
            "abs_ccs_error_im2deep",
            "perc_ccs_error_im2deep",
        ]

    def add_features(self, psm_list: PSMList) -> None:
        """Add IM2Deep-derived features to PSMs"""
        logger.info("Adding IM2Deep-derived features to PSMs")

        # Convert ion mobility to CCS
        psm_list_df = psm_list.to_dataframe()
        # Remove unnecessary columns to save memory #TODO: optimize further?
        psm_list_df = psm_list_df[
            [
                "peptidoform",
                "ion_mobility",
                "precursor_mz",
                "run",
                "qvalue",
                "is_decoy",
                "metadata",
            ]
        ]

        psm_list_df["charge"] = [pep.precursor_charge for pep in psm_list_df["peptidoform"]]
        psm_list_df["sequence"] = psm_list_df["peptidoform"].apply(lambda x: x.proforma)
        psm_list_df["ccs_observed_im2deep"] = im2ccs(
            psm_list_df["ion_mobility"],
            psm_list_df["precursor_mz"],
            psm_list_df["charge"],
        )

        # Make predictions with IM2Deep
        logger.debug("Predicting CCS values...")
        predictions = predict_ccs(psm_list, write_output=False, **self.im2deep_kwargs)
        psm_list_df["ccs_predicted"] = predictions

        # Create dataframe with high confidence hits for calibration
        logger.debug("Calibrating IM2Deep...")
        reference_dataset = pd.read_csv(REFERENCE_DATASET_PATH)
        reference_dataset["charge"] = reference_dataset["peptidoform"].apply(
            lambda x: int(x.split("/")[1]) if isinstance(x, str) else x.precursor_charge
        )
        logger.debug(f"Loaded reference dataset with {len(reference_dataset)} entries")

        run_shift_dict = {}
        for run in psm_list_df["run"].unique():
            cal_run_psm_df = self.make_calibration_df(psm_list_df[psm_list_df["run"] == run])
            # Rename for calculate_ccs_shift compatibility
            cal_run_psm_df = cal_run_psm_df.rename(
                columns={"ccs_observed_im2deep": "ccs_observed"}
            )
            shift = calculate_ccs_shift(
                cal_df=cal_run_psm_df, reference_dataset=reference_dataset, per_charge=True
            )
            run_shift_dict[run] = shift
        shift_df = pd.DataFrame.from_dict(run_shift_dict, orient="index").stack().reset_index()
        shift_df.columns = ["run", "charge", "ccs_shift"]

        # Apply calibration shifts
        psm_list_df = psm_list_df.merge(shift_df, on=["run", "charge"], how="left")
        psm_list_df["ccs_shift"] = psm_list_df["ccs_shift"].fillna(
            0
        )  # Fill missing shifts with 0 (no calibration for that run/charge) #TODO check with ROBBE
        psm_list_df["ccs_predicted_im2deep"] = (
            psm_list_df["ccs_predicted"] + psm_list_df["ccs_shift"]
        )
        psm_list_df["ccs_error_im2deep"] = (
            psm_list_df["ccs_predicted_im2deep"] - psm_list_df["ccs_observed_im2deep"]
        )
        psm_list_df["abs_ccs_error_im2deep"] = np.abs(psm_list_df["ccs_error_im2deep"])
        psm_list_df["perc_ccs_error_im2deep"] = (
            np.abs(psm_list_df["ccs_error_im2deep"]) / psm_list_df["ccs_observed_im2deep"] * 100
        )

        psm_list_feature_dicts = psm_list_df[self.feature_names].to_dict(orient="records")
        # Add features to PSMs
        logger.debug("Adding features to PSMs...")
        for psm, features in zip(psm_list, psm_list_feature_dicts):
            psm.rescoring_features.update(features)

    @staticmethod
    def make_calibration_df(psm_list_df: pd.DataFrame, threshold: float = 0.25) -> pd.DataFrame:
        """
        Make dataframe for calibration of IM2Deep predictions.

        Parameters
        ----------
        psm_list_df
            DataFrame with PSMs.
        threshold
            Percentage of highest scoring identified target PSMs to use for calibration,
            default 0.25.

        Returns
        -------
        pd.DataFrame
            DataFrame with high confidence hits for calibration.

        """
        identified_psms = psm_list_df[
            (psm_list_df["qvalue"] < 0.01)
            & (~psm_list_df["is_decoy"])
            & (psm_list_df["charge"] < 7)  # predictions do not go higher for IM2Deep
            & ([metadata.get("original_psm", True) for metadata in psm_list_df["metadata"]])
        ]
        calibration_psms = identified_psms[
            identified_psms["qvalue"] < identified_psms["qvalue"].quantile(1 - threshold)
        ]
        if isinstance(calibration_psms["peptidoform"].iloc[0], Peptidoform):
            calibration_psms["peptidoform"] = calibration_psms["peptidoform"].apply(
                lambda x: x.proforma
            )
        logger.debug(
            f"Number of high confidence hits for calculating shift: {len(calibration_psms)}"
        )
        return calibration_psms
