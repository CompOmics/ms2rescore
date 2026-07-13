"""
IM2Deep ion mobility-based feature generator.

IM2Deep is a fully modification-aware peptide ion mobility predictor. It uses a deep convolutional
neural network to predict retention times based on the atomic composition of the (modified) amino
acid residues in the peptide. See
`github.com/compomics/IM2Deep <https://github.com/compomics/IM2Deep>`_ for more information.

"""

import logging
from typing import List, Union

import numpy as np
import pandas as pd
from im2deep.calibration import LinearCCSCalibration, get_default_reference
from im2deep.core import predict
from im2deep.utils import im2ccs
from psm_utils import PSMList

from ms2rescore.feature_generators.base import FeatureGeneratorBase
from ms2rescore.parse_spectra import MSDataType

logger = logging.getLogger(__name__)


class IM2DeepFeatureGenerator(FeatureGeneratorBase):
    """IM2Deep collision cross section feature generator."""

    required_ms_data = {MSDataType.ion_mobility}

    def __init__(
        self,
        multi: bool = False,
        calibration_set_size: Union[int, float] = None,
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

        self.multi = multi
        if self.multi:
            raise NotImplementedError(
                "Multi-IM mode is not yet implemented for IM2DeepFeatureGenerator."
            )

        self.calibration_set_size = calibration_set_size

        self.im2deep_kwargs = kwargs or {}
        self._verbose = logger.getEffectiveLevel() <= logging.DEBUG

        self.model = self.im2deep_kwargs.get("model", None)

        # Prepare IM2Deep predict kwargs
        self.predict_kwargs = {
            k: v for k, v in self.im2deep_kwargs.items() if k in ["device", "batch_size"]
        }
        self.predict_kwargs["num_workers"] = processes

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

        psm_list_df = psm_list.to_dataframe()
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

        psm_list_df["sequence"] = psm_list_df["peptidoform"].apply(lambda x: x.modified_sequence)
        psm_list_df["charge"] = [pep.precursor_charge for pep in psm_list_df["peptidoform"]]
        psm_list_df["ccs_observed_im2deep"] = im2ccs(
            psm_list_df["ion_mobility"],
            psm_list_df["precursor_mz"],
            psm_list_df["charge"],
        )

        # Make predictions with IM2Deep
        logger.info("Predicting CCS values with IM2Deep...")
        # float64 so the per-run calibrated write-back below stays dtype-consistent (predict
        # returns float32; calibration.transform returns object).
        psm_list_df["predicted_CCS_uncalibrated"] = np.asarray(
            predict(psm_list, model=self.model, predict_kwargs=self.predict_kwargs),
            dtype="float64",
        )

        # getting reference CCS values for calibration
        source_dataframe = get_default_reference(multi=self.multi)

        # Create dataframe with high confidence hits for calibration
        logger.info("Calibrating predicted CCS values per run...")
        for run in psm_list_df["run"].unique():
            run_df = psm_list_df[psm_list_df["run"] == run].copy()

            calibration_df = self._get_im_calibration_data(run_df)
            if calibration_df.empty:
                raise ValueError(f"Run '{run}' has no target PSMs available for calibration.")

            calibration = LinearCCSCalibration()
            calibration.fit(
                psm_df_target=calibration_df,
                psm_df_source=source_dataframe,
            )

            calibrated_im = calibration.transform(
                run_df[["peptidoform", "predicted_CCS_uncalibrated"]]
            )

            # Update predictions with calibrated values. transform() returns an object-dtype
            # array; cast to float64 so the write-back does not upcast the column (pandas
            # incompatible-dtype FutureWarning).
            psm_list_df.loc[psm_list_df["run"] == run, "predicted_CCS_uncalibrated"] = np.asarray(
                calibrated_im, dtype="float64"
            )

        # Apply calibration shifts
        psm_list_df.rename(
            columns={"predicted_CCS_uncalibrated": "ccs_predicted_im2deep"}, inplace=True
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

    def _get_im_calibration_data(self, run_df) -> pd.DataFrame:
        """Get calibration data (observed and predicted CCS values) from run dataframe.

        Only target (non-decoy) PSMs are used for calibration.

        Parameters
        ----------
        run_df : pd.DataFrame
            Dataframe containing PSMs for a single run, with columns:
            'ccs_observed_im2deep', 'ccs_predicted_im2deep', 'qvalue', 'is_decoy'

        Returns
        -------
        pd.DataFrame
            DataFrame with 'peptidoform' and 'CCS' columns for calibration.
        """
        # Filter to target PSMs only
        target_df = run_df[~run_df["is_decoy"]].copy()
        target_df = target_df.sort_values("qvalue", ascending=True)

        # Determine number of calibration PSMs
        if isinstance(self.calibration_set_size, float):
            if not 0 < self.calibration_set_size <= 1:
                raise ValueError(
                    "If `calibration_set_size` is a float, it cannot be smaller than "
                    "or equal to 0 or larger than 1."
                )
            num_calibration_psms = round(len(target_df) * self.calibration_set_size)
        elif isinstance(self.calibration_set_size, int):
            if self.calibration_set_size > len(target_df):
                logger.warning(
                    f"Requested number of calibration PSMs ({self.calibration_set_size}) "
                    f"is larger than total number of target PSMs ({len(target_df)}). Using "
                    "all target PSMs for calibration."
                )
                num_calibration_psms = len(target_df)
            else:
                num_calibration_psms = self.calibration_set_size
        else:
            # Use PSMs with q-value <= 0.01
            num_calibration_psms = (target_df["qvalue"] <= 0.01).sum()

        logger.debug(f"Using {num_calibration_psms} target PSMs for calibration")

        # Select calibration PSMs (assuming they are sorted by q-value)
        calibration_df = target_df.head(num_calibration_psms)

        return calibration_df[["peptidoform", "ccs_observed_im2deep"]].rename(
            columns={"ccs_observed_im2deep": "CCS"}
        )
