"""
DeepLC retention time-based feature generator.

DeepLC is a fully modification-aware peptide retention time predictor. It uses a deep convolutional
neural network to predict retention times based on the atomic composition of the (modified) amino
acid residues in the peptide. See
`github.com/compomics/deeplc <https://github.com/compomics/deeplc>`_ for more information.

If you use DeepLC through MS²Rescore, please cite:

.. epigraph::
    Bouwmeester, R., Gabriels, R., Hulstaert, N. et al. DeepLC can predict retention times for
    peptides that carry unknown modifications. *Nat Methods* 18, 1363-1369 (2021).
    `doi:10.1038/s41592-021-01301-5 <https://doi.org/10.1038/s41592-021-01301-5>`_

"""

import logging
import warnings
from typing import List, Union

import numpy as np
from deeplc.calibration import SplineTransformerCalibration

# NOTE: `_best_correlating_head` is a private DeepLC function. Imported here intentionally to
# reuse DeepLC's multitask head-selection logic. To be replaced once DeepLC exposes a public API.
from deeplc.core import _best_correlating_head, finetune, predict
from psm_utils import PSMList

from ms2rescore.feature_generators.base import FeatureGeneratorBase
from ms2rescore.parse_spectra import MSDataType
from ms2rescore._utils import get_original_hit_mask

logger = logging.getLogger(__name__)

# Suppress verbose onnx2torch logging
logging.getLogger("onnx2torch").setLevel(logging.WARNING)


class DeepLCFeatureGenerator(FeatureGeneratorBase):
    """DeepLC retention time-based feature generator."""

    required_ms_data = {MSDataType.retention_time}

    def __init__(
        self,
        *args,
        calibration_set_size: Union[int, float, None] = None,
        processes: int = 1,
        **kwargs,
    ) -> None:
        """
        Generate DeepLC-based features for rescoring.

        DeepLC retraining is on by default. Add ``deeplc_retrain: False`` as a keyword argument to
        disable retraining.

        Parameters
        ----------

        calibration_set_size: int or float
            Amount of best PSMs to use for DeepLC calibration. If this value is lower
            than the number of available PSMs, all PSMs will be used. (default: 0.15)
        processes: {int, None}
            Number of processes to use in DeepLC. Defaults to 1.
        kwargs: dict
            Additional keyword arguments are passed to DeepLC.

        Attributes
        ----------
        feature_names: list[str]
            Names of the features that will be added to the PSMs.

        """
        super().__init__(*args, **kwargs)

        self.calibration_set_size = calibration_set_size
        self.deeplc_kwargs = kwargs or {}

        self._verbose = logger.getEffectiveLevel() <= logging.DEBUG

        self.model = self.deeplc_kwargs.get("model", None)

        self.calibration = None

        # Prepare DeepLC predict kwargs
        self.predict_kwargs = {
            k: v
            for k, v in self.deeplc_kwargs.items()
            if k in ["device", "batch_size", "num_threads"]
        }  # getfullargspec(predict).args does not work on this outer predict function
        self.predict_kwargs["num_threads"] = processes if processes > 0 else None

        # Prepare DeepLC finetune kwargs
        if "deeplc_retrain" not in self.deeplc_kwargs:
            self.deeplc_kwargs["deeplc_retrain"] = False
            return  # skip the rest of the init if no retraining

        if self.deeplc_kwargs["deeplc_retrain"]:
            self.finetune_kwargs = {
                k: v
                for k, v in self.deeplc_kwargs.items()
                if k
                in [
                    "epochs",
                    "device",
                    "batch_size",
                    "learning_rate",
                    "patience",
                    "trainable_layers",
                    "validation_split",
                ]
            }
            self.finetune_kwargs["num_threads"] = processes if processes > 0 else None

    @property
    def feature_names(self) -> List[str]:
        return [
            "observed_retention_time",
            "predicted_retention_time",
            "rt_diff",
            "observed_retention_time_best",
            "predicted_retention_time_best",
            "rt_diff_best",
        ]

    def add_features(self, psm_list: PSMList) -> None:
        """Add DeepLC-derived features to PSMs."""
        warnings.filterwarnings("ignore", category=UserWarning, module="deeplc._features")

        logger.info("Adding DeepLC-derived features to PSMs.")
        # Mumble-generated candidate PSMs are unconfirmed mass-shift explanations and must never
        # be used to calibrate or fine-tune DeepLC, only the original search engine hits can.
        original_hit_mask = get_original_hit_mask(psm_list)
        psm_list_df = psm_list.to_dataframe()
        psm_list_df["original_psm"] = original_hit_mask
        psm_list_df = psm_list_df[
            [
                "peptidoform",
                "retention_time",
                "run",
                "qvalue",
                "is_decoy",
                "original_psm",
            ]
        ]
        psm_list_df["sequence"] = psm_list_df["peptidoform"].apply(lambda x: x.modified_sequence)

        if self.deeplc_kwargs["deeplc_retrain"]:
            # Filter high-confidence target PSMs once for transfer learning
            target_mask = (
                (psm_list["qvalue"] <= 0.01) & (~psm_list["is_decoy"]) & original_hit_mask
            )
            target_psms = psm_list[target_mask]

            # Determine best run for transfer learning
            best_run = self._best_run_by_shared_proteoforms(
                target_psms["run"],
                target_psms["peptidoform"],
            )

            # Fine-tune on best run
            best_run_psms = target_psms[target_psms["run"] == best_run]
            logger.debug(
                f"Fine-tuning DeepLC on run '{best_run}'... with {len(best_run_psms)} PSMs"
            )
            self.model = finetune(
                best_run_psms,
                model=self.model,
                train_kwargs=self.finetune_kwargs,
            )

        # Ensure a positional index so DataFrame rows map 1:1 to pred_matrix rows
        psm_list_df = psm_list_df.reset_index(drop=True)

        # Predict retention times for all PSMs at once. The default DeepLC model is a multitask
        # model, so `return_matrix=True` yields the full (n_psms, n_heads) prediction matrix. The
        # best-correlating head is selected per run below, rather than defaulting to head 0.
        logger.info("Predicting retention times with DeepLC...")
        pred_matrix = predict(
            psm_list,
            model=self.model,
            predict_kwargs=self.predict_kwargs,
            return_matrix=True,
        )

        # Calibrate predictions per run, selecting the best-correlating head from each run's
        # calibration reference PSMs.
        logger.info("Selecting best head and calibrating predicted retention times per run...")
        psm_list_df = psm_list_df.sort_values("qvalue")
        psm_list_df["predicted_retention_time"] = np.nan
        for run in psm_list_df["run"].unique():
            run_df = psm_list_df[psm_list_df["run"] == run]

            # Get calibration data (target PSMs only): row indices into pred_matrix + observed RTs
            calibration_idx, observed_rt_calibration = self._get_calibration_data(run_df)
            if len(observed_rt_calibration) == 0:
                raise ValueError(f"Run '{run}' has no target PSMs available for calibration.")

            # Select the head that best correlates with observed RT on the calibration PSMs
            reference_matrix = pred_matrix[calibration_idx]
            head = _best_correlating_head(reference_matrix, observed_rt_calibration)
            logger.debug(f"Run '{run}': selected DeepLC head {head}")

            # Fit calibration on the selected head and transform all predictions for this run
            calibration = SplineTransformerCalibration()
            calibration.fit(
                target=observed_rt_calibration,
                source=reference_matrix[:, head],
            )

            run_idx = run_df.index.values
            calibrated_rt = calibration.transform(pred_matrix[run_idx, head])

            # Update predictions with calibrated values
            psm_list_df.loc[psm_list_df["run"] == run, "predicted_retention_time"] = calibrated_rt

        del pred_matrix

        psm_list_df = (
            psm_list_df.sort_index()
        )  # restore original PSM order after sort_values above
        psm_list_df["observed_retention_time"] = psm_list_df["retention_time"]
        psm_list_df["rt_diff"] = (
            psm_list_df["observed_retention_time"] - psm_list_df["predicted_retention_time"]
        ).abs()
        psm_list_df_best = psm_list_df.loc[
            psm_list_df.groupby(["run", "sequence"])["rt_diff"].idxmin()
        ]

        psm_list_df = psm_list_df.merge(
            psm_list_df_best[
                [
                    "run",
                    "sequence",
                    "observed_retention_time",
                    "predicted_retention_time",
                    "rt_diff",
                ]
            ],
            on=["run", "sequence"],
            how="left",
            suffixes=("", "_best"),
        )

        psm_list_feature_dicts = psm_list_df[self.feature_names].to_dict(orient="records")
        # Add features to PSMs
        logger.debug("Adding features to PSMs...")
        for psm, features in zip(psm_list, psm_list_feature_dicts):
            psm.rescoring_features.update(features)

    def _get_calibration_data(self, run_df) -> tuple[np.ndarray, np.ndarray]:
        """Get calibration data (pred_matrix row indices and observed RTs) from run dataframe.

        Only target (non-decoy), original (non-Mumble) PSMs are used for calibration.

        Parameters
        ----------
        run_df : pd.DataFrame
            Dataframe containing PSMs for a single run, pre-sorted by qvalue ascending, with a
            positional index into the prediction matrix and columns: 'retention_time', 'qvalue',
            'is_decoy', 'original_psm'

        Returns
        -------
        tuple[np.ndarray, np.ndarray]
            Row indices (into pred_matrix) and observed retention times for calibration
        """
        # Filter to target, original-hit PSMs only (run_df is pre-sorted by qvalue). Mumble
        # candidates are unconfirmed mass-shift explanations and must not calibrate the model.
        target_df = run_df[~run_df["is_decoy"] & run_df["original_psm"]]

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

        calibration_df = target_df.head(num_calibration_psms)

        return (
            calibration_df.index.values,
            calibration_df["retention_time"].values,
        )

    @staticmethod
    def _best_run_by_shared_proteoforms(runs, proteoforms):
        """
        Return the run whose proteoform set has the largest total overlap with all other runs:
            score(run_i) = sum_{j != i} |P_i ∩ P_j|

        Tie / degenerate handling (per request):
        - If no run shares anything (all scores == 0): warn and return the first run (by first appearance).
        - If multiple runs are tied for best score: return the first among them (by first appearance).
        """
        logger.debug("Determining best run for transfer learning based on shared proteoforms.")
        runs = np.asarray(runs)
        proteoforms = np.asarray(proteoforms)
        if runs.shape[0] != proteoforms.shape[0]:
            raise ValueError("runs and proteoforms must have the same length.")
        if runs.size == 0:
            raise ValueError("Empty input: runs/proteoforms must contain at least one entry.")

        # Preserve run order by first appearance
        run_to_idx = {}
        run_order = []
        run_idx = np.empty(runs.shape[0], dtype=np.int64)
        for i, r in enumerate(runs):
            if r not in run_to_idx:
                run_to_idx[r] = len(run_order)
                run_order.append(r)
            run_idx[i] = run_to_idx[r]

        # Fast path: sparse incidence matrix
        try:
            from scipy.sparse import coo_matrix
        except ImportError:
            # Fallback (slower): set-based
            run_sets = {}
            for r, p in zip(runs, proteoforms):
                run_sets.setdefault(r, set()).add(p)

            scores = []
            for r in run_order:
                Pi = run_sets[r]
                s = 0
                for r2 in run_order:
                    if r2 is r:
                        continue
                    s += len(Pi & run_sets[r2])
                scores.append(s)

            scores = np.asarray(scores, dtype=np.int64)
            max_score = scores.max()
            best_candidates = np.flatnonzero(scores == max_score)

            if max_score == 0:
                logger.warning(
                    "No runs share any identified proteoforms with other runs; transfer learning "
                    "might not be as effective."
                )
                return run_order[0]

            return run_order[int(best_candidates[0])]

        # Encode proteoforms (order does not matter for correctness)
        _, prot_idx = np.unique(proteoforms, return_inverse=True)

        # De-duplicate (run, proteoform) pairs
        pairs = np.unique(np.stack([run_idx, prot_idx], axis=1), axis=0)
        r = pairs[:, 0]
        p = pairs[:, 1]

        M = coo_matrix(
            (np.ones(len(r), dtype=np.int32), (r, p)),
            shape=(len(run_order), int(prot_idx.max()) + 1),
        ).tocsr()

        # Overlap matrix O[i,j] = |P_i ∩ P_j|
        overlap = (M @ M.T).toarray()
        np.fill_diagonal(overlap, 0)

        scores = overlap.sum(axis=1).astype(np.int64)
        max_score = int(scores.max())
        best_candidates = np.flatnonzero(scores == max_score)

        if max_score == 0:
            logger.warning(
                "No runs share any identified proteoforms with other runs; transfer learning "
                "might not be as effective."
            )
            return run_order[0]

        return run_order[int(best_candidates[0])]
