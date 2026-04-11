"""
Evidence score for targeted MALDI-MSI acquisitions.

When PSM counts are too low for classifier-based rescoring (e.g., iprm-PASEF
with ~10 targeted precursors), this module provides an interpretable score
combining:
  (a) CCS agreement: predicted vs observed, z-scored against a calibrated
      error distribution
  (b) Fragment-ion intensity similarity: spectral angle between observed and
      MS²PIP-predicted MS/MS spectra

Each PSM gets a composite score and a breakdown of contributing evidence,
making the score transparent and auditable.
"""

import logging
from dataclasses import dataclass
from typing import Optional

import numpy as np
import pandas as pd
from psm_utils import PSMList

logger = logging.getLogger(__name__)


@dataclass
class EvidenceBreakdown:
    """Per-PSM evidence breakdown."""

    spectrum_id: str
    peptidoform: str
    composite_score: float

    # CCS evidence
    ccs_observed: Optional[float]
    ccs_predicted: Optional[float]
    ccs_error: Optional[float]
    ccs_z_score: Optional[float]
    ccs_score: float

    # MS/MS evidence
    spectral_angle: Optional[float]
    ms2_score: float

    # Weights used
    ccs_weight: float
    ms2_weight: float


class EvidenceScorer:
    """
    Compute interpretable evidence scores for targeted MALDI-MSI PSMs.

    Parameters
    ----------
    ccs_weight : float
        Weight for CCS agreement in composite score. Default 0.6.
    ms2_weight : float
        Weight for fragment similarity in composite score. Default 0.4.
    ccs_error_mean : float or None
        Mean of the calibrated CCS error distribution (Angstrom^2).
        If None, estimated from the provided data.
    ccs_error_std : float or None
        Std of the calibrated CCS error distribution.
        If None, estimated from the provided data.
    """

    def __init__(
        self,
        ccs_weight: float = 0.6,
        ms2_weight: float = 0.4,
        ccs_error_mean: Optional[float] = None,
        ccs_error_std: Optional[float] = None,
    ):
        if abs(ccs_weight + ms2_weight - 1.0) > 1e-6:
            raise ValueError("ccs_weight + ms2_weight must equal 1.0")

        self.ccs_weight = ccs_weight
        self.ms2_weight = ms2_weight
        self.ccs_error_mean = ccs_error_mean
        self.ccs_error_std = ccs_error_std

    def score_psms(
        self,
        psm_list: PSMList,
        ccs_observed: np.ndarray,
        ccs_predicted: np.ndarray,
        spectral_angles: Optional[np.ndarray] = None,
    ) -> list[EvidenceBreakdown]:
        """
        Compute evidence scores for a list of PSMs.

        Parameters
        ----------
        psm_list
            PSMs to score.
        ccs_observed
            Observed CCS values (from MALDI 1/K0 conversion).
        ccs_predicted
            Predicted CCS values (from calibrated IM2Deep).
        spectral_angles
            Spectral angle similarity scores (0-1, from MS²PIP comparison).
            If None, MS/MS evidence is set to 0 and CCS gets full weight.

        Returns
        -------
        List of EvidenceBreakdown objects, one per PSM.
        """
        n = len(psm_list)

        if len(ccs_observed) != n or len(ccs_predicted) != n:
            raise ValueError("Array lengths must match PSM list length")

        # Compute CCS errors
        ccs_errors = ccs_predicted - ccs_observed

        # Estimate error distribution if not provided
        if self.ccs_error_mean is None:
            self.ccs_error_mean = float(np.nanmean(ccs_errors))
        if self.ccs_error_std is None:
            self.ccs_error_std = float(np.nanstd(ccs_errors))
            if self.ccs_error_std == 0:
                self.ccs_error_std = 1.0  # Prevent division by zero

        # Z-score the CCS errors
        ccs_z_scores = (ccs_errors - self.ccs_error_mean) / self.ccs_error_std

        # Convert z-scores to [0, 1] scores using survival function
        # Low |z| → high score (good agreement)
        ccs_scores = np.exp(-0.5 * ccs_z_scores**2)

        # Handle MS/MS evidence
        if spectral_angles is not None:
            ms2_scores = np.clip(spectral_angles, 0, 1)
            effective_ccs_weight = self.ccs_weight
            effective_ms2_weight = self.ms2_weight
        else:
            ms2_scores = np.zeros(n)
            effective_ccs_weight = 1.0
            effective_ms2_weight = 0.0

        # Composite score
        composite = (
            effective_ccs_weight * ccs_scores + effective_ms2_weight * ms2_scores
        )

        # Build breakdowns
        breakdowns = []
        for i, psm in enumerate(psm_list):
            breakdowns.append(EvidenceBreakdown(
                spectrum_id=psm.spectrum_id,
                peptidoform=str(psm.peptidoform),
                composite_score=float(composite[i]),
                ccs_observed=float(ccs_observed[i]) if not np.isnan(ccs_observed[i]) else None,
                ccs_predicted=float(ccs_predicted[i]) if not np.isnan(ccs_predicted[i]) else None,
                ccs_error=float(ccs_errors[i]) if not np.isnan(ccs_errors[i]) else None,
                ccs_z_score=float(ccs_z_scores[i]) if not np.isnan(ccs_z_scores[i]) else None,
                ccs_score=float(ccs_scores[i]),
                spectral_angle=float(ms2_scores[i]) if spectral_angles is not None else None,
                ms2_score=float(ms2_scores[i]),
                ccs_weight=effective_ccs_weight,
                ms2_weight=effective_ms2_weight,
            ))

        logger.info(
            f"Scored {n} PSMs: mean composite={np.mean(composite):.3f}, "
            f"mean CCS score={np.mean(ccs_scores):.3f}"
        )

        return breakdowns

    @staticmethod
    def breakdowns_to_dataframe(breakdowns: list[EvidenceBreakdown]) -> pd.DataFrame:
        """Convert evidence breakdowns to a DataFrame."""
        records = []
        for b in breakdowns:
            records.append({
                "spectrum_id": b.spectrum_id,
                "peptidoform": b.peptidoform,
                "composite_score": b.composite_score,
                "ccs_observed": b.ccs_observed,
                "ccs_predicted": b.ccs_predicted,
                "ccs_error": b.ccs_error,
                "ccs_z_score": b.ccs_z_score,
                "ccs_score": b.ccs_score,
                "spectral_angle": b.spectral_angle,
                "ms2_score": b.ms2_score,
                "ccs_weight": b.ccs_weight,
                "ms2_weight": b.ms2_weight,
            })
        return pd.DataFrame(records)

    @staticmethod
    def format_report(breakdowns: list[EvidenceBreakdown]) -> str:
        """Format a human-readable evidence report."""
        lines = [
            "Evidence Score Report",
            "=" * 60,
            "",
            f"{'Peptidoform':<30s} {'Composite':>10s} {'CCS':>6s} {'MS2':>6s} "
            f"{'CCS err':>8s} {'z':>6s}",
            "-" * 70,
        ]
        for b in sorted(breakdowns, key=lambda x: x.composite_score, reverse=True):
            pf = b.peptidoform[:29]
            ccs_err = f"{b.ccs_error:.1f}" if b.ccs_error is not None else "N/A"
            z = f"{b.ccs_z_score:.2f}" if b.ccs_z_score is not None else "N/A"
            sa = f"{b.spectral_angle:.3f}" if b.spectral_angle is not None else "N/A"
            lines.append(
                f"{pf:<30s} {b.composite_score:>10.4f} {b.ccs_score:>6.3f} "
                f"{sa:>6s} {ccs_err:>8s} {z:>6s}"
            )
        lines.extend([
            "",
            f"CCS error distribution: mean={breakdowns[0].ccs_weight:.2f} weight",
            f"MS/MS weight: {breakdowns[0].ms2_weight:.2f}",
        ])
        return "\n".join(lines)
