"""Utility functions for MS²Rescore report generation."""

import logging
from collections import defaultdict
from csv import DictReader
from pathlib import Path
from typing import List, Optional

import pandas as pd
import psm_utils

from ms2rescore.constants import CHARGE_PATTERN
from mokapot import LinearPsmDataset, read_fasta

logger = logging.getLogger(__name__)

# Stat card background color, one per identification level
_STAT_CARD_COLORS = {"psms": "card-bg-blue", "peptides": "card-bg-green", "proteins": "card-bg-red"}


def read_feature_names(feature_names_path: Optional[Path]) -> dict:
    """Read feature names and mapping with feature generator from file."""
    feature_names = defaultdict(list)
    if not feature_names_path or not feature_names_path.is_file():
        return feature_names

    try:
        with open(feature_names_path) as f:
            reader = DictReader(f, delimiter="\t")
            for line in reader:
                feature_names[line["feature_generator"]].append(line["feature_name"])
    except (FileNotFoundError, KeyError, ValueError) as e:
        logger.warning(f"Could not read feature names file: {e}")

    return feature_names


def compute_protein_stats(psm_df: pd.DataFrame, fasta_file: Optional[str]) -> Optional[List[dict]]:
    """
    Compute protein-group-level statistics before and after rescoring.

    This is the only report function that relies on mokapot's protein inference. It builds a
    :py:class:`~mokapot.dataset.LinearPsmDataset`, adds proteins from ``fasta_file``, and assigns
    confidence on the before- and after-rescoring scores to count accepted protein groups at 1%
    FDR. Returns a single-element list with a stat card, or ``None`` if no fasta is given or
    mokapot fails to assign confidence.
    """
    if not fasta_file:
        return None
    if "score_before" not in psm_df.columns or "score_after" not in psm_df.columns:
        logger.warning("Before/after scores not found. Skipping protein-group statistics.")
        return None

    # Build mokapot dataset with protein inference from the fasta
    peptide = psm_df["peptidoform"].astype(str).str.replace(CHARGE_PATTERN, "", n=1, regex=True)
    psms = pd.DataFrame({"peptide": peptide, "is_target": ~psm_df["is_decoy"]}).reset_index()
    dataset = LinearPsmDataset(
        psms=psms,
        target_column="is_target",
        spectrum_columns="index",
        peptide_column="peptide",
    )
    try:
        dataset.add_proteins(read_fasta(fasta_file))
    except (FileNotFoundError, ValueError) as e:
        logger.warning("Could not add proteins from fasta: %s. Skipping protein statistics.", e)
        return None

    # Count accepted protein groups on before- and after-rescoring scores
    counts = {}
    for when, score_column in [("before", "score_before"), ("after", "score_after")]:
        try:
            confidence = dataset.assign_confidence(scores=list(psm_df[score_column].astype(float)))
            counts[when] = confidence.accepted["proteins"]
        except (RuntimeError, IndexError, KeyError):
            logger.warning("Could not assign protein confidence for %s rescoring.", when)
            counts[when] = None

    before, after = counts["before"], counts["after"]
    if not before or not after:
        return None

    return [build_stat_card("Protein groups", "proteins", before, after)]


def build_stat_card(item: str, level: str, before: int, after: int) -> dict:
    """Build a single overview stat card comparing before/after counts."""
    increase = (after - before) / before * 100
    return {
        "item": item,
        "card_color": _STAT_CARD_COLORS[level],
        "number": after,
        "diff": f"({after - before:+})",
        "percentage": f"{increase:.1f}%",
        "is_increase": increase > 0,
        "bar_percentage": before / after * 100 if increase > 0 else after / before * 100,
        "bar_color": "#24a143" if increase > 0 else "#a12424",
    }


def create_psm_dataframe(psm_list: psm_utils.PSMList) -> pd.DataFrame:
    """
    Create a comprehensive dataframe from PSM list with all necessary information.

    This dataframe includes:
    - Basic PSM information (peptide, score, qvalue, is_decoy, etc.)
    - Before rescoring scores from provenance data
    - All rescoring features

    Parameters
    ----------
    psm_list
        PSM list to convert to dataframe.

    Returns
    -------
    pd.DataFrame
        Dataframe with all PSM information.
    """
    # Start with basic PSM dataframe
    psm_df = psm_list.to_dataframe()

    # Add before rescoring scores from provenance data
    try:
        provenance_df = pd.DataFrame.from_records(psm_list["provenance_data"])
        if "before_rescoring_score" in provenance_df.columns:
            psm_df["score_before"] = provenance_df["before_rescoring_score"].astype(float)
        if "before_rescoring_qvalue" in provenance_df.columns:
            psm_df["qvalue_before"] = provenance_df["before_rescoring_qvalue"].astype(float)
    except (KeyError, ValueError) as e:
        logger.warning("Could not extract before rescoring scores from provenance data: %s", e)
        psm_df["score_before"] = None
        psm_df["qvalue_before"] = None

    # Add rescoring features - vectorized extraction
    if psm_list[0].rescoring_features:
        # Extract all rescoring_features dicts at once (much faster than looping)
        features_df = pd.DataFrame.from_records(psm_list["rescoring_features"]).astype("float32")
        # Merge features with PSM dataframe (they should have same index)
        psm_df = pd.concat([psm_df, features_df], axis=1)
        # Remove duplicate columns (keep last, i.e., from features_df)
        psm_df = psm_df.loc[:, ~psm_df.columns.duplicated(keep="last")]

    # Rename current score/qvalue to score_after/qvalue_after for clarity
    psm_df["score_after"] = psm_df["score"]
    psm_df["qvalue_after"] = psm_df["qvalue"]

    return psm_df
