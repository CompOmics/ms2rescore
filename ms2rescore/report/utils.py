"""Utility functions for MS²Rescore report generation."""

import logging
from collections import defaultdict
from csv import DictReader
from pathlib import Path
from typing import List, Optional

import pandas as pd
import psm_utils
from ristretto import RescoreResult

logger = logging.getLogger(__name__)

# Stat card background color, one per identification level
_STAT_CARD_COLORS = {
    "psms": "card-bg-blue",
    "peptides": "card-bg-green",
    "proteins": "card-bg-red",
}

_FDR_THRESHOLD = 0.01


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


def _n_identified(df: pd.DataFrame, fdr_threshold: float) -> int:
    """
    Count target rows passing the FDR threshold (decoys are never identifications).

    Mumble-generated alternate candidates are already excluded from ``before`` upstream, in
    :py:func:`ms2rescore.rescoring.evaluate_before`, so no masking is needed here.

    """
    return int(((df["qvalue"] <= fdr_threshold) & ~df["is_decoy"]).sum())


def compute_protein_stats(
    before: RescoreResult, after: RescoreResult, fdr_threshold: float = _FDR_THRESHOLD
) -> Optional[List[dict]]:
    """
    Compare protein-group-level identifications before and after rescoring.

    Returns a single-element list with a stat card, or ``None`` if either `RescoreResult` has no
    protein-level rollup (i.e. no ``protein_col`` was available) or no proteins pass the FDR
    threshold.
    """
    if before.proteins is None or after.proteins is None:
        return None

    n_before = _n_identified(before.proteins, fdr_threshold)
    n_after = _n_identified(after.proteins, fdr_threshold)
    if n_before == 0 or n_after == 0:
        return None

    return [build_stat_card("Protein groups", "proteins", n_before, n_after)]


def compute_id_stats(
    before: RescoreResult, after: RescoreResult, fdr_threshold: float = _FDR_THRESHOLD
) -> List[dict]:
    """Build the PSM/peptide/(optional) protein overview stat cards from before/after results."""
    stats = []

    n_before_psms = _n_identified(before.psms, fdr_threshold)
    n_after_psms = _n_identified(after.psms, fdr_threshold)
    if n_before_psms > 0:
        stats.append(build_stat_card("PSMs", "psms", n_before_psms, n_after_psms))

    n_before_peptides = _n_identified(before.peptidoforms, fdr_threshold)
    n_after_peptides = _n_identified(after.peptidoforms, fdr_threshold)
    if n_before_peptides > 0:
        stats.append(build_stat_card("Peptides", "peptides", n_before_peptides, n_after_peptides))

    protein_stats = compute_protein_stats(before, after, fdr_threshold)
    if protein_stats:
        stats.extend(protein_stats)

    return stats


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
    Create a comprehensive dataframe from a PSM list with all information needed for the report.

    Includes basic PSM information (peptidoform, score, qvalue, is_decoy, rank, etc.) and all
    rescoring features. Rows are exactly ``psm_list``'s rows -- may include multiple ranks per
    spectrum if ``max_psm_rank_output > 1`` was configured; charts that need a single row per
    spectrum (or a before/after comparison) collapse or look it up themselves, from the
    before/after ``RescoreResult``s directly, rather than from a merge performed here.

    """
    psm_df = psm_list.to_dataframe()

    # Add rescoring features - vectorized extraction
    if psm_list[0].rescoring_features:
        features_df = pd.DataFrame.from_records(psm_list["rescoring_features"]).astype("float32")
        psm_df = pd.concat([psm_df.reset_index(drop=True), features_df], axis=1)
        # Remove duplicate columns (keep last, i.e., from features_df)
        psm_df = psm_df.loc[:, ~psm_df.columns.duplicated(keep="last")]

    return psm_df
