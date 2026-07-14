"""In-memory data container for MS²Rescore report generation."""

import json
import logging
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd
import psm_utils
import psm_utils.io

from ms2rescore.feature_generators import FEATURE_GENERATORS
from ms2rescore.report.utils import (
    compute_protein_stats,
    create_psm_dataframe,
    read_feature_names,
)

logger = logging.getLogger(__name__)

# PSM dataframe columns that are not rescoring features
_NON_FEATURE_COLUMNS = {
    "spectrum_id",
    "run",
    "collection",
    "spectrum",
    "peptidoform",
    "precursor_mz",
    "retention_time",
    "ion_mobility",
    "protein_list",
    "rank",
    "source",
    "provenance_data",
    "metadata",
    "rescoring_features",
    "qvalue",
    "pep",
    "score",
    "precursor_charge",
    "is_decoy",
    "score_before",
    "qvalue_before",
    "score_after",
    "qvalue_after",
}


@dataclass
class ReportData:
    """
    Everything needed to render an MS²Rescore report, held in memory.

    Instances are built either from an in-memory pipeline run
    (:py:meth:`from_run`) or from the files written by a previous run
    (:py:meth:`from_files`). The report generator consumes these fields directly
    and performs no file reading or dataframe construction of its own.
    """

    psm_df: pd.DataFrame
    feature_names: Dict[str, List[str]]
    config: dict = field(default_factory=lambda: {"ms2rescore": {}})
    feature_weights: Optional[pd.DataFrame] = None
    protein_stats: Optional[List[dict]] = None
    log_html: Optional[str] = None

    @classmethod
    def from_run(
        cls,
        psm_list: psm_utils.PSMList,
        feature_names: Optional[Dict[str, set]] = None,
        config: Optional[dict] = None,
        feature_weights: Optional[pd.DataFrame] = None,
    ) -> "ReportData":
        """Build report data from an in-memory MS²Rescore run."""
        config = config or {"ms2rescore": {}}
        psm_df = create_psm_dataframe(psm_list)
        feature_names = _normalize_feature_names(feature_names) or _infer_feature_names(psm_df)
        protein_stats = compute_protein_stats(psm_df, _resolve_fasta(config))
        return cls(
            psm_df=psm_df,
            feature_names=feature_names,
            config=config,
            feature_weights=feature_weights,
            protein_stats=protein_stats,
        )

    @classmethod
    def from_files(cls, output_path_prefix: str) -> "ReportData":
        """Build report data by reading the files written by a previous run."""
        psm_file = Path(output_path_prefix + ".psms.tsv")
        if not psm_file.is_file():
            raise FileNotFoundError(f"PSM file not found: {psm_file.as_posix()}")

        logger.info("Reading PSMs from %s...", psm_file.as_posix())
        psm_list = psm_utils.io.read_file(psm_file, filetype="tsv", show_progressbar=True)
        psm_df = create_psm_dataframe(psm_list)

        config = _read_config(Path(output_path_prefix + ".full-config.json"))
        feature_names = _read_feature_names_or_infer(
            Path(output_path_prefix + ".feature_names.tsv"), psm_df
        )
        feature_weights = _read_feature_weights(Path(output_path_prefix + ".mokapot.weights.tsv"))
        protein_stats = compute_protein_stats(psm_df, _resolve_fasta(config))
        return cls(
            psm_df=psm_df,
            feature_names=feature_names,
            config=config,
            feature_weights=feature_weights,
            protein_stats=protein_stats,
        )


def _normalize_feature_names(feature_names: Optional[Dict[str, set]]) -> Dict[str, List[str]]:
    """Convert a generator -> feature-name mapping to plain lists, dropping empties."""
    if not feature_names:
        return {}
    return {gen: list(features) for gen, features in feature_names.items() if features}


def _infer_feature_names(psm_df: pd.DataFrame) -> Dict[str, List[str]]:
    """Infer the generator -> feature-name mapping from the PSM dataframe columns."""
    feature_columns = [col for col in psm_df.columns if col not in _NON_FEATURE_COLUMNS]
    if not feature_columns:
        return {}

    # Map each feature to its generator by instantiating the known generators
    feature_to_generator = {}
    for generator_name, generator_class in FEATURE_GENERATORS.items():
        try:
            for feature in generator_class().feature_names:
                feature_to_generator[feature] = generator_name
        except Exception:
            continue

    feature_names = defaultdict(list)
    for feature in feature_columns:
        feature_names[feature_to_generator.get(feature, "other")].append(feature)
    return dict(feature_names)


def _read_feature_names_or_infer(path: Path, psm_df: pd.DataFrame) -> Dict[str, List[str]]:
    """Read feature names from file, falling back to inference from the dataframe."""
    feature_names = read_feature_names(path)
    if feature_names:
        return {gen: list(features) for gen, features in feature_names.items()}
    logger.info("Feature names file not found. Inferring from feature generators...")
    return _infer_feature_names(psm_df)


def _read_config(path: Path) -> dict:
    """Read the full-config JSON, returning an empty config if it is missing."""
    if not path.is_file():
        logger.info("No configuration file found. Proceeding without it.")
        return {"ms2rescore": {}}
    try:
        return json.loads(path.read_text())
    except (json.JSONDecodeError, OSError):
        logger.warning("Could not read configuration file. Proceeding without it.")
        return {"ms2rescore": {}}


def _read_feature_weights(path: Path) -> Optional[pd.DataFrame]:
    """Read mokapot feature weights (one row per fold), or None if the file is absent."""
    if not path.is_file():
        logger.info("Feature weights file not found. Skipping feature weights.")
        return None
    try:
        return pd.read_csv(path, sep="\t")
    except (OSError, pd.errors.ParserError) as e:
        logger.warning("Could not read feature weights file: %s", e)
        return None


def _resolve_fasta(config: dict) -> Optional[str]:
    """Resolve the fasta path from config, checking the known locations in order."""
    ms2rescore_config = config.get("ms2rescore", {})
    engine_config = ms2rescore_config.get("rescoring_engine", {}).get("mokapot", {})
    return (
        ms2rescore_config.get("fasta_file")
        or config.get("fasta_file")
        or engine_config.get("fasta_file")
    )
