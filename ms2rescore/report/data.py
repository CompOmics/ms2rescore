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
from ristretto import RescoreResult

from ms2rescore.feature_generators import FEATURE_GENERATORS
from ms2rescore.report.utils import (
    compute_id_stats,
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
    "pep_before",
    "score_after",
    "qvalue_after",
    "pep_after",
}

_EMPTY_PSMS_COLUMNS = ["spectrum_id", "run", "is_decoy", "peptidoform", "score", "qvalue", "pep"]
_EMPTY_ROLLUP_COLUMNS = ["peptidoform", "score", "qvalue", "pep", "is_decoy", "n_psms"]


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
    id_stats: List[dict] = field(default_factory=list)
    log_html: Optional[str] = None
    before: RescoreResult = field(default_factory=lambda: _empty_result())
    after: RescoreResult = field(default_factory=lambda: _empty_result())
    rescoring_tables_unavailable: bool = False

    @classmethod
    def from_run(
        cls,
        psm_list: psm_utils.PSMList,
        feature_names: Optional[Dict[str, set]] = None,
        config: Optional[dict] = None,
        before: Optional[RescoreResult] = None,
        after: Optional[RescoreResult] = None,
    ) -> "ReportData":
        """Build report data from an in-memory MS²Rescore run."""
        config = config or {"ms2rescore": {}}
        rescoring_tables_unavailable = before is None or after is None
        before = before or _empty_result()
        after = after or _empty_result()
        psm_df = create_psm_dataframe(psm_list, before, after)
        feature_names = _normalize_feature_names(feature_names) or _infer_feature_names(psm_df)
        return cls(
            psm_df=psm_df,
            feature_names=feature_names,
            config=config,
            feature_weights=after.feature_weights,
            id_stats=compute_id_stats(before, after),
            before=before,
            after=after,
            rescoring_tables_unavailable=rescoring_tables_unavailable,
        )

    @classmethod
    def from_files(cls, output_path_prefix: str) -> "ReportData":
        """Build report data by reading the files written by a previous run."""
        psm_file = Path(output_path_prefix + ".psms.tsv")
        if not psm_file.is_file():
            raise FileNotFoundError(f"PSM file not found: {psm_file.as_posix()}")

        logger.info("Reading PSMs from %s...", psm_file.as_posix())
        psm_list = psm_utils.io.read_file(psm_file, filetype="tsv", show_progressbar=True)

        # The report is always a rank-1 view: prefer the separate after_report tables (written
        # when max_psm_rank_output > 1), falling back to after (already rank-1 when it was ==1).
        after_report_path = Path(f"{output_path_prefix}.ristretto.psms_after_report.parquet")
        after_suffix = "after_report" if after_report_path.is_file() else "after"

        before_raw = _read_rescore_result(output_path_prefix, "before")
        after_raw = _read_rescore_result(output_path_prefix, after_suffix)
        rescoring_tables_unavailable = before_raw is None or after_raw is None
        if rescoring_tables_unavailable:
            logger.warning(
                "Ristretto before/after tables not found for '%s'. They were likely disabled "
                "(write_rescoring_tables=false) during the original run. Before/after "
                "identification statistics and comparison charts will be unavailable.",
                output_path_prefix,
            )
        before = before_raw or _empty_result()
        after = after_raw or _empty_result()
        psm_df = create_psm_dataframe(psm_list, before, after)

        config = _read_config(Path(output_path_prefix + ".full-config.json"))
        feature_names = _read_feature_names_or_infer(
            Path(output_path_prefix + ".feature_names.tsv"), psm_df
        )
        return cls(
            psm_df=psm_df,
            feature_names=feature_names,
            config=config,
            feature_weights=after.feature_weights,
            id_stats=compute_id_stats(before, after),
            before=before,
            after=after,
            rescoring_tables_unavailable=rescoring_tables_unavailable,
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


def _empty_result() -> RescoreResult:
    """A placeholder `RescoreResult` used when a before/after result is unavailable."""
    return RescoreResult(
        psms=pd.DataFrame(columns=_EMPTY_PSMS_COLUMNS),
        peptidoforms=pd.DataFrame(columns=_EMPTY_ROLLUP_COLUMNS),
        peptides=None,
        proteins=None,
        pi0=float("nan"),
        n_iterations=[],
        feature_weights=pd.DataFrame(),
    )


def _read_rescore_result(output_path_prefix: str, suffix: str) -> Optional[RescoreResult]:
    """
    Reconstruct a `RescoreResult` from the Parquet tables written by `write_rescoring_tables`.

    Returns ``None`` (rather than a placeholder) when the tables are missing, so callers can
    tell "genuinely unavailable" apart from "loaded, but empty".

    """
    psms_path = Path(f"{output_path_prefix}.ristretto.psms_{suffix}.parquet")
    if not psms_path.is_file():
        return None

    peptidoforms_path = Path(f"{output_path_prefix}.ristretto.peptidoforms_{suffix}.parquet")
    peptides_path = Path(f"{output_path_prefix}.ristretto.peptides_{suffix}.parquet")
    proteins_path = Path(f"{output_path_prefix}.ristretto.proteins_{suffix}.parquet")
    weights_path = Path(f"{output_path_prefix}.ristretto.weights.parquet")

    return RescoreResult(
        psms=pd.read_parquet(psms_path),
        peptidoforms=(
            pd.read_parquet(peptidoforms_path)
            if peptidoforms_path.is_file()
            else pd.DataFrame(columns=_EMPTY_ROLLUP_COLUMNS)
        ),
        peptides=pd.read_parquet(peptides_path) if peptides_path.is_file() else None,
        proteins=pd.read_parquet(proteins_path) if proteins_path.is_file() else None,
        pi0=float("nan"),
        n_iterations=[],
        feature_weights=pd.read_parquet(weights_path)
        if weights_path.is_file()
        else pd.DataFrame(),
    )
