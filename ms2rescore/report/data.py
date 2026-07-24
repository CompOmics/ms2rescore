"""In-memory data container for MS²Rescore report generation."""

import json
import logging
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd
import psm_utils
import psm_utils.io
from ristretto import RescoreResult

from ms2rescore import _ristretto_utils
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
    feature_names: dict[str, list[str]]
    before: RescoreResult
    after: RescoreResult
    config: dict = field(default_factory=lambda: {"ms2rescore": {}})
    feature_weights: pd.DataFrame | None = None
    id_stats: list[dict] = field(default_factory=list)
    log_html: str | None = None
    fdr_threshold: float = 0.01

    @classmethod
    def from_run(
        cls,
        psm_list: psm_utils.PSMList,
        feature_names: dict[str, set] | None = None,
        config: dict | None = None,
        before: RescoreResult | None = None,
        after: RescoreResult | None = None,
        fdr_threshold: float = 0.01,
    ) -> "ReportData":
        """Build report data from an in-memory MS²Rescore run."""
        config = config or {"ms2rescore": {}}
        psm_df = create_psm_dataframe(psm_list)
        feature_names = _normalize_feature_names(feature_names) or _infer_feature_names(psm_df)
        return cls(
            psm_df=psm_df,
            feature_names=feature_names,
            before=before,
            after=after,
            config=config,
            feature_weights=after.feature_weights,
            id_stats=compute_id_stats(before, after, fdr_threshold),
            fdr_threshold=fdr_threshold,
        )

    @classmethod
    def from_files(
        cls, output_path_prefix: str, fdr_threshold: float | None = None
    ) -> "ReportData":
        """
        Build report data by reading the files written by a previous run.

        Only the main PSM list file is read. Before/after ``RescoreResult``s (PSM-level
        scores plus peptidoform/peptide/protein rollups) are fully reconstructed from it --
        the pre-rescoring score is stashed in each PSM's ``provenance_data``, and the
        post-rescoring score/identity is simply the PSM's current state -- so no separate
        rescoring-result table needs to be read.

        ``fdr_threshold``, if given, overrides the ``report_fdr`` stored in the run's
        ``full-config.json`` -- lets a report be regenerated at a different FDR threshold
        without rerunning rescoring.

        """
        psm_file = Path(output_path_prefix + ".tsv")
        if not psm_file.is_file():
            raise FileNotFoundError(f"PSM file not found: {psm_file.as_posix()}")

        config = _read_config(Path(output_path_prefix + ".full-config.json"))
        ms2rescore_config = config.get("ms2rescore", {})
        if fdr_threshold is None:
            fdr_threshold = ms2rescore_config.get("report_fdr", 0.01)

        logger.info("Reading PSMs from %s...", psm_file.as_posix())
        psm_list = psm_utils.io.read_file(psm_file, filetype="tsv", show_progressbar=True)

        before = _ristretto_utils.evaluate_before_from_provenance(psm_list, ms2rescore_config)
        after = _ristretto_utils.evaluate_after_from_psm_list(psm_list, ms2rescore_config)
        psm_df = create_psm_dataframe(psm_list)

        feature_names = _read_feature_names_or_infer(
            Path(output_path_prefix + ".feature_names.tsv"), psm_df
        )
        return cls(
            psm_df=psm_df,
            feature_names=feature_names,
            before=before,
            after=after,
            config=config,
            feature_weights=after.feature_weights,
            id_stats=compute_id_stats(before, after, fdr_threshold),
            fdr_threshold=fdr_threshold,
        )


def _normalize_feature_names(feature_names: dict[str, set] | None) -> dict[str, list[str]]:
    """Convert a generator -> feature-name mapping to plain lists, dropping empties."""
    if not feature_names:
        return {}
    return {gen: list(features) for gen, features in feature_names.items() if features}


def _infer_feature_names(psm_df: pd.DataFrame) -> dict[str, list[str]]:
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
            logger.exception("Could not instantiate feature generator `%s`", generator_name)
            continue

    feature_names = defaultdict(list)
    for feature in feature_columns:
        feature_names[feature_to_generator.get(feature, "other")].append(feature)
    return dict(feature_names)


def _read_feature_names_or_infer(path: Path, psm_df: pd.DataFrame) -> dict[str, list[str]]:
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
