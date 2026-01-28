"""Generate an HTML report with various QC charts for of MS²Rescore results."""

import importlib.resources
import json
import logging
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Dict, Optional

import pandas as pd
import plotly.express as px
import psm_utils.io
from jinja2 import Environment, FileSystemLoader
from plotly.offline import get_plotlyjs_version
from psm_utils.psm_list import PSMList

try:
    import tomllib
except ImportError:
    import tomli as tomllib

import ms2rescore
import ms2rescore.report.charts as charts
import ms2rescore.report.templates as templates
from ms2rescore.feature_generators import FEATURE_GENERATORS
from ms2rescore.report.utils import (
    calculate_fdr_stats,
    create_psm_dataframe,
    get_confidence_estimates,
    get_feature_values,
    infer_feature_names_from_psm_list,
    read_feature_names,
)

logger = logging.getLogger(__name__)

PLOTLY_HTML_KWARGS = {
    "full_html": False,
    "include_plotlyjs": False,
    "include_mathjax": False,
    "config": {
        "displayModeBar": True,
        "displaylogo": False,
    },
}


TEXTS = tomllib.loads(importlib.resources.read_text(templates, "texts.toml"))


def generate_report(
    output_path_prefix: str,
    psm_list: Optional[psm_utils.PSMList] = None,
    feature_names: Optional[Dict[str, set]] = None,
    use_txt_log: bool = False,
    output_file: Optional[Path] = None,
    use_mokapot: bool = False,
):
    """
    Generate the report.

    Parameters
    ----------
    output_path_prefix
        Prefix of the MS²Rescore output file names. For example, if the output PSM file is
        ``/path/to/file.ms2rescore.psms.tsv``, the prefix is ``/path/to/file.ms2rescore``.
    psm_list
        PSMs to be used for the report. If not provided, the PSMs will be read from the
        PSM file that matches the ``output_path_prefix``.
    feature_names
        Feature names to be used for the report. If not provided, the feature names will be
        read from the feature names file that matches the ``output_path_prefix``.
    use_txt_log
        If True, the log file will be read from ``output_path_prefix + ".log.txt"`` instead of
        ``output_path_prefix + ".log.html"``.
    output_file
        Path to the output HTML file. If not provided, will be ``output_path_prefix + ".report.html"``.
    use_mokapot
        If True, use mokapot LinearConfidence objects for overview charts (legacy mode).
        If False (default), use PSM dataframe directly.

    """
    files = _collect_files(output_path_prefix, use_txt_log=use_txt_log)

    # Read PSMs
    if not psm_list:
        if files["PSMs"]:
            logger.info("Reading PSMs...")
            psm_list = psm_utils.io.read_file(files["PSMs"], filetype="tsv", show_progressbar=True)
        else:
            raise FileNotFoundError("PSM file not found and no PSM list provided.")

    # Create comprehensive dataframe from PSM list
    logger.debug("Creating PSM dataframe...")
    psm_df = create_psm_dataframe(psm_list)

    # Pre-compute commonly used filtered subsets for performance
    targets = psm_df[~psm_df["is_decoy"]]
    is_decoy = psm_df["is_decoy"]  # Extract once for reuse
    if "qvalue_before" in psm_df.columns and "qvalue_after" in psm_df.columns:
        targets_before_fdr = targets[targets["qvalue_before"] <= 0.01]
        targets_after_fdr = targets[targets["qvalue_after"] <= 0.01]
    else:
        targets_before_fdr = None
        targets_after_fdr = None

    # Try to read config, but don't fail if it doesn't exist
    config = None
    if files["configuration"] and files["configuration"].is_file():
        try:
            config = json.loads(files["configuration"].read_text())
        except (json.JSONDecodeError, FileNotFoundError):
            logger.warning("Could not read configuration file. Proceeding without it.")
            config = {"ms2rescore": {}}
    else:
        logger.info("No configuration file found. Proceeding without it.")
        config = {"ms2rescore": {}}

    # Generate overview context
    if use_mokapot and config.get("ms2rescore", {}).get("fasta_file"):
        logger.debug("Recalculating confidence estimates with mokapot...")
        fasta_file = config["ms2rescore"]["fasta_file"]
        confidence_before, confidence_after = get_confidence_estimates(psm_list, fasta_file)
        overview_context = _get_overview_context(confidence_before, confidence_after)
    else:
        logger.debug("Generating overview from PSM dataframe...")
        overview_context = _get_overview_context_df(
            psm_df,
            targets=targets,
            targets_before_fdr=targets_before_fdr,
            targets_after_fdr=targets_after_fdr,
        )

    target_decoy_context = _get_target_decoy_context(psm_df)
    features_context = _get_features_context(
        psm_df, files, is_decoy=is_decoy, feature_names=feature_names
    )
    config_context = _get_config_context(config)
    log_context = _get_log_context(files)

    # Get PSM filename(s) for metadata
    if config.get("ms2rescore", {}).get("psm_file"):
        psm_filenames = "\n".join(
            [Path(id_file).name for id_file in config["ms2rescore"]["psm_file"]]
        )
    elif files["PSMs"]:
        psm_filenames = files["PSMs"].name
    else:
        psm_filenames = "Unknown"

    context = {
        "plotlyjs_version": get_plotlyjs_version(),
        "metadata": {
            "generated_on": datetime.now().strftime("%d/%m/%Y %H:%M:%S"),
            "ms2rescore_version": ms2rescore.__version__,  # TODO: Write during run?
            "psm_filename": psm_filenames,
        },
        "main_tabs": [
            {
                "id": "main_tab_comparison",
                "title": "Overview",
                "template": "overview.html",
                "context": overview_context,
            },
            {
                "id": "main_tab_target_decoy",
                "title": "Target/decoy evaluation",
                "template": "target-decoy.html",
                "context": target_decoy_context,
            },
            {
                "id": "main_tab_features",
                "title": "Rescoring features",
                "template": "features.html",
                "context": features_context,
            },
            {
                "id": "main_tab_config",
                "title": "Full configuration",
                "template": "config.html",
                "context": config_context,
            },
            {
                "id": "main_tab_log",
                "title": "Log",
                "template": "log.html",
                "context": log_context,
            },
        ],
    }

    _render_and_write(output_path_prefix, output_file=output_file, **context)


def _collect_files(output_path_prefix, use_txt_log=False):
    """Collect all files generated by MS²Rescore."""
    logger.debug("Collecting files...")
    files = {
        "PSMs": Path(output_path_prefix + ".psms.tsv").resolve(),
        "configuration": Path(output_path_prefix + ".full-config.json").resolve(),
        "feature names": Path(output_path_prefix + ".feature_names.tsv").resolve(),
        "feature weights": Path(output_path_prefix + ".mokapot.weights.tsv").resolve(),
        "log": (
            Path(output_path_prefix + ".log.txt").resolve()
            if use_txt_log
            else Path(output_path_prefix + ".log.html").resolve()
        ),
    }
    for file, path in files.items():
        if Path(path).is_file():
            logger.debug("✅ Found %s: '%s'", file, path.as_posix())
        else:
            logger.warning("❌ %s: '%s'", file, path.as_posix())
            files[file] = None
    return files


def _get_stats_context(confidence_before, confidence_after):
    """Return context for overview statistics pane."""
    stats = []
    levels = ["psms", "peptides", "proteins"]
    level_names = ["PSMs", "Peptides", "Protein groups"]
    card_colors = ["card-bg-blue", "card-bg-green", "card-bg-red"]

    # Cannot report stats if confidence estimates are not present
    if not confidence_before or not confidence_after:
        return stats

    for level, level_name, card_color in zip(levels, level_names, card_colors):
        try:
            before = confidence_before.accepted[level.lower()]
            after = confidence_after.accepted[level.lower()]
        except KeyError:
            continue  # Level not present (e.g. no fasta provided)
        if not before or not after:
            continue
        increase = (after - before) / before * 100
        stats.append(
            {
                "item": level_name,
                "card_color": card_color,
                "number": after,
                "diff": f"({after - before:+})",
                "percentage": f"{increase:.1f}%",
                "is_increase": increase > 0,
                "bar_percentage": before / after * 100 if increase > 0 else after / before * 100,
                "bar_color": "#24a143" if increase > 0 else "#a12424",
            }
        )
    return stats


def _get_stats_context_df(
    psm_df: pd.DataFrame,
    targets: Optional[pd.DataFrame] = None,
    targets_before_fdr: Optional[pd.DataFrame] = None,
    targets_after_fdr: Optional[pd.DataFrame] = None,
    fdr_threshold: float = 0.01,
) -> list:
    """Return context for overview statistics pane from dataframe."""
    stats = []

    if "qvalue_before" not in psm_df.columns or "qvalue_after" not in psm_df.columns:
        return stats

    # Use pre-computed subsets if available, otherwise compute now
    if targets is None:
        targets = psm_df[~psm_df["is_decoy"]]

    # PSM level stats - use pre-computed subsets if available
    if targets_before_fdr is not None and targets_after_fdr is not None:
        psms_before = len(targets_before_fdr)
        psms_after = len(targets_after_fdr)
    else:
        psms_before = len(targets[targets["qvalue_before"] <= fdr_threshold])
        psms_after = len(targets[targets["qvalue_after"] <= fdr_threshold])

    if psms_before > 0:
        increase = (psms_after - psms_before) / psms_before * 100
        stats.append(
            {
                "item": "PSMs",
                "card_color": "card-bg-blue",
                "number": psms_after,
                "diff": f"({psms_after - psms_before:+})",
                "percentage": f"{increase:.1f}%",
                "is_increase": increase > 0,
                "bar_percentage": psms_before / psms_after * 100
                if increase > 0
                else psms_after / psms_before * 100,
                "bar_color": "#24a143" if increase > 0 else "#a12424",
            }
        )

    # Peptide level stats
    if "peptidoform" in psm_df.columns:
        if targets_before_fdr is not None and targets_after_fdr is not None:
            peptides_before = targets_before_fdr["peptidoform"].nunique()
            peptides_after = targets_after_fdr["peptidoform"].nunique()
        else:
            peptides_before = targets[targets["qvalue_before"] <= fdr_threshold][
                "peptidoform"
            ].nunique()
            peptides_after = targets[targets["qvalue_after"] <= fdr_threshold][
                "peptidoform"
            ].nunique()

        if peptides_before > 0:
            increase = (peptides_after - peptides_before) / peptides_before * 100
            stats.append(
                {
                    "item": "Peptides",
                    "card_color": "card-bg-green",
                    "number": peptides_after,
                    "diff": f"({peptides_after - peptides_before:+})",
                    "percentage": f"{increase:.1f}%",
                    "is_increase": increase > 0,
                    "bar_percentage": peptides_before / peptides_after * 100
                    if increase > 0
                    else peptides_after / peptides_before * 100,
                    "bar_color": "#24a143" if increase > 0 else "#a12424",
                }
            )

    return stats


def _get_overview_context_df(
    psm_df: pd.DataFrame,
    targets: Optional[pd.DataFrame] = None,
    targets_before_fdr: Optional[pd.DataFrame] = None,
    targets_after_fdr: Optional[pd.DataFrame] = None,
) -> dict:
    """Return context for overview tab from dataframe."""
    logger.debug("Generating overview charts from dataframe...")
    return {
        "stats": _get_stats_context_df(
            psm_df,
            targets=targets,
            targets_before_fdr=targets_before_fdr,
            targets_after_fdr=targets_after_fdr,
        ),
        "charts": [
            {
                "title": TEXTS["charts"]["score_comparison"]["title"],
                "description": TEXTS["charts"]["score_comparison"]["description"],
                "chart": charts.score_scatter_plot_df(psm_df).to_html(**PLOTLY_HTML_KWARGS),
            },
            {
                "title": TEXTS["charts"]["fdr_comparison"]["title"],
                "description": TEXTS["charts"]["fdr_comparison"]["description"],
                "chart": charts.fdr_plot_comparison_df(psm_df).to_html(**PLOTLY_HTML_KWARGS),
            },
            {
                "title": TEXTS["charts"]["identification_overlap"]["title"],
                "description": TEXTS["charts"]["identification_overlap"]["description"],
                "chart": charts.identification_overlap_df(psm_df).to_html(**PLOTLY_HTML_KWARGS),
            },
        ],
    }


def _get_overview_context(confidence_before, confidence_after) -> dict:
    """Return context for overview tab."""
    logger.debug("Generating overview charts...")
    return {
        "stats": _get_stats_context(confidence_before, confidence_after),
        "charts": [
            {
                "title": TEXTS["charts"]["score_comparison"]["title"],
                "description": TEXTS["charts"]["score_comparison"]["description"],
                "chart": charts.score_scatter_plot(
                    confidence_before,
                    confidence_after,
                ).to_html(**PLOTLY_HTML_KWARGS),
            },
            {
                "title": TEXTS["charts"]["fdr_comparison"]["title"],
                "description": TEXTS["charts"]["fdr_comparison"]["description"],
                "chart": charts.fdr_plot_comparison(
                    confidence_before,
                    confidence_after,
                ).to_html(**PLOTLY_HTML_KWARGS),
            },
            {
                "title": TEXTS["charts"]["identification_overlap"]["title"],
                "description": TEXTS["charts"]["identification_overlap"]["description"],
                "chart": charts.identification_overlap(
                    confidence_before,
                    confidence_after,
                ).to_html(**PLOTLY_HTML_KWARGS),
            },
        ],
    }


def _get_target_decoy_context(psm_df: pd.DataFrame) -> dict:
    logger.debug("Generating target-decoy charts...")
    return {
        "charts": [
            {
                "title": TEXTS["charts"]["score_histogram"]["title"],
                "description": TEXTS["charts"]["score_histogram"]["description"],
                "chart": charts.score_histogram(psm_df).to_html(**PLOTLY_HTML_KWARGS),
            },
            {
                "title": TEXTS["charts"]["pp_plot"]["title"],
                "description": TEXTS["charts"]["pp_plot"]["description"],
                "chart": charts.pp_plot(psm_df).to_html(**PLOTLY_HTML_KWARGS),
            },
        ]
    }


def _get_features_context(
    psm_df: pd.DataFrame,
    files: Dict[str, Path],
    is_decoy: Optional[pd.Series] = None,
    feature_names: Optional[Dict[str, set]] = None,
) -> dict:
    """Return context for features tab."""
    logger.debug("Generating feature-related charts...")

    # Use pre-computed is_decoy if provided, otherwise extract
    if is_decoy is None:
        is_decoy = psm_df["is_decoy"]
    context: dict[str, list] = {"charts": []}

    # Get feature names, mapping with generator, and flat list
    if feature_names is None or not feature_names:
        # Try to read from file first
        feature_names = read_feature_names(files.get("feature names"))

        # If file doesn't exist or is empty, infer from feature generators
        if not feature_names:
            logger.info("Feature names file not found. Inferring from feature generators...")
            # Get feature columns (exclude standard PSM columns)
            standard_columns = {
                "spectrum_id",
                "run",
                "collection",
                "spectrum",
                "peptidoform",
                "precursor_mz",
                "retention_time",
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
            feature_cols = [col for col in psm_df.columns if col not in standard_columns]

            if feature_cols:
                # Build mapping of feature name -> generator by checking each generator's feature_names
                feature_to_generator = {}
                for fgen_name, fgen_class in FEATURE_GENERATORS.items():
                    try:
                        # Instantiate with empty config to get feature names
                        fgen = fgen_class()
                        fgen_features = fgen.feature_names
                        for fname in fgen_features:
                            feature_to_generator[fname] = fgen_name
                    except Exception:
                        # If instantiation fails, skip this generator
                        continue

                # Categorize features based on generator mapping
                feature_names = defaultdict(list)
                for fname in feature_cols:
                    if fname in feature_to_generator:
                        feature_names[feature_to_generator[fname]].append(fname)
                    else:
                        # Unknown feature - put in "other" category
                        feature_names["other"].append(fname)

    # If still no features, return empty context
    if not feature_names:
        logger.warning("No features found. Skipping feature charts.")
        return context

    # Convert sets to lists if needed (for compatibility with both sources)
    feature_names = {k: list(v) if isinstance(v, set) else v for k, v in feature_names.items()}

    feature_names_flat = [f_name for f_list in feature_names.values() for f_name in f_list]
    feature_names_inv = {name: gen for gen, f_list in feature_names.items() for name in f_list}

    # Fixed color map for feature generators (Okabe-Ito colorblind-safe palette)
    FEATURE_GENERATOR_COLORS = {
        "ms2pip": "#009E73",  # Bluish green
        "deeplc": "#E69F00",  # Orange
        "im2deep": "#0072B2",  # Blue
        "ms2": "#56B4E9",  # Sky blue
        "basic": "#000000",  # Black
        "psm_file": "#F0E442",  # Yellow
        "other": "#CC79A7",  # Pink
    }
    color_map = {fg: FEATURE_GENERATOR_COLORS.get(fg, "#FFFFFF") for fg in feature_names.keys()}

    # feature weights
    if not files.get("feature weights") or not files["feature weights"].is_file():
        logger.info("Feature weights file not found. Skipping feature weights plot.")
    else:
        try:
            feature_weights = pd.read_csv(files["feature weights"], sep="\t").melt(
                var_name="feature", value_name="weight"
            )
            feature_weights["feature"] = feature_weights["feature"].str.replace(
                r"^(feature:)?", "", regex=True
            )
            feature_weights["feature_generator"] = feature_weights["feature"].map(
                feature_names_inv
            )

            context["charts"].append(
                {
                    "title": TEXTS["charts"]["feature_usage"]["title"],
                    "description": TEXTS["charts"]["feature_usage"]["description"],
                    "chart": charts.feature_weights_by_generator(
                        feature_weights, color_discrete_map=color_map
                    ).to_html(**PLOTLY_HTML_KWARGS)
                    + charts.feature_weights(
                        feature_weights, color_discrete_map=color_map
                    ).to_html(**PLOTLY_HTML_KWARGS),
                }
            )
        except Exception as e:
            logger.warning(f"Could not generate feature weights plot: {e}")

    # Individual feature performance - extract features from dataframe
    features = psm_df[feature_names_flat].copy()
    _, feature_ecdf_auc = charts.calculate_feature_qvalues(features, is_decoy)
    feature_ecdf_auc["feature_generator"] = feature_ecdf_auc["feature"].map(feature_names_inv)

    context["charts"].append(
        {
            "title": TEXTS["charts"]["feature_performance"]["title"],
            "description": TEXTS["charts"]["feature_performance"]["description"],
            "chart": charts.feature_ecdf_auc_bar(
                feature_ecdf_auc, color_discrete_map=color_map
            ).to_html(**PLOTLY_HTML_KWARGS),
        }
    )

    # MS²PIP specific charts
    if "ms2pip" in feature_names and "spec_pearson_norm" in feature_names["ms2pip"]:
        context["charts"].append(
            {
                "title": TEXTS["charts"]["ms2pip_pearson"]["title"],
                "description": TEXTS["charts"]["ms2pip_pearson"]["description"],
                "chart": charts.ms2pip_correlation(features, is_decoy, psm_df["qvalue"]).to_html(
                    **PLOTLY_HTML_KWARGS
                ),
            }
        )

    # Pre-compute filtered subset for feature-specific charts (high-confidence targets)
    high_conf_mask = (~is_decoy) & (psm_df["qvalue"] <= 0.01)
    high_conf_features = features[high_conf_mask]

    # DeepLC specific charts
    if "deeplc" in feature_names:
        scatter_chart = charts.rt_scatter(
            df=high_conf_features,
            predicted_column="predicted_retention_time_best",
            observed_column="observed_retention_time_best",
        )
        baseline_chart = charts.rt_distribution_baseline(
            df=high_conf_features,
            predicted_column="predicted_retention_time_best",
            observed_column="observed_retention_time_best",
        )
        context["charts"].append(
            {
                "title": TEXTS["charts"]["deeplc_performance"]["title"],
                "description": TEXTS["charts"]["deeplc_performance"]["description"],
                "chart": scatter_chart.to_html(**PLOTLY_HTML_KWARGS)
                + baseline_chart.to_html(**PLOTLY_HTML_KWARGS),
            }
        )

    # IM2Deep specific charts
    if "im2deep" in feature_names:
        scatter_chart = charts.rt_scatter(
            df=high_conf_features,
            predicted_column="ccs_predicted_im2deep",
            observed_column="ccs_observed_im2deep",
            xaxis_label="Observed CCS",
            yaxis_label="Predicted CCS",
            plot_title="Predicted vs. observed CCS - IM2Deep",
        )

        context["charts"].append(
            {
                "title": TEXTS["charts"]["im2deep_performance"]["title"],
                "description": TEXTS["charts"]["im2deep_performance"]["description"],
                "chart": scatter_chart.to_html(**PLOTLY_HTML_KWARGS),
            }
        )

    return context


def _get_config_context(config: dict) -> dict:
    """Return context for config tab."""
    return {
        "description": TEXTS["configuration"]["description"],
        "config": json.dumps(config, indent=4),
    }


def _get_log_context(files: Dict[str, Path]) -> dict:
    """Return context for log tab."""
    if not files["log"]:
        return {"log": "<i>Log file could not be found.</i>"}

    if files["log"].suffix == ".html":
        return {"log": files["log"].read_text(encoding="utf-8")}

    if files["log"].suffix == ".txt":
        return {"log": "<pre><code>" + files["log"].read_text(encoding="utf-8") + "</code></pre>"}

    return {"log": "<i>Log file format not recognized.</i>"}


def _render_and_write(output_path_prefix: str, output_file: Optional[Path] = None, **context):
    """Render template with context and write to HTML file."""
    if output_file:
        report_path = Path(output_file).resolve()
    else:
        report_path = Path(output_path_prefix + ".report.html").resolve()

    logger.info("Writing report to %s", report_path.as_posix())

    # Use importlib.resources for PyInstaller compatibility
    template_dir = importlib.resources.files(templates)
    env = Environment(loader=FileSystemLoader(str(template_dir), encoding="utf-8"))
    template = env.get_template("base.html")
    with open(report_path, "w", encoding="utf-8") as f:
        f.write(template.render(**context))
