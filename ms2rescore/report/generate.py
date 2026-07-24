"""Generate an HTML report with various QC charts for MS²Rescore results."""

import importlib.resources
import json
import logging
import tomllib
from datetime import datetime
from pathlib import Path

import pandas as pd
from jinja2 import Environment, FileSystemLoader
from plotly.offline import get_plotlyjs_version
from ristretto import RescoreResult

import ms2rescore
from ms2rescore.report import charts, templates
from ms2rescore.report.data import ReportData

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

# Fixed color per feature generator, defined alongside the charts so the generator-specific
# charts (DeepLC, IM2Deep, MS²PIP) reuse the same color as the feature-generator overview charts.
FEATURE_GENERATOR_COLORS = charts.FEATURE_GENERATOR_COLORS

TEXTS = tomllib.loads(importlib.resources.files(templates).joinpath("texts.toml").read_text())


def generate_report(
    output_path_prefix: str,
    data: ReportData,
    output_file: Path | None = None,
):
    """
    Generate the HTML report from an in-memory :py:class:`~ms2rescore.report.data.ReportData`.

    Parameters
    ----------
    output_path_prefix
        Prefix of the MS²Rescore output file names, used to locate the log file and to derive the
        default report path. For example, if the output PSM file is
        ``/path/to/file.ms2rescore.tsv``, the prefix is ``/path/to/file.ms2rescore``.
    data
        Fully-populated report data. Build it with :py:meth:`ReportData.from_run` (in-memory run)
        or :py:meth:`ReportData.from_files` (standalone from output files).
    output_file
        Path to the output HTML file. Defaults to ``output_path_prefix + ".report.html"``.
    """
    psm_df = data.psm_df
    is_decoy = psm_df["is_decoy"]

    context = {
        "plotlyjs_version": get_plotlyjs_version(),
        "metadata": {
            "generated_on": datetime.now().strftime("%d/%m/%Y %H:%M:%S"),
            "ms2rescore_version": ms2rescore.__version__,
            "psm_filename": _get_psm_filenames(data),
        },
        "main_tabs": [
            {
                "id": "main_tab_comparison",
                "title": "Overview",
                "template": "overview.html",
                "context": _get_overview_context(
                    data.id_stats, data.before, data.after, data.fdr_threshold
                ),
            },
            {
                "id": "main_tab_target_decoy",
                "title": "Target/decoy evaluation",
                "template": "target-decoy.html",
                "context": _get_target_decoy_context(psm_df, data.fdr_threshold),
            },
            {
                "id": "main_tab_features",
                "title": "Rescoring features",
                "template": "features.html",
                "context": _get_features_context(
                    psm_df, data.feature_names, data.feature_weights, is_decoy, data.fdr_threshold
                ),
            },
            {
                "id": "main_tab_config",
                "title": "Full configuration",
                "template": "config.html",
                "context": _get_config_context(data.config),
            },
            {
                "id": "main_tab_log",
                "title": "Log",
                "template": "log.html",
                "context": _get_log_context(output_path_prefix, data.log_html),
            },
        ],
    }

    _render_and_write(output_path_prefix, output_file=output_file, **context)


def _get_psm_filenames(data: ReportData) -> str:
    """Return the input PSM filename(s) for the report metadata."""
    psm_files = data.config.get("ms2rescore", {}).get("psm_file")
    if psm_files:
        return "\n".join(Path(psm_file).name for psm_file in psm_files)
    return "Unknown"


def _get_overview_context(
    id_stats: list, before: RescoreResult, after: RescoreResult, fdr_threshold: float
) -> dict:
    """Return context for the overview tab."""
    logger.debug("Generating overview charts...")
    return {
        "stats": id_stats,
        "charts": [
            {
                "title": TEXTS["charts"]["score_comparison"]["title"],
                "description": TEXTS["charts"]["score_comparison"]["description"].format(
                    fdr_threshold=fdr_threshold
                ),
                "chart": charts.score_scatter_plot(before, after, fdr_threshold).to_html(
                    **PLOTLY_HTML_KWARGS
                ),
            },
            {
                "title": TEXTS["charts"]["fdr_comparison"]["title"],
                "description": TEXTS["charts"]["fdr_comparison"]["description"].format(
                    fdr_threshold=fdr_threshold
                ),
                "chart": charts.fdr_plot_comparison(before, after, fdr_threshold).to_html(
                    **PLOTLY_HTML_KWARGS
                ),
            },
            {
                "title": TEXTS["charts"]["identification_overlap"]["title"],
                "description": TEXTS["charts"]["identification_overlap"]["description"],
                "chart": charts.identification_overlap(before, after, fdr_threshold).to_html(
                    **PLOTLY_HTML_KWARGS
                ),
            },
        ],
    }


def _get_target_decoy_context(psm_df: pd.DataFrame, fdr_threshold: float) -> dict:
    """Return context for the target/decoy tab."""
    logger.debug("Generating target-decoy charts...")
    return {
        "charts": [
            {
                "title": TEXTS["charts"]["score_histogram"]["title"],
                "description": TEXTS["charts"]["score_histogram"]["description"].format(
                    fdr_threshold=fdr_threshold
                ),
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
    feature_names: dict,
    feature_weights: pd.DataFrame | None,
    is_decoy: pd.Series,
    fdr_threshold: float,
) -> dict:
    """Return context for the rescoring-features tab."""
    logger.debug("Generating feature-related charts...")
    context: dict = {"charts": []}

    if not feature_names:
        logger.warning("No features found. Skipping feature charts.")
        return context

    feature_names_flat = [name for names in feature_names.values() for name in names]
    feature_names_inv = {name: gen for gen, names in feature_names.items() for name in names}
    color_map = {gen: FEATURE_GENERATOR_COLORS.get(gen, "#FFFFFF") for gen in feature_names}

    # Feature weights (empty only when rescoring was skipped)
    if feature_weights is not None and not feature_weights.empty:
        _add_feature_weights_chart(context, feature_weights, feature_names_inv, color_map)

    # Individual feature performance
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

    # MS²PIP correlation
    if "spec_pearson_norm" in feature_names.get("ms2pip", []):
        context["charts"].append(
            {
                "title": TEXTS["charts"]["ms2pip_pearson"]["title"],
                "description": TEXTS["charts"]["ms2pip_pearson"]["description"].format(
                    fdr_threshold=fdr_threshold
                ),
                "chart": charts.ms2pip_correlation(
                    features, is_decoy, psm_df["qvalue"], color=color_map.get("ms2pip")
                ).to_html(**PLOTLY_HTML_KWARGS),
            }
        )

    # Retention-time and ion-mobility charts on high-confidence targets
    high_conf_features = features[(~is_decoy) & (psm_df["qvalue"] <= fdr_threshold)]
    if "deeplc" in feature_names:
        try:
            _add_deeplc_chart(
                context, high_conf_features, fdr_threshold, color=color_map.get("deeplc")
            )
        except Exception as e:
            logger.warning("Could not generate DeepLC performance plot: %s", e)
    if "im2deep" in feature_names:
        try:
            _add_im2deep_chart(context, high_conf_features, color=color_map.get("im2deep"))
        except Exception as e:
            logger.warning("Could not generate IM2Deep performance plot: %s", e)

    return context


def _add_feature_weights_chart(context, feature_weights, feature_names_inv, color_map):
    """Append the feature-weights charts to the features context."""
    try:
        # `feature_weights` is indexed by feature name, one column per CV fold
        weights = feature_weights.reset_index(names="feature").melt(
            id_vars="feature", var_name="fold", value_name="weight"
        )
        weights["feature_generator"] = weights["feature"].map(feature_names_inv)
        context["charts"].append(
            {
                "title": TEXTS["charts"]["feature_usage"]["title"],
                "description": TEXTS["charts"]["feature_usage"]["description"],
                "chart": charts.feature_weights_by_generator(
                    weights, color_discrete_map=color_map
                ).to_html(**PLOTLY_HTML_KWARGS)
                + charts.feature_weights(weights, color_discrete_map=color_map).to_html(
                    **PLOTLY_HTML_KWARGS
                ),
            }
        )
    except Exception as e:
        logger.warning("Could not generate feature weights plot: %s", e)


def _add_deeplc_chart(context, high_conf_features, fdr_threshold, color=None):
    """Append the DeepLC retention-time charts to the features context."""
    scatter = charts.rt_scatter(
        df=high_conf_features,
        predicted_column="predicted_retention_time_best",
        observed_column="observed_retention_time_best",
        marker_color=color,
    )
    baseline = charts.rt_distribution_baseline(
        df=high_conf_features,
        predicted_column="predicted_retention_time_best",
        observed_column="observed_retention_time_best",
        highlight_color=color,
    )
    context["charts"].append(
        {
            "title": TEXTS["charts"]["deeplc_performance"]["title"],
            "description": TEXTS["charts"]["deeplc_performance"]["description"].format(
                fdr_threshold=fdr_threshold
            ),
            "chart": scatter.to_html(**PLOTLY_HTML_KWARGS)
            + baseline.to_html(**PLOTLY_HTML_KWARGS),
        }
    )


def _add_im2deep_chart(context, high_conf_features, color=None):
    """Append the IM2Deep CCS chart to the features context."""
    scatter = charts.rt_scatter(
        df=high_conf_features,
        predicted_column="ccs_predicted_im2deep",
        observed_column="ccs_observed_im2deep",
        xaxis_label="Observed CCS",
        yaxis_label="Predicted CCS",
        plot_title="Predicted vs. observed CCS - IM2Deep",
        marker_color=color,
    )
    context["charts"].append(
        {
            "title": TEXTS["charts"]["im2deep_performance"]["title"],
            "description": TEXTS["charts"]["im2deep_performance"]["description"],
            "chart": scatter.to_html(**PLOTLY_HTML_KWARGS),
        }
    )


def _get_config_context(config: dict) -> dict:
    """Return context for the config tab."""
    return {
        "description": TEXTS["configuration"]["description"],
        "config": json.dumps(config, indent=4),
    }


def _get_log_context(output_path_prefix: str, log_html: str | None) -> dict:
    """Return context for the log tab, reading the log file when not provided in memory."""
    if log_html is not None:
        return {"log": log_html}

    # Locate the log written during the run, preferring the HTML variant
    for suffix in (".log.html", ".log.txt"):
        log_path = Path(output_path_prefix + suffix)
        if log_path.is_file():
            content = log_path.read_text(encoding="utf-8")
            if suffix == ".log.txt":
                content = "<pre><code>" + content + "</code></pre>"
            return {"log": content}

    return {"log": "<i>Log file could not be found.</i>"}


def _render_and_write(output_path_prefix: str, output_file: Path | None = None, **context):
    """Render the base template with context and write it to the HTML report file."""
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
