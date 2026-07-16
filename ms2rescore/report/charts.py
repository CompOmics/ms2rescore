"""Collection of Plotly-based charts for reporting results of MS²Rescore."""

import importlib.resources
import warnings
from collections import defaultdict
from typing import Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import plotly.subplots
import pyteomics.auxiliary
from numpy.typing import ArrayLike
from psm_utils.psm_list import PSMList
from ristretto import RescoreResult

# Fixed color per feature generator (ColorBrewer "Dark2", colorblind-safe and mutually distinct).
# Used both for the feature-generator overview charts and for the generator-specific charts, so a
# generator keeps the same color everywhere in the report.
FEATURE_GENERATOR_COLORS = {
    "ms2pip": "#1B9E77",  # Teal
    "deeplc": "#D95F02",  # Orange
    "im2deep": "#3C93C2",  # Blue
    "ms2": "#7570B3",  # Violet
    "basic": "#666666",  # Gray
    "psm_file": "#E6AB02",  # Gold
    "other": "#66A61E",  # Olive
}

# Semantic colors reused across charts.
_COLOR_TARGET = "#2c6fbb"  # Blue
_COLOR_DECOY = "#c0392b"  # Red
_COLOR_REFERENCE = "#7a7a7a"  # Neutral gray for reference/identity lines
_COLOR_NEUTRAL = "#a7b3bf"  # Muted slate for background distributions

# Categorical color sequence (Dark2) for charts without an explicit mapping.
_COLORWAY = list(FEATURE_GENERATOR_COLORS.values())

# Shared Plotly template giving every chart the same typographic and grid style as the report.
_TEMPLATE = go.layout.Template(
    layout=go.Layout(
        font=dict(family="Lato, sans-serif", size=13, color="#2b2b2b"),
        title=dict(
            font=dict(family="Oswald, sans-serif", size=18, color="#1a1a2e"),
            x=0.02,
            xanchor="left",
        ),
        paper_bgcolor="white",
        plot_bgcolor="white",
        colorway=_COLORWAY,
        margin=dict(l=60, r=30, t=60, b=50),
        xaxis=dict(
            gridcolor="#ececec",
            zeroline=False,
            showline=True,
            linecolor="#cfcfcf",
            ticks="outside",
            tickcolor="#cfcfcf",
            ticklen=4,
            automargin=True,
        ),
        yaxis=dict(
            gridcolor="#ececec",
            zeroline=False,
            showline=True,
            linecolor="#cfcfcf",
            ticks="outside",
            tickcolor="#cfcfcf",
            ticklen=4,
            automargin=True,
        ),
        legend=dict(
            bgcolor="rgba(255, 255, 255, 0.7)",
            bordercolor="#e0e0e0",
            borderwidth=1,
        ),
        hoverlabel=dict(font=dict(family="Lato, sans-serif", size=12), bordercolor="white"),
    )
)


def _style(fig: go.Figure) -> go.Figure:
    """Apply the shared MS²Rescore chart template to a figure and return it."""
    fig.update_layout(template=_TEMPLATE)
    return fig


class _ECDF:
    """
    Return the Empirical CDF of an array as a step function.

    Parameters
    ----------
    x : array_like
        Observations
    """

    def __init__(self, x):
        # Get ECDF
        x = np.array(x, copy=True)
        x.sort()
        nobs = len(x)
        y = np.linspace(1.0 / nobs, 1, nobs)

        # Make into step function
        _x = np.asarray(x)
        _y = np.asarray(y)

        if _x.shape != _y.shape:
            msg = "x and y do not have the same shape"
            raise ValueError(msg)
        if len(_x.shape) != 1:
            msg = "x and y must be 1-dimensional"
            raise ValueError(msg)

        self.x = np.r_[-np.inf, _x]
        self.y = np.r_[0.0, _y]
        self.n = self.x.shape[0]

    def __call__(self, time):
        tind = np.searchsorted(self.x, time, side="right") - 1
        return self.y[tind]


def score_histogram(psms: Union[PSMList, pd.DataFrame]) -> go.Figure:
    """
    Plot histogram of scores for a single PSM dataset.

    Parameters
    ----------
    psms
        PSMs to plot, as :py:class:`psm_utils.PSMList` or :py:class:`pandas.DataFrame` generated
        with :py:meth:`psm_utils.PSMList.to_dataframe`.

    """
    if isinstance(psms, PSMList):
        psm_df = psms.to_dataframe()
    else:
        psm_df = psms

    is_decoy = psm_df["is_decoy"].map({True: "decoy", False: "target"})

    fig = px.histogram(
        psm_df,
        x="score",
        color=is_decoy,
        barmode="overlay",
        histnorm="",
        labels={"is_decoy": "PSM type", "False": "target", "True": "decoy"},
        opacity=0.6,
        color_discrete_map={"target": _COLOR_TARGET, "decoy": _COLOR_DECOY},
    )

    # Get score thresholds
    if all(psm_df["qvalue"]):
        try:
            score_threshold = (
                psm_df[psm_df["qvalue"] <= 0.01]
                .sort_values("qvalue", ascending=False)["qvalue"]
                .iloc[0]
            )
        except IndexError:  # No PSMs below threshold
            pass
        else:
            fig.add_vline(x=score_threshold, line_dash="dash", line_color=_COLOR_REFERENCE)

    return _style(fig)


def pp_plot(psms: Union[PSMList, pd.DataFrame]) -> go.Figure:
    """
    Generate PP plot of target and decoy score distributions.

    Parameters
    ----------
    psms
        PSMs to plot, as :py:class:`psm_utils.PSMList` or :py:class:`pandas.DataFrame` generated
        with :py:meth:`psm_utils.PSMList.to_dataframe`.

    """
    if isinstance(psms, PSMList):
        psm_df = psms.to_dataframe()
    else:
        psm_df = psms

    n_decoys = np.count_nonzero(psm_df["is_decoy"])
    n_targets = len(psm_df) - n_decoys
    pi_zero = n_decoys / n_targets
    if n_decoys == 0:
        raise ValueError("No decoy PSMs found in PSM file.")
    target_scores = psm_df["score"][~psm_df["is_decoy"]]
    decoy_scores = psm_df["score"][psm_df["is_decoy"]]
    if len(psm_df) > 1000:
        target_scores_quantiles = psm_df["score"][~psm_df["is_decoy"]].quantile(
            np.linspace(0, 1, 1000)
        )
    else:
        target_scores_quantiles = target_scores
    target_ecdf = _ECDF(target_scores)(target_scores_quantiles)
    decoy_ecdf = _ECDF(decoy_scores)(target_scores_quantiles)

    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=decoy_ecdf,
            y=target_ecdf,
            mode="markers",
            marker=dict(color=_COLOR_TARGET),
        )
    )
    fig.add_trace(
        go.Scatter(
            x=[0, 1],
            y=[0, pi_zero],
            mode="lines",
            line=go.scatter.Line(color=_COLOR_REFERENCE, dash="dash"),
            showlegend=True,
            name="pi0",
        )
    )
    fig.update_layout(
        xaxis_title="Fdp",
        yaxis_title="Ftp",
        showlegend=False,
    )
    return _style(fig)


def fdr_plot(
    psms: Union[PSMList, pd.DataFrame],
    fdr_thresholds: Optional[List[float]] = None,
    log: bool = True,
) -> go.Figure:
    """
    Plot number of identifications in function of FDR threshold.

    Parameters
    ----------
    psms
        PSMs to plot, as :py:class:`psm_utils.PSMList` or :py:class:`pandas.DataFrame` generated
        with :py:meth:`psm_utils.PSMList.to_dataframe`.
    fdr_thresholds
        List of FDR thresholds to draw as vertical lines.
    log
        Whether to plot the x-axis on a log scale. Defaults to ``True``.

    """
    if isinstance(psms, PSMList):
        psm_df = psms.to_dataframe()
    else:
        psm_df = psms

    df = (
        psm_df[~psm_df["is_decoy"]]
        .reset_index(drop=True)
        .sort_values("qvalue", ascending=True)
        .copy()
    )
    df["count"] = (~df["is_decoy"]).cumsum()
    fig = px.line(
        df,
        x="qvalue",
        y="count",
        log_x=log,
        labels={"count": "Number of identified target PSMs", "qvalue": "FDR threshold"},
        color_discrete_sequence=[_COLOR_TARGET],
    )
    if fdr_thresholds:
        for threshold in fdr_thresholds:
            fig.add_vline(x=threshold, line_dash="dash", line_color=_COLOR_REFERENCE)
    return _style(fig)


def feature_weights(
    feature_weights: pd.DataFrame, color_discrete_map: Optional[Dict[str, str]] = None
) -> go.Figure:
    """
    Plot bar chart of feature weights.

    Parameters
    ----------
    feature_weights
        Data frame with columns ``feature``, ``feature_generator``, and ``weight``.
    color_discrete_map
        Mapping of feature generator names to colors for plotting.

    """
    bar_data = (
        feature_weights.groupby(["feature", "feature_generator"])
        .median(numeric_only=True)
        .abs()
        .sort_values("weight")
        .reset_index()
    )

    fig = px.bar(
        data_frame=bar_data,
        x="weight",
        y="feature",
        color="feature_generator",
        orientation="h",
        hover_name="feature",
        title="Absolute median weights by feature",
        labels={
            "weight": "Absolute median weight",
            "feature_generator": "Feature generator",
            "feature": "Feature",
        },
        color_discrete_map=color_discrete_map,
    )
    return _style(fig)


def feature_weights_by_generator(
    feature_weights: pd.DataFrame, color_discrete_map: Optional[Dict[str, str]] = None
) -> go.Figure:
    """
    Plot bar chart of feature weights, summed by feature generator.

    Parameters
    ----------
    feature_weights
        Data frame with columns "feature", "feature_generator", and "weight".
    color_discrete_map
        Mapping of feature generator names to colors for plotting.

    """
    bar_data = (
        feature_weights.groupby(["feature", "feature_generator"])
        .median()
        .abs()
        .reset_index()
        .groupby("feature_generator")
        .sum(numeric_only=True)
        .reset_index()
        .sort_values("weight")
    )

    fig = px.bar(
        data_frame=bar_data,
        x="weight",
        y="feature_generator",
        color="feature_generator",
        orientation="h",
        hover_name="feature_generator",
        title="Absolute median weights, summed by feature generator",
        labels={
            "weight": "Absolute median weight",
            "feature_generator": "Feature generator",
            "feature": "Feature",
        },
        color_discrete_map=color_discrete_map,
    )
    return _style(fig)


def ms2pip_correlation(
    features: pd.DataFrame,
    is_decoy: Union[pd.Series, np.ndarray],
    qvalue: Union[pd.Series, np.ndarray],
    color: Optional[str] = None,
) -> go.Figure:
    """
    Plot MS²PIP correlation for target PSMs with q-value <= 0.01.

    Parameters
    ----------
    features
        Data frame with features. Must contain the column ``spec_pearson_norm``.
    is_decoy
        Boolean array indicating whether each PSM is a decoy.
    qvalue
        Array of q-values for each PSM.
    color
        Bar color. Defaults to the MS²PIP feature-generator color.

    """
    data = features["spec_pearson_norm"][(qvalue < 0.01) & (~is_decoy)]
    fig = px.histogram(
        x=data,
        labels={"x": "Pearson correlation"},
        color_discrete_sequence=[color or FEATURE_GENERATOR_COLORS["ms2pip"]],
    )
    # Draw vertical line at median
    fig.add_vline(
        x=data.median(),
        line_width=3,
        line_dash="dash",
        line_color=_COLOR_REFERENCE,
        annotation_text=f"Median: {data.median():.2f}",
        annotation_position="top left",
    )
    return _style(fig)


def calculate_feature_qvalues(
    features: pd.DataFrame,
    is_decoy: ArrayLike,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Calculate q-values and ECDF AUC for all rescoring features.

    Q-values are calculated for each feature as if it was directly used PSM score. For each q-value
    distribution, the ECDF AUC is calculated as a measure of overall individual performance of the
    feature.

    As it is not known whether higher or lower values are better for each feature, q-values are
    calculated for both the original and reversed scores. The q-values and ECDF AUC are returned
    for the calculation with the highest ECDF AUC.

    Parameters
    ----------
    features
        Data frame with features. Must contain the column ``spec_pearson_norm``.
    is_decoy
        Boolean array indicating whether each PSM is a decoy.

    Returns
    -------
    feature_qvalues
        Wide-form data frame with q-values for each feature.
    feature_ecdf_auc
        Long-form data frame with ECDF AUC for each feature.

    """
    feature_qvalues = dict()
    feature_ecdf_auc = dict()
    for fname in features:
        # Calculate q-values for reversed and non-reversed scores
        q_values = []
        for reverse in [False, True]:
            with warnings.catch_warnings():  # Ignore divide by zero warning
                warnings.simplefilter("ignore", category=RuntimeWarning)
                q_values.append(
                    pyteomics.auxiliary.qvalues(
                        features,
                        key=fname,
                        is_decoy=is_decoy,
                        remove_decoy=True,
                        reverse=reverse,
                        full_output=False,
                    )["q"]
                )

        # Calculate ECDF AUC as measure of overall individual performance of feature
        ecdf_aucs = []
        for q in q_values:
            sorted_q = np.sort(q)
            y_vals = np.max(q) - sorted_q
            if hasattr(np, "trapezoid"):
                auc = np.trapezoid(y_vals)  # Numpy 2.0 and later
            else:
                auc = np.trapz(y_vals)  # type: ignore[reportAttributeAccessIssue] # Numpy 1.x
            ecdf_aucs.append(auc)

        # Select and save q-value calculation with best AUC (score reversed or not)
        idx_best = np.argmax(ecdf_aucs)
        feature_qvalues[fname] = q_values[idx_best]
        feature_ecdf_auc[fname] = ecdf_aucs[idx_best]

    # Restructure as data frames
    feature_qvalues = pd.DataFrame(feature_qvalues)
    feature_ecdf_auc = (
        pd.DataFrame([feature_ecdf_auc])
        .transpose()
        .reset_index()
        .rename(columns={"index": "feature", 0: "ecdf_auc"})
    )

    return feature_qvalues, feature_ecdf_auc


def feature_ecdf_auc_bar(
    feature_ecdf_auc: pd.DataFrame, color_discrete_map: Optional[Dict[str, str]] = None
) -> go.Figure:
    """
    Plot bar chart of feature q-value ECDF AUCs.

    Parameters
    ----------
    feature_ecdf_auc
        Data frame with columns ``feature``, ``feature_generator``, and ``ecdf_auc``.
    color_discrete_map
        Mapping of feature generator names to colors for plotting.

    """
    fig = px.bar(
        data_frame=feature_ecdf_auc.sort_values("ecdf_auc", ascending=True),
        x="ecdf_auc",
        y="feature",
        color="feature_generator",
        orientation="h",
        hover_name="feature",
        labels={
            "ecdf_auc": "Q-value ECDF AUC",
            "feature_generator": "Feature generator",
            "feature": "Feature",
        },
        color_discrete_map=color_discrete_map,
    )
    return _style(fig)


def rt_scatter(
    df: pd.DataFrame,
    predicted_column: str = "Predicted retention time",
    observed_column: str = "Observed retention time",
    xaxis_label: str = "Observed retention time",
    yaxis_label: str = "Predicted retention time",
    plot_title: str = "Predicted vs. observed retention times",
    marker_color: Optional[str] = None,
) -> go.Figure:
    """
    Plot a scatter plot of the predicted vs. observed retention times.

    Parameters
    ----------
    df : pd.DataFrame
        Dataframe containing the predicted and observed retention times.
    predicted_column : str, optional
        Name of the column containing the predicted retention times, by default
        ``Predicted retention time``.
    observed_column : str, optional
        Name of the column containing the observed retention times, by default
        ``Observed retention time``.
    xaxis_label : str, optional
        X-axis label, by default ``Observed retention time``.
    yaxis_label : str, optional
        Y-axis label, by default ``Predicted retention time``.
    plot_title : str, optional
        Scatter plot title, by default ``Predicted vs. observed retention times``
    marker_color : str, optional
        Color of the scatter points. Defaults to the Plotly template color. Pass the feature
        generator color to match the point color to the rest of the report.

    """
    # Draw scatter
    fig = px.scatter(
        df,
        x=observed_column,
        y=predicted_column,
        opacity=0.3,
        color_discrete_sequence=[marker_color] if marker_color else None,
    )

    # Draw diagonal reference line
    fig.add_scatter(
        x=[min(df[observed_column]), max(df[observed_column])],
        y=[min(df[observed_column]), max(df[observed_column])],
        mode="lines",
        line=dict(color=_COLOR_REFERENCE, width=2, dash="dash"),
    )

    # Hide legend
    fig.update_layout(
        title=plot_title,
        showlegend=False,
        xaxis_title=xaxis_label,
        yaxis_title=yaxis_label,
    )

    return _style(fig)


def rt_distribution_baseline(
    df: pd.DataFrame,
    predicted_column: str = "Predicted retention time",
    observed_column: str = "Observed retention time",
    highlight_color: Optional[str] = None,
) -> go.Figure:
    """
    Plot a distribution plot of the relative mean absolute error of the current
    DeepLC performance compared to the baseline performance.

    Parameters
    ----------
    df : pd.DataFrame
        Dataframe containing the predicted and observed retention times.
    predicted_column : str, optional
        Name of the column containing the predicted retention times, by default
        ``Predicted retention time``.
    observed_column : str, optional
        Name of the column containing the observed retention times, by default
        ``Observed retention time``.
    highlight_color : str, optional
        Color of the current-performance line. Defaults to the DeepLC feature-generator color.

    """
    # Get baseline data from deeplc package
    try:
        import deeplc.package_data

        baseline_ref = (
            importlib.resources.files(deeplc.package_data)
            / "baseline_performance"
            / "baseline_predictions.csv"
        )
        with importlib.resources.as_file(baseline_ref) as baseline_path:
            baseline_df = pd.read_csv(baseline_path)
    except (ImportError, FileNotFoundError):
        # If deeplc is not installed or baseline data not found, return empty figure
        fig = go.Figure()
        fig.add_annotation(
            text="DeepLC baseline data not available. Install DeepLC to view performance comparison.",
            showarrow=False,
        )
        return _style(fig)

    baseline_df["rel_mae_best"] = baseline_df[
        ["rel_mae_transfer_learning", "rel_mae_new_model", "rel_mae_calibrate"]
    ].min(axis=1)
    baseline_df.fillna(0.0, inplace=True)

    # Calculate current RMAE and percentile compared to baseline
    mae = sum(abs(df[observed_column] - df[predicted_column])) / len(df.index)
    mae_rel = (mae / max(df[observed_column])) * 100
    percentile = round((baseline_df["rel_mae_transfer_learning"] < mae_rel).mean() * 100, 1)

    # Calculate x-axis range with 5% padding
    all_values = np.append(baseline_df["rel_mae_transfer_learning"].values, mae_rel)
    padding = (all_values.max() - all_values.min()) / 20  # 5% padding
    x_min = all_values.min() - padding
    x_max = all_values.max() + padding

    # Make labels human-readable
    hover_label_mapping = {
        "train_number": "Training dataset size",
        "rel_mae_transfer_learning": "RMAE with transfer learning",
        "rel_mae_new_model": "RMAE with new model from scratch",
        "rel_mae_calibrate": "RMAE with calibrating existing model",
        "rel_mae_best": "RMAE with best method",
    }
    label_mapping = hover_label_mapping.copy()
    label_mapping.update({"Unnamed: 0": "Dataset"})

    # Generate plot
    fig = px.histogram(
        data_frame=baseline_df,
        x="rel_mae_best",
        marginal="rug",
        hover_data=hover_label_mapping.keys(),
        hover_name="Unnamed: 0",
        labels=label_mapping,
        opacity=0.8,
        color_discrete_sequence=[_COLOR_NEUTRAL],
    )
    fig.add_vline(
        x=mae_rel,
        line_width=3,
        line_dash="dash",
        line_color=highlight_color or FEATURE_GENERATOR_COLORS["deeplc"],
        annotation_text=f"Current performance (percentile {percentile}%)",
        annotation_position="top left",
        name="Current performance",
        row=1,
    )
    fig.update_xaxes(range=[x_min, x_max])
    fig.update_layout(
        title=(f"Current DeepLC performance compared to {len(baseline_df.index)} datasets"),
        xaxis_title="Relative mean absolute error (%)",
    )

    return _style(fig)


def _best_per_spectrum(psms: pd.DataFrame) -> pd.DataFrame:
    """
    Collapse a ``RescoreResult.psms`` table to the single best-scoring row per spectrum.

    Groups by ``(run, spectrum_id)``, not just ``spectrum_id`` (which can collide across
    input files). Needed wherever a chart wants one point/row per physical spectrum: under
    ``max_psm_rank_output > 1``, a spectrum can have multiple candidate rows.

    """
    return psms.sort_values("score", ascending=False).drop_duplicates(
        subset=["run", "spectrum_id"], keep="first"
    )


def score_scatter_plot(
    before: RescoreResult,
    after: RescoreResult,
    fdr_threshold: float = 0.01,
) -> go.Figure:
    """
    Plot PSM scores before and after rescoring, best-scoring candidate per spectrum.

    Collapses ``before``/``after`` to one row per ``(run, spectrum_id)`` -- the
    best-scoring candidate on each side, independently -- before comparing. Under
    ``max_psm_rank_output > 1`` either side can have multiple candidate rows per spectrum,
    and rescoring can legitimately promote a different peptidoform as a spectrum's winner,
    so this is a spectrum-level comparison, not a peptidoform-level one.

    Parameters
    ----------
    before
        Result of evaluating the PSMs' pre-rescoring score with ristretto.
    after
        Result of rescoring the PSMs with ristretto.
    fdr_threshold
        FDR threshold for drawing threshold lines.

    Returns
    -------
    go.Figure
        Plotly figure with score comparison.

    """
    if before.psms.empty or after.psms.empty:
        figure = go.Figure()
        figure.add_annotation(
            text="No before/after score data available for comparison.",
            showarrow=False,
        )
        return _style(figure)

    before_best = _best_per_spectrum(before.psms)[["run", "spectrum_id", "score", "qvalue"]].rename(
        columns={"score": "score_before", "qvalue": "qvalue_before"}
    )
    after_best = _best_per_spectrum(after.psms)[
        ["run", "spectrum_id", "score", "qvalue", "is_decoy"]
    ].rename(columns={"score": "score_after", "qvalue": "qvalue_after"})
    plot_df = before_best.merge(after_best, on=["run", "spectrum_id"], how="inner")
    plot_df["PSM type"] = plot_df["is_decoy"].map({True: "decoy", False: "target"})

    # Get score thresholds
    try:
        score_threshold_before = (
            plot_df[plot_df["qvalue_before"] <= fdr_threshold]
            .sort_values("qvalue_before", ascending=False)["score_before"]
            .iloc[0]
        )
    except IndexError:
        score_threshold_before = None

    try:
        score_threshold_after = (
            plot_df[plot_df["qvalue_after"] <= fdr_threshold]
            .sort_values("qvalue_after", ascending=False)["score_after"]
            .iloc[0]
        )
    except IndexError:
        score_threshold_after = None

    # Plot
    fig = px.scatter(
        data_frame=plot_df,
        x="score_before",
        y="score_after",
        color="PSM type",
        marginal_x="histogram",
        marginal_y="histogram",
        opacity=0.1,
        labels={
            "score_before": "Spectrum score (before rescoring)",
            "score_after": "Spectrum score (after rescoring)",
        },
        color_discrete_map={"target": _COLOR_TARGET, "decoy": _COLOR_DECOY},
    )

    # Draw FDR thresholds
    if score_threshold_before:
        fig.add_vline(x=score_threshold_before, line_dash="dash", line_color=_COLOR_REFERENCE, row=1, col=1)
        fig.add_vline(x=score_threshold_before, line_dash="dash", line_color=_COLOR_REFERENCE, row=2, col=1)
    if score_threshold_after:
        fig.add_hline(y=score_threshold_after, line_dash="dash", line_color=_COLOR_REFERENCE, row=1, col=1)
        fig.add_hline(y=score_threshold_after, line_dash="dash", line_color=_COLOR_REFERENCE, row=1, col=2)

    return _style(fig)


def fdr_plot_comparison(
    before: RescoreResult, after: RescoreResult, fdr_threshold: float = 0.01
) -> go.Figure:
    """
    Plot number of identified spectra as a function of FDR threshold, before vs. after.

    Collapses each of ``before``/``after`` independently to one row per ``(run,
    spectrum_id)`` -- the best-scoring candidate -- before counting, so a spectrum with
    multiple ambiguous candidate rows (``max_psm_rank_output > 1``) isn't counted more than
    once.

    Parameters
    ----------
    before
        Result of evaluating the PSMs' pre-rescoring score with ristretto.
    after
        Result of rescoring the PSMs with ristretto.
    fdr_threshold
        FDR threshold to draw as a reference line.

    Returns
    -------
    go.Figure
        Plotly figure with FDR comparison.

    """
    if before.psms.empty or after.psms.empty:
        figure = go.Figure()
        figure.add_annotation(
            text="No before/after q-value data available for comparison.",
            showarrow=False,
        )
        return _style(figure)

    before_best = _best_per_spectrum(before.psms)
    after_best = _best_per_spectrum(after.psms)

    plot_data = pd.concat(
        [
            before_best.loc[~before_best["is_decoy"], ["qvalue"]]
            .rename(columns={"qvalue": "q-value"})
            .assign(**{"before/after": "before rescoring"}),
            after_best.loc[~after_best["is_decoy"], ["qvalue"]]
            .rename(columns={"qvalue": "q-value"})
            .assign(**{"before/after": "after rescoring"}),
        ]
    )

    # Plot
    fig = px.ecdf(
        data_frame=plot_data,
        x="q-value",
        color="before/after",
        log_x=True,
        ecdfnorm=None,
        labels={
            "q-value": "FDR threshold",
            "before/after": "",
        },
        color_discrete_map={
            "before rescoring": _COLOR_NEUTRAL,
            "after rescoring": "#24a143",
        },
    )
    fig.add_vline(x=fdr_threshold, line_dash="dash", line_color=_COLOR_REFERENCE)
    fig.update_layout(yaxis_title="Identified spectra")
    return _style(fig)


def _group_keys(df: pd.DataFrame, group_cols: Union[str, List[str]]) -> list:
    """Build hashable group keys from one column, or a compound key from several."""
    if isinstance(group_cols, str):
        return list(df[group_cols])
    return list(zip(*(df[c] for c in group_cols)))


def identification_overlap(
    before: RescoreResult,
    after: RescoreResult,
    fdr_threshold: float = 0.01,
) -> go.Figure:
    """
    Plot stacked bar charts of removed, retained, and gained IDs at each rollup level.

    Compares ristretto's own before/after rollup tables directly -- spectrum, peptidoform,
    peptide, and (optionally) protein -- rather than re-deriving sets from a merged per-PSM
    dataframe. The latter would only be correct at the peptidoform/peptide/protein level if
    every spectrum kept the same winning peptidoform between before and after, which is
    exactly what rescoring is expected to change for at least some spectra.

    Parameters
    ----------
    before
        Result of evaluating the PSMs' pre-rescoring score with ristretto.
    after
        Result of rescoring the PSMs with ristretto.
    fdr_threshold
        FDR threshold for counting identifications.

    Returns
    -------
    go.Figure
        Plotly figure with identification overlap.

    """
    levels = [
        ("spectra", before.psms, after.psms, ["run", "spectrum_id"]),
        ("peptidoforms", before.peptidoforms, after.peptidoforms, "peptidoform"),
    ]
    if before.peptides is not None and after.peptides is not None:
        levels.append(("peptides", before.peptides, after.peptides, "peptide"))
    if before.proteins is not None and after.proteins is not None:
        levels.append(("protein groups", before.proteins, after.proteins, "protein"))

    overlap_data = defaultdict(dict)
    for level_name, before_df, after_df, group_cols in levels:
        if before_df.empty or after_df.empty:
            continue
        before_mask = (before_df["qvalue"] <= fdr_threshold) & ~before_df["is_decoy"]
        after_mask = (after_df["qvalue"] <= fdr_threshold) & ~after_df["is_decoy"]

        ids_before = set(_group_keys(before_df[before_mask], group_cols))
        ids_after = set(_group_keys(after_df[after_mask], group_cols))
        overlap_data["removed"][level_name] = -len(ids_before - ids_after)
        overlap_data["retained"][level_name] = len(ids_after & ids_before)
        overlap_data["gained"][level_name] = len(ids_after - ids_before)

    if not overlap_data["retained"]:
        figure = go.Figure()
        figure.add_annotation(
            text="No before/after data available for comparison.",
            showarrow=False,
        )
        return _style(figure)

    colors = [_COLOR_DECOY, _COLOR_TARGET, "#24a143"]
    level_names = list(overlap_data["retained"].keys())
    fig = plotly.subplots.make_subplots(rows=len(level_names), cols=1)

    for i, level_name in enumerate(level_names):
        for (item, data), color in zip(overlap_data.items(), colors):
            if level_name not in data:
                continue
            fig.add_trace(
                go.Bar(
                    y=[level_name],
                    x=[data[level_name]],
                    marker={"color": color},
                    orientation="h",
                    width=0.4,
                    name=item,
                    showlegend=True if i == 0 else False,
                ),
                row=i + 1,
                col=1,
            )
    fig.update_layout(barmode="relative")

    return _style(fig)
