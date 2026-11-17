"""
Rescore PSMs with ristretto.

:py:mod:`ristretto` is a lean, dependency-light (numpy/scikit-learn/pandas only)
reimplementation of the Percolator/Käll 2007 semi-supervised rescoring algorithm, purpose-built
to serve MS²Rescore. It takes a plain :py:class:`pandas.DataFrame` in and returns a
:py:class:`~ristretto.result.RescoreResult` of plain DataFrames out.

"""

import logging
import re
from concurrent.futures import BrokenExecutor
from dataclasses import replace

import numpy as np
import pandas as pd
import ristretto
from psm_utils import PSMList
from ristretto import RescoreResult
from ristretto._cv import spectrum_grouped_kfold
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

from ms2rescore._ristretto_utils import (
    _build_features_dataframe,
    _is_original_psm,
    _trim_and_evaluate,
)
from ms2rescore.exceptions import RescoringError
from ms2rescore.localisation import _pep, _residue_masses, _same_fragments, localise
from ms2rescore.parse_psms import infer_score_direction

logger = logging.getLogger(__name__)


def rescore(
    psm_list: PSMList,
    config: dict,
    output_file_root: str,
    decoy_sites: PSMList | None = None,
    peaks: dict | None = None,
) -> tuple[PSMList, RescoreResult]:
    """
    Rescore PSMs with ristretto and write the new scores, q-values, and PEPs back to ``psm_list``.

    Always scores every candidate rank with ristretto's semi-supervised model
    (``multi_rank_rescoring=True``), then runs the actual competition/FDR/PEP step
    separately via :py:func:`~ms2rescore._ristretto_utils._trim_and_evaluate`, trimming to
    ``max_psm_rank_output`` PSMs per spectrum.

    If ``decoy_sites`` is given, :py:func:`rank_sites` runs on every candidate of each spectrum
    before that trim, so the rank cut keeps the candidates the site ranker prefers.

    Returns
    -------
    tuple[PSMList, RescoreResult]
        The (possibly filtered) PSMList with updated scores, and the competition/FDR result.

    """
    feature_names = {f for psm in psm_list for f in psm.rescoring_features}
    lower_score_is_better = infer_score_direction(
        psm_list, config["rescoring"].get("train_fdr", 0.01)
    )
    features_df = _build_features_dataframe(psm_list, feature_names, lower_score_is_better)

    peptide_col = "peptide"
    protein_col = "protein" if "protein" in features_df.columns else None
    decoy_pattern = config["id_decoy_pattern"]

    try:
        ml_result = ristretto.rescore(
            features_df,
            run_col="run",
            peptide_col=peptide_col,
            protein_col=protein_col,
            feature_cols=sorted(feature_names),
            decoy_pattern=decoy_pattern,
            n_jobs=int(config["processes"]),
            multi_rank_rescoring=True,
            # train_fdr/model; ristretto's own kwarg defaults apply to any key a partial
            # user-provided rescoring dict omitted.
            **config["rescoring"],
        )
        scored = ml_result.psms
        if decoy_sites is not None:
            psm_list["score"] = scored["score"].reindex(features_df.index).to_numpy()
            logger.info("Ranking modification-site candidates within spectra...")
            psm_list = rank_sites(psm_list, decoy_sites, config, output_file_root, peaks)
            scored = scored.assign(
                score=pd.Series(np.asarray(psm_list["score"], dtype=float), index=features_df.index)
            )
        final_result = _trim_and_evaluate(
            scored,
            config["max_psm_rank_output"],
            run_col="run",
            peptide_col=peptide_col,
            protein_col=protein_col,
            decoy_pattern=decoy_pattern,
        )
        # _trim_and_evaluate() calls ristretto.evaluate(), which always returns empty
        # feature_weights/n_iterations (no training occurs there) -- carry over the real
        # values from ml_result, the actual ristretto.rescore() call above.
        final_result = replace(
            final_result,
            feature_weights=ml_result.feature_weights,
            n_iterations=ml_result.n_iterations,
        )
    except (RuntimeError, BrokenExecutor, ValueError) as e:
        raise RescoringError("Ristretto could not be run. Please check the input data.") from e

    keep_mask = np.zeros(len(psm_list), dtype=bool)
    keep_mask[final_result.psms.index.to_numpy()] = True
    psm_list = psm_list[keep_mask]

    psm_list["score"] = final_result.psms["score"].to_numpy()
    psm_list["qvalue"] = final_result.psms["qvalue"].to_numpy()
    psm_list["pep"] = final_result.psms["pep"].to_numpy()
    psm_list.set_ranks(lower_score_better=False)

    psm_list, final_result = _fix_constant_pep(
        psm_list, final_result, peptide_col, protein_col, decoy_pattern
    )

    _write_group_metadata(psm_list, final_result.psms, final_result.peptidoforms, "peptidoform")
    if final_result.peptides is not None:
        _write_group_metadata(psm_list, final_result.psms, final_result.peptides, "peptide")
    if final_result.proteins is not None:
        if decoy_pattern:
            logger.debug("Protein-level FDR used picked-protein competition (decoy_pattern set).")
        else:
            logger.warning(
                "id_decoy_pattern is not set: protein-level FDR falls back to a plain rollup "
                "instead of picked-protein competition. Protein q-values/PEPs may be less "
                "accurate. Set id_decoy_pattern to enable picked-protein competition."
            )
        _write_group_metadata(
            psm_list, final_result.psms, final_result.proteins, "protein", decoy_pattern
        )

    return psm_list, final_result


def _write_group_metadata(
    psm_list: PSMList,
    psms: pd.DataFrame,
    rollup: pd.DataFrame,
    group_col: str,
    decoy_pattern: str | None = None,
) -> None:
    """
    Write a rollup's score/qvalue/pep onto each PSM's metadata, keyed by ``group_col``.

    For protein-level rollups computed with ``decoy_pattern`` set, ristretto's
    ``picked_rollup`` strips the decoy pattern from each group key (picked-protein
    competition, Savitski et al. 2015), so it must be stripped from ``psms[group_col]``
    the same way before the lookup.

    """
    lookup = rollup.set_index(group_col)[["score", "qvalue", "pep"]].to_dict(orient="index")
    strip_pattern = re.compile(decoy_pattern) if group_col == "protein" and decoy_pattern else None
    for psm, group_key in zip(psm_list, psms[group_col]):
        if strip_pattern is not None:
            group_key = strip_pattern.sub("", group_key)
        values = lookup[group_key]
        psm.metadata.update(
            {
                f"{group_col}_score": values["score"],
                f"{group_col}_qvalue": values["qvalue"],
                f"{group_col}_pep": values["pep"],
            }
        )


def _fix_constant_pep_result(
    result: RescoreResult,
    peptide_col: str | None,
    protein_col: str | None,
    decoy_pattern: str | None,
) -> tuple[RescoreResult, np.ndarray | None]:
    """
    Detect and fix constant PEP (all 1.0) on a single ``RescoreResult``.

    Removes decoy PSMs that score higher than the best target, then recomputes q-values,
    PEPs, and the peptidoform/peptide/protein rollups from scratch on the corrected
    population -- the original per-row values (and rollups derived from them) are products
    of the broken calibration and must not be carried forward.

    Returns
    -------
    tuple[RescoreResult, numpy.ndarray | None]
        The (possibly fixed) result, and the boolean keep-mask applied to ``result.psms`` --
        or ``None`` if no fix was needed/possible, in which case ``result`` is unchanged.

    """
    psms = result.psms
    if not (psms["pep"] == 1.0).all():
        return result, None

    logger.warning(
        "Attempting to fix constant PEP values by removing decoy PSMs that score higher than the "
        "best target PSM."
    )
    max_target_score = psms["score"][~psms["is_decoy"]].max()
    higher_scoring_decoys = psms["is_decoy"].to_numpy() & (
        psms["score"].to_numpy() > max_target_score
    )

    if not higher_scoring_decoys.any():
        logger.warning("No decoys scoring higher than the best target found. Skipping fix.")
        return result, None

    keep = ~higher_scoring_decoys
    logger.warning(f"Removed {higher_scoring_decoys.sum()} decoy PSMs.")
    filtered_psms = psms[keep]
    if filtered_psms["is_decoy"].nunique() < 2:
        # No decoys (or no targets) remain to recompute a competition against -- an extreme
        # edge case of an already-degenerate input. Keep the filtered rows as-is rather than
        # erroring on ristretto.evaluate()'s target/decoy requirement.
        logger.warning(
            "Cannot recompute q-values/PEP/rollups after removing decoys: no decoys remain "
            "in the corrected population. Keeping filtered rows with their prior values."
        )
        return replace(result, psms=filtered_psms), keep

    fixed = ristretto.evaluate(
        filtered_psms,
        run_col="run" if "run" in filtered_psms.columns else None,
        peptide_col=peptide_col,
        protein_col=protein_col,
        decoy_pattern=decoy_pattern,
        multi_rank_rescoring=True,
    )
    # ristretto.evaluate() always returns empty feature_weights/n_iterations (no training
    # occurs there) -- carry over the real values already present on the input result.
    fixed = replace(
        fixed, feature_weights=result.feature_weights, n_iterations=result.n_iterations
    )
    return fixed, keep


def _fix_constant_pep(
    psm_list: PSMList,
    result: RescoreResult,
    peptide_col: str | None = None,
    protein_col: str | None = None,
    decoy_pattern: str | None = None,
) -> tuple[PSMList, RescoreResult]:
    """
    Workaround for broken PEP calculation if the best-scoring PSM is a decoy.

    Filters ``psm_list`` and ``result`` together (same length, same row order) so they never
    desync -- callers (result tables, report, identification counts) all keep seeing the same
    PSM population as the returned ``psm_list``, with the same freshly recomputed rollups.

    """
    fixed_result, keep = _fix_constant_pep_result(result, peptide_col, protein_col, decoy_pattern)
    if keep is None:
        return psm_list, fixed_result

    psm_list = psm_list[keep]
    psm_list["score"] = fixed_result.psms["score"].to_numpy()
    psm_list["qvalue"] = fixed_result.psms["qvalue"].to_numpy()
    psm_list["pep"] = fixed_result.psms["pep"].to_numpy()
    return psm_list, fixed_result


# features not provided to the ranker, so they are not used for ranking
RANKER_EXCLUDED_FEATURES = {
    "abs_ms1_error_ppm",
    "mass_error",
    "theoretical_mass",
    "experimental_mass",
    "search_engine_score",
}


LOCALISATION_KEYS = [
    "localisation_set",
    "localisation_alternatives",
    "localisation_pep",
    "modification_pep",
    "localisation_confidence",
    "localisation_qvalue",
]


def _rank_explanations(rows: pd.DataFrame, feature_cols: list[str], n_folds: int = 3,
                       max_rounds: int = 5, tol: float = 0.01, C: float = 0.1, seed: int = 42):
    """
    Rank the explanations of each spectrum (semi-supervised, group-wise cross-validation).

    The top row per group is the positive and every other row of the group a negative; the start
    is the rescoring score ``score``. Features are centred within group. Only target groups with
    more than one explanation are used for training. Returns the scores and the feature weights
    per fold.
    """
    groups = pd.factorize(rows["group"])[0]
    X = np.nan_to_num(rows[feature_cols].to_numpy(dtype=np.float64))
    X -= pd.DataFrame(X).groupby(groups).transform("mean").to_numpy()
    trainable = rows["is_target"].to_numpy(bool) & (np.bincount(groups)[groups] > 1)
    if not trainable.any():  # one explanation per spectrum: nothing to rank
        return np.zeros(len(rows)), pd.DataFrame({"fold_1": 0.0}, index=feature_cols)
    top = lambda s: pd.Series(s).groupby(groups).idxmax().to_numpy()  # noqa: E731
    positives = top(np.nan_to_num(rows["score"].to_numpy(dtype=np.float64), nan=-np.inf))
    scores, weights = np.zeros(len(rows)), []
    for _ in range(max_rounds):
        y = np.zeros(len(rows), int)
        y[positives] = 1
        weights = []
        for train_idx, test_idx in spectrum_grouped_kfold(groups, n_folds, seed):
            train_idx = train_idx[trainable[train_idx]]
            scaler = StandardScaler().fit(X[train_idx])
            est = LogisticRegression(C=C, max_iter=2000).fit(scaler.transform(X[train_idx]), y[train_idx])
            scores[test_idx] = est.decision_function(scaler.transform(X[test_idx]))
            weights.append(est.coef_.ravel())
        new = top(scores)
        changed, positives = np.mean(new != positives), new
        if changed < tol:
            break
    w = pd.DataFrame({f"fold_{i}": v for i, v in enumerate(weights, start=1)}, index=feature_cols)
    return scores, w


def _prefer_common_isobars(d: pd.DataFrame, rep: pd.DataFrame, tol: float) -> None:
    """
    Among explanations with the same fragment masses as the step-1 winner, move the most frequent
    one to the top (in place, ``explanation_score``).

    Frequency = number of target spectra where the explanation wins step 1 with a hyperscore gain
    over the original hit that is unlikely under the decoy-peptide null (PEP < 0.05). Explanations
    whose fragments differ from the winner's are not affected.
    """
    if "hyperscore" not in d:
        return
    win = rep.assign(s=[d.loc[i, "explanation_score"] for i in rep.index])
    win = win.sort_values("s", ascending=False).drop_duplicates("group")
    original_hs = d[d["original"]].groupby("group")["hyperscore"].first()
    win = win.assign(gain=win["hyperscore"] - original_hs.reindex(win["group"]).to_numpy())
    mod = win[~win["original"]].dropna(subset=["gain"])
    g_t, g_0 = mod.loc[mod["is_target"], "gain"], mod.loc[~mod["is_target"], "gain"]
    if not len(g_t) or not len(g_0):
        return
    pi0 = min(1.0, (g_t <= 0).mean() / max((g_0 <= 0).mean(), 1e-9))
    confident = mod["is_target"].to_numpy() & (_pep(g_t, g_0, pi0)(mod["gain"]) < 0.05)
    counts = mod.loc[confident, "explanation"].value_counts()
    masses = {}
    mass = lambda pf: masses.setdefault(pf, _residue_masses(pf))  # noqa: E731
    real = dict(tuple(d[~d["decoy_site"] & ~d["original"]].groupby("group")))
    for _, top in mod.iterrows():
        cand = real.get(top["group"], d.iloc[:0])
        alts = {
            e for e, pf in zip(cand["explanation"], cand["peptidoform"])
            if e != top["explanation"] and _same_fragments(mass(pf), mass(top["peptidoform"]), tol)
        }
        best = max([top["explanation"]] + sorted(alts), key=lambda e: counts.get(e, 0))
        if best != top["explanation"]:
            d.loc[cand.index[cand["explanation"] == best], "explanation_score"] = top["s"] + 1e-6


def rank_sites(
    psm_list: PSMList,
    decoy_sites: PSMList,
    config: dict,
    output_file_root: str | None = None,
    peaks: dict | None = None,
) -> PSMList:
    """
    Rerank the candidate explanations of each spectrum after rescoring.

    Two steps, both learning from within-group feature differences (group-level features cancel;
    precursor-mass features and the search engine score are excluded, see
    ``RANKER_EXCLUDED_FEATURES``):

    1. explanation: the original search engine hit and each modification set of the spectrum,
       each represented by its best candidate by rescoring score, are ranked against each other
       (:py:func:`_rank_explanations`; the top explanation is the positive, the others negatives);
    2. site: within each modification set, the candidates are ranked with the mumble decoy sites
       (the same modification on a residue it cannot occupy) as known negatives
       (:py:func:`ristretto.rank_within_groups`).

    Candidates are ordered by explanation score, then site score. The group's rescoring scores,
    q-values and PEPs are reassigned in that order, so the spectrum keeps its FDR status but
    reports the best-supported explanation. Decoy sites are not returned. Target and decoy
    peptides never share a group; only target groups are used for training.

    Every returned PSM gets ``metadata["site_rank"]`` (1 = best in its group), ``"site_score"``
    (step-2 log-odds, 0 for the original hit) and the localisation outcome of its spectrum
    (:py:func:`ms2rescore.localisation.localise`): ``"localisation_status"``,
    ``"localisation_set"`` and ``"localisation_alternatives"`` (peptidoforms joined by ``;``),
    ``"localisation_pep"``, ``"modification_pep"``, ``"localisation_confidence"`` and
    ``"localisation_qvalue"``. PSMs without competitors get rank 1 and status ``"unranked"``.

    Parameters
    ----------
    psm_list
        Rescored PSMs (scores, q-values and PEPs set), without decoy sites.
    decoy_sites
        Decoy-site PSMs with rescoring features, set aside before rescoring.
    config
        MS²Rescore configuration.
    output_file_root
        If given, the learned feature weights are written to
        ``<output_file_root>.explanation_weights.tsv`` and ``<output_file_root>.site_weights.tsv``.
    peaks
        ``(run, spectrum_id) -> (mz, intensity)`` of the spectra, for the site sets.

    """
    for psm in psm_list:  # defaults: no competitor; every PSM gets every key (one TSV column each)
        psm.metadata.update(
            {"site_rank": "1", "site_score": "0.0000", "localisation_status": "unranked"}
            | dict.fromkeys(LOCALISATION_KEYS, "")
        )
    if not len(decoy_sites):
        logger.warning(
            "Site ranking requested but no mumble decoy sites are present; enable "
            "`include_mumble_decoys` in the mumble configuration. Skipping."
        )
        return psm_list

    all_psms = PSMList(psm_list=list(psm_list) + list(decoy_sites))
    n_scored = len(psm_list)
    feature_names = {f for psm in all_psms for f in psm.rescoring_features}
    df = _build_features_dataframe(all_psms, feature_names, False)
    df["decoy_site"] = np.r_[np.zeros(n_scored, bool), np.ones(len(decoy_sites), bool)]
    df["is_target"] = ~df["is_decoy"].fillna(False).astype(bool)
    df["modset"] = (
        df["peptidoform"].str.findall(r"\[([^\]]+)\]").apply(lambda m: "+".join(sorted(m)))
    )
    # One group per spectrum; target and decoy peptides never compete with each other.
    df["group"] = pd.factorize(
        df["run"].astype(str)
        + "|"
        + df["spectrum_id"].astype(str)
        + "|"
        + df["is_decoy"].astype(str)
    )[0]
    in_group = df.groupby("group")["group"].transform("size") >= 2
    if not in_group.any():
        logger.warning("No spectra with multiple candidates found; skipping site ranking.")
        return psm_list

    feature_cols = sorted(feature_names - RANKER_EXCLUDED_FEATURES)
    d = df[in_group].copy()
    d["original"] = [_is_original_psm(all_psms[i]) for i in d.index]
    d["explanation"] = np.where(d["original"], "original", d["modset"])

    # Step 1: which explanation (original hit or modification set) fits the spectrum best,
    # judged on the best real candidate of each explanation (by rescoring score); no decoy sites.
    rep = d[~d["decoy_site"]].sort_values("score", ascending=False)
    rep = rep.drop_duplicates(["group", "explanation"])
    s1, w1 = _rank_explanations(rep, feature_cols)
    s1 = dict(zip(zip(rep["group"], rep["explanation"]), s1))
    # an explanation with only decoy sites left (no real candidate) is never chosen
    d["explanation_score"] = [s1.get(k, -1e9) for k in zip(d["group"], d["explanation"])]
    _prefer_common_isobars(d, rep, float(config.get("tolerance_value", 0.02)))

    # Step 2: which site within each modification set, with the decoy sites as negatives
    mod = d[~d["original"]].assign(site_group=lambda x: x["group"].astype(str) + "|" + x["explanation"])
    if mod.empty:
        logger.warning("No mumble candidates found; skipping site ranking.")
        return psm_list
    result = ristretto.rank_within_groups(
        mod,
        group_col="site_group",
        negative_col="decoy_site",
        initial_score_col="score",
        feature_cols=feature_cols,
        train_col="is_target",
    )
    d["site_score"] = result.scores["score"].reindex(d.index).fillna(0.0)
    order = np.lexsort((-d["site_score"], -d["explanation_score"], d["group"]))
    d["rank"] = 0
    d.iloc[order, d.columns.get_loc("rank")] = (
        np.arange(len(d)) - np.searchsorted(d["group"].to_numpy()[order], d["group"].to_numpy()[order])
    ) + 1
    logger.info(
        f"Site ranking: {len(d)} candidates in {d['group'].nunique()} groups; "
        f"top candidate is a decoy site in {result.negative_top_rate:.1%} of modification sets."
    )
    for name, weights in [("explanation", w1), ("site", result.feature_weights)]:
        weights = weights.assign(mean=lambda w: w.mean(axis=1))
        top = weights["mean"].reindex(weights["mean"].abs().sort_values(ascending=False).index)
        logger.info(
            f"{name.capitalize()} model feature weights (mean over folds, centred within group): "
            + ", ".join(f"{f} {v:+.2f}" for f, v in top.head(8).items())
        )
        if output_file_root:
            weights.sort_values("mean", key=abs, ascending=False).to_csv(
                f"{output_file_root}.{name}_weights.tsv", sep="\t"
            )

    chosen = d.loc[d.groupby("group")["explanation_score"].idxmax(), ["group", "explanation"]]
    chosen = dict(zip(chosen["group"], chosen["explanation"]))
    grp = d[["group", "decoy_site", "rank", "site_score"]]
    real = grp[~grp["decoy_site"]]
    cand = d[["group", "run", "spectrum_id", "peptidoform", "modset", "decoy_site", "original"]]
    cand = cand.assign(
        score=d["site_score"],
        chosen=[chosen[g] == e for g, e in zip(d["group"], d["explanation"])],
        is_decoy=d["is_decoy"].fillna(False).astype(bool),
        charge=[all_psms[i].peptidoform.precursor_charge or 2 for i in d.index],
        hyperscore=d["hyperscore"] if "hyperscore" in d else np.nan,
    )
    loc = localise(cand, peaks or {}, float(config.get("tolerance_value", 0.02)))
    called = loc[~loc["is_decoy"]]["status"].value_counts(normalize=True)
    logger.info(
        "Localisation outcome of target spectra: "
        + ", ".join(f"{k} {v:.1%}" for k, v in called.items())
    )
    for idx, psm in zip(grp.index, (all_psms[i] for i in grp.index)):
        psm.metadata["site_rank"] = str(int(grp.at[idx, "rank"]))
        psm.metadata["site_score"] = f"{grp.at[idx, 'site_score']:.4f}"
        out = loc.loc[grp.at[idx, "group"]]
        psm.metadata["localisation_status"] = out["status"]
        psm.metadata["localisation_set"] = ";".join(out["site_set"])
        psm.metadata["localisation_alternatives"] = ";".join(out["alternatives"])
        for key, col in [("localisation_pep", "localisation_pep"), ("modification_pep", "modification_pep"),
                         ("localisation_confidence", "confidence"), ("localisation_qvalue", "qvalue")]:
            psm.metadata[key] = "" if pd.isna(out[col]) else f"{out[col]:.4f}"

    # Reassign the rescoring outcome within each group in ranker order (scored PSMs only)
    score, qvalue, pep = (np.asarray(psm_list[c], dtype=float) for c in ("score", "qvalue", "pep"))
    for _, members in real.groupby("group"):
        idx = members.index.to_numpy()
        by_rank = idx[np.argsort(members["rank"].to_numpy())]
        by_svm = idx[np.argsort(-score[idx])]
        score[by_rank], qvalue[by_rank], pep[by_rank] = score[by_svm], qvalue[by_svm], pep[by_svm]
    psm_list["score"], psm_list["qvalue"], psm_list["pep"] = score, qvalue, pep
    psm_list.set_ranks(lower_score_better=False)
    return psm_list
