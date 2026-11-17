"""
Localisation statistics for the site ranker output.

Per spectrum, the top-ranked real candidate is either the original PSM or a modified candidate.
For a modified candidate:

- **site set**: the top candidate plus every candidate with the same modification set for which
  the observed fragments give no net support to the top candidate (site-determining ions,
  b/y ions at charge 1 to max(1, z - 1));
- **alternatives**: candidates with another modification set but the same fragment masses;
- **localisation PEP**: target-decoy competition on the margin of the site set over the best
  candidate outside it, against the margin of the best decoy site; isotonic fit;
- **modification PEP**: target-decoy competition on the hyperscore gain of the top candidate over
  the original PSM, with the same gain in decoy-peptide spectra as null;
- **confidence** = (1 - modification PEP) x (1 - localisation PEP), and a q-value from it.

Decoy sites with the same fragment masses as a real candidate do not enter the competition.
"""

import numpy as np
import pandas as pd
from psm_utils import Peptidoform
from pyteomics.mass import std_aa_mass
from scipy.special import logsumexp
from sklearn.isotonic import IsotonicRegression

PROTON = 1.007276466812
WATER = 18.010564684


def _residue_masses(peptidoform: str) -> np.ndarray:
    pf = Peptidoform(peptidoform)
    masses = [
        std_aa_mass.get(aa, 0.0) + sum(float(m.mass) for m in mods or [])
        for aa, mods in pf.parsed_sequence
    ]
    masses[0] += sum(float(m.mass) for m in pf.properties.get("n_term") or [])
    masses[-1] += sum(float(m.mass) for m in pf.properties.get("c_term") or [])
    return np.asarray(masses)


def _fragments(masses: np.ndarray, charge: int) -> np.ndarray:
    """b and y ion m/z, charges 1 to max(1, charge - 1)."""
    cum = np.cumsum(masses)
    neutral = np.r_[cum[:-1], cum[-1] - cum[:-1] + WATER]
    return np.concatenate([(neutral + z * PROTON) / z for z in range(1, max(1, charge - 1) + 1)])


def _unique(a: np.ndarray, b: np.ndarray, tol: float) -> np.ndarray:
    """Ions of a without an ion of b within tol."""
    return a[np.abs(a[:, None] - b[None, :]).min(axis=1) > tol]


def _n_observed(ions: np.ndarray, mz: np.ndarray, tol: float) -> int:
    i = np.clip(np.searchsorted(mz, ions), 1, len(mz) - 1)
    return int((np.minimum(np.abs(mz[i] - ions), np.abs(mz[i - 1] - ions)) <= tol).sum())


def _same_fragments(a: np.ndarray, b: np.ndarray, tol: float) -> bool:
    return len(a) == len(b) and bool(np.all(np.abs(np.cumsum(a) - np.cumsum(b)) <= tol))


def _pep(target, null, pi0=1.0):
    """Isotonic target-decoy PEP for a statistic where higher is better."""
    target, null = np.asarray(target, float), np.asarray(null, float)
    if not len(null):
        return lambda v: np.zeros(len(v))
    iso = IsotonicRegression(increasing=False, y_min=1e-6, y_max=1 - 1e-6, out_of_bounds="clip")
    iso.fit(np.r_[target, null], np.r_[np.zeros(len(target)), np.ones(len(null))])
    ratio = pi0 * len(target) / len(null)
    return lambda v: np.clip(iso.predict(np.asarray(v, float)) / (1 - iso.predict(np.asarray(v, float))) * ratio, 0, 1)


def localise(df: pd.DataFrame, peaks: dict, tol: float = 0.02) -> pd.DataFrame:
    """
    Localisation outcome per group.

    Parameters
    ----------
    df
        One row per candidate with columns ``group``, ``run``, ``spectrum_id``, ``peptidoform``
        (without charge), ``charge``, ``modset``, ``decoy_site``, ``original``, ``is_decoy``,
        ``hyperscore`` and ``score`` (ranker log-odds); optional ``chosen`` (bool) marks the
        explanation that takes the top position before the score is compared.
    peaks
        ``(run, spectrum_id) -> (mz, intensity)``; spectra without peaks get a single-site set.
    tol
        Fragment tolerance in Da.

    Returns
    -------
    pandas.DataFrame
        Indexed by group: ``status`` (``confident``, ``ambiguous``, ``original``, ``decoy``,
        ``no_gain``, ``low_confidence``), ``site_set`` and ``alternatives`` (lists of
        peptidoforms), ``localisation_pep``, ``modification_pep``, ``confidence``, ``qvalue``.

    """
    masses = {pf: _residue_masses(pf) for pf in df["peptidoform"].unique()}
    rows = []
    order = ["chosen", "score"] if "chosen" in df else ["score"]
    for gid, g in df.sort_values(order, ascending=False).groupby("group", sort=False):
        real, top = g[~g["decoy_site"]], None
        if not len(real) or real["original"].iloc[0]:
            rows.append((gid, g["is_decoy"].iloc[0], "original", [], [], np.nan, np.nan, False, np.nan))
            continue
        top = real.iloc[0]
        same = real[(real["modset"] == top["modset"]) & ~real["original"]]
        decoys = g[g["decoy_site"] & (g["modset"] == top["modset"])]
        decoys = decoys[np.array(
            [not any(_same_fragments(masses[d], masses[r], tol) for r in same["peptidoform"])
             for d in decoys["peptidoform"]], dtype=bool
        )]
        site = [True] * len(same)
        spectrum = peaks.get((top["run"], top["spectrum_id"]))
        if spectrum is not None:
            mz = np.sort(np.asarray(spectrum[0], float))
            ions = {pf: _fragments(masses[pf], int(top["charge"])) for pf in same["peptidoform"]}
            t = ions[top["peptidoform"]]
            for i, pf in enumerate(same["peptidoform"].iloc[1:], start=1):
                b = ions[pf]
                site[i] = _n_observed(_unique(t, b, tol), mz, tol) <= _n_observed(_unique(b, t, tol), mz, tol)
        else:
            site = [True] + [False] * (len(same) - 1)
        site = np.asarray(site)
        others = pd.concat([same["score"][~site], decoys["score"]])
        margin = logsumexp(same["score"][site]) - (others.max() if len(others) else 50.0)
        decoy_margin = np.nan
        if len(decoys):
            rest = pd.concat([same["score"], decoys["score"].iloc[1:]])
            decoy_margin = decoys["score"].iloc[0] - rest.max()
        alternatives = [
            pf for pf in real.loc[real["modset"] != top["modset"], "peptidoform"]
            if _same_fragments(masses[pf], masses[top["peptidoform"]], tol)
        ]
        original = g.loc[g["original"], "hyperscore"]
        gain = top["hyperscore"] - original.iloc[0] if len(original) else np.nan
        decoy_wins = len(decoys) > 0 and decoys["score"].iloc[0] > top["score"]
        rows.append((gid, g["is_decoy"].iloc[0], "decoy" if decoy_wins else "", list(same["peptidoform"][site]),
                     alternatives, margin, decoy_margin, len(decoys) == 0 and len(same) == 1, gain))
    r = pd.DataFrame(rows, columns=["group", "is_decoy", "status", "site_set", "alternatives", "margin",
                                    "decoy_margin", "single", "gain"]).set_index("group")

    modified = (r["status"] != "original") & ~r["is_decoy"]
    target = r[modified & r["decoy_margin"].notna()]
    loc = _pep(target["margin"].clip(-50, 50), target["decoy_margin"].dropna().clip(-50, 50))
    r["localisation_pep"] = np.where(r["single"], 0.0, loc(r["margin"].fillna(-50).clip(-50, 50)))

    g_t, g_0 = r.loc[modified, "gain"].dropna(), r.loc[r["is_decoy"], "gain"].dropna()
    pi0 = min(1.0, (g_t <= 0).mean() / (g_0 <= 0).mean()) if len(g_0) and (g_0 <= 0).any() else 1.0
    r["modification_pep"] = np.where(r["gain"].notna(), _pep(g_t, g_0, pi0)(r["gain"].fillna(0)), 1.0)

    r["confidence"] = (1 - r["modification_pep"]) * (1 - r["localisation_pep"])
    r.loc[r["status"] == "original", ["localisation_pep", "modification_pep", "confidence"]] = np.nan
    open_ = r["status"] == ""
    r.loc[open_ & (r["modification_pep"] >= 0.5), "status"] = "no_gain"
    r.loc[(r["status"] == "") & (r["confidence"] < 0.5), "status"] = "low_confidence"
    r.loc[r["status"] == "", "status"] = np.where(r.loc[r["status"] == "", "site_set"].str.len() > 1,
                                                  "ambiguous", "confident")
    # q-value: running mean of (1 - confidence) over target spectra sorted by confidence
    called = r[~r["is_decoy"] & r["confidence"].notna()].sort_values("confidence", ascending=False)
    q = (1 - called["confidence"]).expanding().mean()[::-1].cummin()[::-1]
    r["qvalue"] = q.reindex(r.index)
    return r.drop(columns=["margin", "decoy_margin", "single", "gain"])
