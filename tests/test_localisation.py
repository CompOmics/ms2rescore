"""Tests for ms2rescore.localisation.localise()."""

import numpy as np
import pandas as pd

from ms2rescore.localisation import _fragments, _residue_masses, localise

S4, S5, DECOY = "PEPS[Phospho]SIDEK", "PEPSS[Phospho]IDEK", "P[Phospho]EPSSIDEK"


def _ions(pf):
    return _fragments(_residue_masses(pf), 2)


def test_localise_reports_site_sets_from_peaks():
    shared = np.intersect1d(np.round(_ions(S4), 4), np.round(_ions(S5), 4))
    s4_only = np.setdiff1d(np.round(_ions(S4), 4), np.round(_ions(S5), 4))
    peaks = {
        ("r", "ambiguous"): (np.sort(shared), np.ones(len(shared))),  # no site-determining ion
        ("r", "s4"): (np.sort(np.r_[shared, s4_only[:2]]), np.ones(len(shared) + 2)),
    }
    rows = []
    for sid in ("ambiguous", "s4", "original"):
        top_score = -1.0 if sid == "original" else 2.0
        for pf, decoy_site, original, score, hs in [
            ("PEPSSIDEK", False, True, 0.0, 10.0),
            (S4, False, False, top_score, 20.0),
            (S5, False, False, top_score - 0.1, 19.0),
            (DECOY, True, False, -3.0, 12.0),
        ]:
            rows.append(dict(run="r", spectrum_id=sid, peptidoform=pf, decoy_site=decoy_site,
                             original=original, score=score, hyperscore=hs, is_decoy=False))
    for i in range(20):  # decoy-peptide spectra: no gain over their original
        rows += [dict(run="r", spectrum_id=f"d{i}", peptidoform=pf, decoy_site=False, original=o,
                      score=s, hyperscore=10.0, is_decoy=True)
                 for pf, o, s in [("KEDISSPEP", True, 0.0), ("KEDIS[Phospho]SPEP", False, 1.0)]]
    df = pd.DataFrame(rows).assign(charge=2)
    df["modset"] = np.where(df.peptidoform.str.contains("Phospho"), "Phospho", "")
    df["group"] = pd.factorize(df.spectrum_id)[0]

    out = localise(df, peaks)
    out.index = df.groupby("group").spectrum_id.first().reindex(out.index)

    assert sorted(out.loc["ambiguous", "site_set"]) == sorted([S4, S5])
    assert out.loc["s4", "site_set"] == [S4]
    assert out.loc["original", "status"] == "original"
    assert out.loc["s4", "modification_pep"] < 0.5
