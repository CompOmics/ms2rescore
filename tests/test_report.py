"""Tests for in-memory report data assembly and HTML generation."""

import re

import pandas as pd
import psm_utils.io
import pytest
from psm_utils import PSM, PSMList
from ristretto import RescoreResult

from ms2rescore.report.data import ReportData
from ms2rescore.report.generate import generate_report
from ms2rescore.report.utils import build_stat_card, compute_id_stats, compute_protein_stats

FEATURE_NAMES = {"basic": {"feature_a", "feature_b"}}
_CHARGE_PATTERN = re.compile(r"(/\d+$)")


def _make_result(psm_list: PSMList, scores, qvalues) -> RescoreResult:
    """Build a RescoreResult whose psms join onto psm_list's spectrum_id/peptidoform."""
    psms = pd.DataFrame(
        {
            "spectrum_id": list(psm_list["spectrum_id"]),
            "is_decoy": list(psm_list["is_decoy"]),
            "peptidoform": [_CHARGE_PATTERN.sub("", str(p)) for p in psm_list["peptidoform"]],
            "score": scores,
            "qvalue": qvalues,
            "pep": [0.1] * len(psm_list),
        }
    )
    peptidoforms = psms.copy()
    peptidoforms["n_psms"] = 1
    return RescoreResult(
        psms=psms,
        peptidoforms=peptidoforms,
        peptides=None,
        proteins=None,
        pi0=0.1,
        n_iterations=[3],
        feature_weights=pd.DataFrame(
            {"fold_1": [0.5, -0.3]}, index=pd.Index(["feature_a", "feature_b"], name="feature")
        ),
    )


@pytest.fixture
def psm_list():
    """PSMList with rescoring features (10 targets, 4 decoys)."""
    residues = "ACDEFGHIKLMNPQ"  # one unique residue per PSM
    psms = []
    for i in range(14):
        is_decoy = i >= 10
        psm = PSM(
            peptidoform=f"PEPTIDE{residues[i]}K/2",
            spectrum_id=str(i),
            run="run1",
            is_decoy=is_decoy,
            score=float(20 - i),
            qvalue=0.001 * (i + 1),
            pep=0.5,
            rescoring_features={"feature_a": float(i), "feature_b": float(i * 2)},
        )
        psms.append(psm)
    return PSMList(psm_list=psms)


@pytest.fixture
def before_after(psm_list):
    n = len(psm_list)
    before = _make_result(
        psm_list,
        scores=[float(10 - i) for i in range(n)],
        qvalues=[0.002 * (i + 1) for i in range(n)],
    )
    after = _make_result(
        psm_list,
        scores=[float(20 - i) for i in range(n)],
        qvalues=[0.001 * (i + 1) for i in range(n)],
    )
    return before, after


def test_from_run_builds_dataframe_in_memory(psm_list, before_after):
    before, after = before_after
    data = ReportData.from_run(psm_list, feature_names=FEATURE_NAMES, before=before, after=after)

    expected_columns = {
        "score_before",
        "score_after",
        "qvalue_before",
        "qvalue_after",
        "is_decoy",
        "feature_a",
        "feature_b",
    }
    assert expected_columns.issubset(set(data.psm_df.columns))
    assert len(data.psm_df) == 14
    assert set(data.feature_names) == {"basic"}
    assert set(data.feature_names["basic"]) == {"feature_a", "feature_b"}
    assert data.feature_weights is after.feature_weights


def test_from_run_infers_feature_names_when_absent(psm_list, before_after):
    before, after = before_after
    data = ReportData.from_run(psm_list, feature_names=None, before=before, after=after)

    inferred = {name for names in data.feature_names.values() for name in names}
    assert inferred == {"feature_a", "feature_b"}


def test_generate_report_writes_html_without_input_files(psm_list, before_after, tmp_path):
    before, after = before_after
    prefix = str(tmp_path / "test.ms2rescore")
    data = ReportData.from_run(psm_list, feature_names=FEATURE_NAMES, before=before, after=after)

    generate_report(prefix, data)

    report_file = tmp_path / "test.ms2rescore.report.html"
    assert report_file.is_file()
    content = report_file.read_text(encoding="utf-8")
    for tab_title in ["Overview", "Target/decoy evaluation", "Rescoring features"]:
        assert tab_title in content


def test_from_files_roundtrip(psm_list, before_after, tmp_path):
    before, after = before_after
    prefix = str(tmp_path / "test.ms2rescore")
    psm_utils.io.write_file(psm_list, prefix + ".psms.tsv", filetype="tsv")
    before.psms.to_parquet(prefix + ".ristretto.psms_before.parquet")
    before.peptidoforms.to_parquet(prefix + ".ristretto.peptidoforms_before.parquet")
    after.psms.to_parquet(prefix + ".ristretto.psms_after.parquet")
    after.peptidoforms.to_parquet(prefix + ".ristretto.peptidoforms_after.parquet")
    after.feature_weights.to_parquet(prefix + ".ristretto.weights.parquet")

    data = ReportData.from_files(prefix)

    assert len(data.psm_df) == 14
    assert {"score_before", "score_after"}.issubset(set(data.psm_df.columns))


def test_from_files_missing_psm_file_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        ReportData.from_files(str(tmp_path / "does_not_exist"))


def test_from_files_without_ristretto_tables_degrades_gracefully(psm_list, tmp_path):
    prefix = str(tmp_path / "test.ms2rescore")
    psm_utils.io.write_file(psm_list, prefix + ".psms.tsv", filetype="tsv")

    data = ReportData.from_files(prefix)

    assert len(data.psm_df) == 14
    assert data.psm_df["score_before"].isna().all()
    assert data.id_stats == []


def test_compute_protein_stats_without_protein_col_returns_none(before_after):
    before, after = before_after
    assert compute_protein_stats(before, after) is None


def test_compute_id_stats_psm_and_peptide_cards(before_after):
    before, after = before_after
    stats = compute_id_stats(before, after)
    items = {s["item"] for s in stats}
    assert "PSMs" in items
    assert "Peptides" in items


def test_build_stat_card_reports_increase():
    card = build_stat_card("PSMs", "psms", before=100, after=150)
    assert card["number"] == 150
    assert card["diff"] == "(+50)"
    assert card["is_increase"] is True
    assert card["card_color"] == "card-bg-blue"
