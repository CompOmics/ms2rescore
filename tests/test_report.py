"""Tests for in-memory report data assembly and HTML generation."""

import psm_utils.io
import pytest
from psm_utils import PSM, PSMList

from ms2rescore.report.data import ReportData
from ms2rescore.report.generate import generate_report
from ms2rescore.report.utils import build_stat_card, compute_protein_stats

FEATURE_NAMES = {"basic": {"feature_a", "feature_b"}}


@pytest.fixture
def psm_list():
    """PSMList with before/after scores, q-values, and rescoring features (10 targets, 4 decoys)."""
    residues = "ACDEFGHIKLMNPQ"  # one unique residue per PSM
    psms = []
    for i in range(14):
        is_decoy = i >= 10
        score_after = float(20 - i)
        qvalue = 0.001 * (i + 1)
        psm = PSM(
            peptidoform=f"PEPTIDE{residues[i]}K/2",
            spectrum_id=str(i),
            run="run1",
            is_decoy=is_decoy,
            score=score_after,
            qvalue=qvalue,
            pep=0.5,
            rescoring_features={"feature_a": float(i), "feature_b": float(i * 2)},
        )
        psm.provenance_data = {
            "before_rescoring_score": float(10 - i),
            "before_rescoring_qvalue": 0.002 * (i + 1),
        }
        psms.append(psm)
    return PSMList(psm_list=psms)


def test_from_run_builds_dataframe_in_memory(psm_list):
    data = ReportData.from_run(psm_list, feature_names=FEATURE_NAMES)

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
    assert data.protein_stats is None  # No fasta provided


def test_from_run_infers_feature_names_when_absent(psm_list):
    data = ReportData.from_run(psm_list, feature_names=None)

    inferred = {name for names in data.feature_names.values() for name in names}
    assert inferred == {"feature_a", "feature_b"}


def test_generate_report_writes_html_without_input_files(psm_list, tmp_path):
    prefix = str(tmp_path / "test.ms2rescore")
    data = ReportData.from_run(psm_list, feature_names=FEATURE_NAMES)

    generate_report(prefix, data)

    report_file = tmp_path / "test.ms2rescore.report.html"
    assert report_file.is_file()
    content = report_file.read_text(encoding="utf-8")
    for tab_title in ["Overview", "Target/decoy evaluation", "Rescoring features"]:
        assert tab_title in content


def test_from_files_roundtrip(psm_list, tmp_path):
    prefix = str(tmp_path / "test.ms2rescore")
    psm_utils.io.write_file(psm_list, prefix + ".psms.tsv", filetype="tsv")

    data = ReportData.from_files(prefix)

    assert len(data.psm_df) == 14
    assert {"score_before", "score_after"}.issubset(set(data.psm_df.columns))


def test_from_files_missing_psm_file_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        ReportData.from_files(str(tmp_path / "does_not_exist"))


def test_compute_protein_stats_without_fasta_returns_none(psm_list):
    data = ReportData.from_run(psm_list, feature_names=FEATURE_NAMES)
    assert compute_protein_stats(data.psm_df, fasta_file=None) is None


def test_build_stat_card_reports_increase():
    card = build_stat_card("PSMs", "psms", before=100, after=150)
    assert card["number"] == 150
    assert card["diff"] == "(+50)"
    assert card["is_increase"] is True
    assert card["card_color"] == "card-bg-blue"
