"""Tests for in-memory report data assembly and HTML generation."""

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
import psm_utils.io
import pytest
from psm_utils import PSM, PSMList
from ristretto import RescoreResult

from click.testing import CliRunner

from ms2rescore.report import charts
from ms2rescore.report.__main__ import main as report_main
from ms2rescore.report.data import ReportData
from ms2rescore.report.generate import generate_report
from ms2rescore.report.utils import (
    _n_identified,
    build_stat_card,
    compute_id_stats,
    compute_protein_stats,
)

FEATURE_NAMES = {"basic": {"feature_a", "feature_b"}}
_CHARGE_PATTERN = re.compile(r"(/\d+$)")


def _make_result(psm_list: PSMList, scores, qvalues) -> RescoreResult:
    """Build a RescoreResult whose psms join onto psm_list's run/spectrum_id."""
    psms = pd.DataFrame(
        {
            "spectrum_id": list(psm_list["spectrum_id"]),
            "run": list(psm_list["run"]),
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
            rank=1,
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

    # psm_df is just psm_list + features -- no before/after merge performed here anymore;
    # charts that need before/after go straight to data.before/data.after.
    expected_columns = {"is_decoy", "score", "qvalue", "feature_a", "feature_b"}
    assert expected_columns.issubset(set(data.psm_df.columns))
    assert len(data.psm_df) == 14
    assert set(data.feature_names) == {"basic"}
    assert set(data.feature_names["basic"]) == {"feature_a", "feature_b"}
    assert data.feature_weights is after.feature_weights
    assert data.before is before
    assert data.after is after


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


def _write_regen_fixture(tmp_path, prefix_name="test.ms2rescore"):
    """
    Write a PSM list + full-config.json for a standalone report-regen test.

    No separate rescoring-result tables are written at all -- ``from_files`` reconstructs
    ``before``/``after`` entirely from the PSM list's ``provenance_data`` and current
    score/identity columns.

    """
    residues = "ACDEFGHIKLMNPQ"
    psms = []
    for i in range(14):
        is_decoy = i >= 10
        psm = PSM(
            peptidoform=f"PEPTIDE{residues[i]}K/2",
            spectrum_id=str(i),
            run="run1",
            is_decoy=is_decoy,
            score=float(20 - i),  # "after" (final, rescored) score
            qvalue=0.001 * (i + 1),
            pep=0.5,
            rank=1,
            protein_list=[f"PROT{i % 3}"],
            provenance_data={"before_rescoring_score": str(float(10 - i))},
            rescoring_features={"feature_a": float(i), "feature_b": float(i * 2)},
        )
        psms.append(psm)
    psm_list = PSMList(psm_list=psms)

    prefix = str(tmp_path / prefix_name)
    psm_utils.io.write_file(psm_list, prefix + ".psms.tsv", filetype="tsv")
    config = {
        "ms2rescore": {
            "id_decoy_pattern": None,
            "max_psm_rank_output": 1,
            "rescoring": {"train_fdr": 0.1},
        }
    }
    with open(prefix + ".full-config.json", "w") as f:
        json.dump(config, f)
    return prefix


def test_from_files_reconstructs_before_and_after_without_any_table_files(tmp_path):
    """
    No *.ristretto.*.tsv/parquet files are written at all -- from_files must still fully
    reconstruct before/after (scores, q-values, rollups) from just the PSM list.
    """
    prefix = _write_regen_fixture(tmp_path)

    data = ReportData.from_files(prefix)

    assert len(data.psm_df) == 14
    assert len(data.before.psms) == 14
    assert len(data.after.psms) == 14
    # Reconstructed "before" used provenance_data's score (10-i), not the current one (20-i).
    assert not np.allclose(
        data.before.psms.sort_values("spectrum_id")["score"].to_numpy(),
        data.after.psms.sort_values("spectrum_id")["score"].to_numpy(),
    )
    # Rollups also reconstructed, not just the PSM-level table.
    assert len(data.before.peptidoforms) == 14
    assert len(data.after.peptidoforms) == 14


def test_from_files_missing_psm_file_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        ReportData.from_files(str(tmp_path / "does_not_exist"))


def test_from_run_default_fdr_threshold(psm_list, before_after):
    before, after = before_after
    data = ReportData.from_run(psm_list, feature_names=FEATURE_NAMES, before=before, after=after)
    assert data.fdr_threshold == 0.01


def test_from_run_respects_custom_fdr_threshold(psm_list, before_after):
    before, after = before_after
    data = ReportData.from_run(
        psm_list, feature_names=FEATURE_NAMES, before=before, after=after, fdr_threshold=0.05
    )
    assert data.fdr_threshold == 0.05
    assert data.id_stats == compute_id_stats(before, after, 0.05)


def test_from_files_reads_report_fdr_from_config(tmp_path):
    prefix = _write_regen_fixture(tmp_path)
    with open(prefix + ".full-config.json") as f:
        config = json.load(f)
    config["ms2rescore"]["report_fdr"] = 0.05
    with open(prefix + ".full-config.json", "w") as f:
        json.dump(config, f)

    data = ReportData.from_files(prefix)

    assert data.fdr_threshold == 0.05


def test_from_files_fdr_override_takes_precedence_over_config(tmp_path):
    prefix = _write_regen_fixture(tmp_path)
    with open(prefix + ".full-config.json") as f:
        config = json.load(f)
    config["ms2rescore"]["report_fdr"] = 0.05
    with open(prefix + ".full-config.json", "w") as f:
        json.dump(config, f)

    data = ReportData.from_files(prefix, fdr_threshold=0.2)

    assert data.fdr_threshold == 0.2


def test_cli_fdr_option_overrides_config(tmp_path):
    prefix = _write_regen_fixture(tmp_path)
    psm_file = prefix + ".psms.tsv"

    result = CliRunner().invoke(report_main, [psm_file, "--fdr", "0.2"])

    assert result.exit_code == 0, result.output
    report_file = Path(prefix + ".report.html")
    assert report_file.is_file()


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


def test_n_identified_counts_targets_below_threshold():
    # Mumble-generated candidates are excluded upstream, in rescoring.evaluate_before, so
    # _n_identified itself does no masking -- just qvalue/is_decoy.
    df = pd.DataFrame(
        {"qvalue": [0.001, 0.001, 0.02], "is_decoy": [False, True, False]}
    )
    assert _n_identified(df, 0.01) == 1


_EMPTY_ROLLUP = pd.DataFrame(columns=["peptidoform", "score", "qvalue", "pep", "is_decoy", "n_psms"])


def test_identification_overlap_disambiguates_spectrum_id_by_run():
    """Two runs reusing the same native spectrum_id must not collide in the overlap counts."""
    before = RescoreResult(
        psms=pd.DataFrame(
            {
                "spectrum_id": ["1", "2", "1", "2"],
                "run": ["runA", "runA", "runB", "runB"],
                "is_decoy": [False, False, False, False],
                "qvalue": [0.001, 0.001, 0.001, 0.001],
            }
        ),
        peptidoforms=_EMPTY_ROLLUP,
        peptides=None,
        proteins=None,
        pi0=0.1,
        n_iterations=[],
        feature_weights=pd.DataFrame(),
    )
    # runA's spectrum "1" drops out after rescoring; everything else survives.
    after = RescoreResult(
        psms=pd.DataFrame(
            {
                "spectrum_id": ["2", "1", "2"],
                "run": ["runA", "runB", "runB"],
                "is_decoy": [False, False, False],
                "qvalue": [0.001, 0.001, 0.001],
            }
        ),
        peptidoforms=_EMPTY_ROLLUP,
        peptides=None,
        proteins=None,
        pi0=0.1,
        n_iterations=[],
        feature_weights=pd.DataFrame(),
    )

    fig = charts.identification_overlap(before, after)
    values = {(trace.name, trace.y[0]): trace.x[0] for trace in fig.data}

    # With the bug (bare spectrum_id), runA's "1"/"2" and runB's "1"/"2" would collide into
    # only 2 distinct keys instead of 4, undercounting "removed" and overcounting "retained".
    assert values[("removed", "spectra")] == -1
    assert values[("retained", "spectra")] == 3
    assert values[("gained", "spectra")] == 0
