"""End-to-end tests for the full ``core.rescore`` pipeline (PSM file in, ristretto rescoring,
report/flashlfq output), plus intermediate-output recovery on mid-pipeline errors."""

import json

import numpy as np
import psm_utils.io
import pytest
from psm_utils import PSM, PSMList

from ms2rescore import core
from ms2rescore.config_parser import parse_configurations
from ms2rescore.feature_generators.base import FeatureGeneratorBase

_RESIDUES = "ACDEFGHIKLMNPQRSTVWY"


def _peptide_seq(i: int) -> str:
    a = _RESIDUES[i % len(_RESIDUES)]
    b = _RESIDUES[(i // len(_RESIDUES)) % len(_RESIDUES)]
    c = _RESIDUES[(i // len(_RESIDUES) ** 2) % len(_RESIDUES)]
    return f"PEPT{a}{b}{c}IDEK"


def _make_psm_list(n_spectra=200, seed=0) -> PSMList:
    """Separable synthetic PSMList: targets score higher than decoys, so rescoring should
    recover more identifications at a fixed FDR than the raw search engine score."""
    rng = np.random.default_rng(seed)
    psms = []
    for i in range(n_spectra):
        is_decoy = rng.random() < 0.3
        shift = 0.0 if is_decoy else 8.0
        psms.append(
            PSM(
                peptidoform=f"{_peptide_seq(i)}/2",
                spectrum_id=str(i),
                run="run1",
                is_decoy=is_decoy,
                score=float(rng.normal(0, 1) + shift),
                precursor_mz=500.0 + i,
                protein_list=["PROT" + str(i % 5)],
                rescoring_features={},
            )
        )
    return PSMList(psm_list=psms)


def _write_psm_file(tmp_path, psm_list) -> str:
    psm_path = str(tmp_path / "psms.tsv")
    psm_utils.io.write_file(psm_list, psm_path, filetype="tsv")
    return psm_path


def _make_config(tmp_path, psm_path, **overrides):
    ms2rescore_config = {
        "psm_file": psm_path,
        "psm_file_type": "tsv",
        "output_path": str(tmp_path / "out"),
        "feature_generators": {"basic": {}},
        "rescoring": {"train_fdr": 0.1, "model": "svm"},
        "write_report": False,
        "write_flashlfq": False,
        "log_level": "info",
    }
    ms2rescore_config.update(overrides)
    return parse_configurations([{"ms2rescore": ms2rescore_config}])


class _FailingFeatureGenerator(FeatureGeneratorBase):
    """Stands in for a real feature generator, but always fails in ``add_features``."""

    @property
    def feature_names(self):
        return ["broken_feature"]

    def add_features(self, psm_list):
        raise RuntimeError("simulated feature generator failure")


def test_rescore_end_to_end_writes_expected_outputs(tmp_path):
    psm_path = _write_psm_file(tmp_path, _make_psm_list())
    config = _make_config(tmp_path, psm_path, write_report=True)

    core.rescore(config)

    output_root = tmp_path / "out"

    # Final PSM output, with rescored q-values.
    result_psm_list = psm_utils.io.read_file(str(output_root) + ".tsv", filetype="tsv")
    assert len(result_psm_list) == 200
    assert all(psm.qvalue is not None for psm in result_psm_list)

    # Separable synthetic data: rescoring must recover a substantial number of IDs at 1% FDR.
    n_identified = sum(1 for psm in result_psm_list if not psm.is_decoy and psm.qvalue <= 0.01)
    assert n_identified > 0

    # Full configuration (defaults + overrides) is dumped alongside the results.
    full_config = json.loads((tmp_path / "out.full-config.json").read_text())
    assert full_config["ms2rescore"]["feature_generators"] == {"basic": {}}

    # Feature names are logged per generator.
    feature_names_content = (tmp_path / "out.feature_names.tsv").read_text()
    assert "basic\tsearch_engine_score" in feature_names_content

    # Rescoring tables (from ristretto's RescoreResult) are written.
    for suffix in ("psms", "peptidoforms", "proteins", "weights"):
        assert (tmp_path / f"out.{suffix}.tsv").is_file()

    # Report generation was requested and must produce an HTML report.
    report_path = tmp_path / "out.report.html"
    assert report_path.is_file()
    assert "<html" in report_path.read_text().lower()


def test_rescore_end_to_end_writes_flashlfq_when_enabled(tmp_path):
    psm_path = _write_psm_file(tmp_path, _make_psm_list(seed=1))
    config = _make_config(tmp_path, psm_path, write_flashlfq=True)

    core.rescore(config)

    flashlfq_path = tmp_path / "out.flashlfq.tsv"
    assert flashlfq_path.is_file()
    # write_file(..., only_target=True) must exclude decoys from the FlashLFQ export.
    flashlfq_psm_list = psm_utils.io.read_file(str(flashlfq_path), filetype="flashlfq")
    assert len(flashlfq_psm_list) > 0
    assert not any(psm.is_decoy for psm in flashlfq_psm_list)


def test_rescore_writes_intermediate_output_on_feature_generator_error(tmp_path, monkeypatch):
    monkeypatch.setitem(core.FEATURE_GENERATORS, "basic", _FailingFeatureGenerator)

    psm_path = _write_psm_file(tmp_path, _make_psm_list(seed=2))
    config = _make_config(tmp_path, psm_path)

    with pytest.raises(RuntimeError, match="simulated feature generator failure"):
        core.rescore(config)

    intermediate_path = tmp_path / "out.intermediate.tsv"
    assert intermediate_path.is_file()
    intermediate_psm_list = psm_utils.io.read_file(str(intermediate_path), filetype="tsv")
    assert len(intermediate_psm_list) == 200


def test_rescore_writes_intermediate_output_on_rescoring_error(tmp_path, monkeypatch):
    def _raise(*args, **kwargs):
        raise RuntimeError("simulated rescoring failure")

    monkeypatch.setattr(core.rescoring, "rescore", _raise)

    psm_path = _write_psm_file(tmp_path, _make_psm_list(seed=3))
    config = _make_config(tmp_path, psm_path)

    with pytest.raises(RuntimeError, match="simulated rescoring failure"):
        core.rescore(config)

    intermediate_path = tmp_path / "out.intermediate.tsv"
    assert intermediate_path.is_file()
    # Feature generation completed before the (simulated) rescoring failure, so the
    # intermediate output must already carry the basic-generator features.
    intermediate_psm_list = psm_utils.io.read_file(str(intermediate_path), filetype="tsv")
    assert len(intermediate_psm_list) == 200
    assert "search_engine_score" in intermediate_psm_list[0].rescoring_features
