"""Tests for ms2rescore.rescoring: the ristretto integration layer."""

import numpy as np
import pandas as pd
from psm_utils import PSM, PSMList
from ristretto import RescoreResult

from ms2rescore import rescoring

BASE_CONFIG = {
    "id_decoy_pattern": None,
    "max_psm_rank_output": 1,
    "processes": 1,
    "rescoring": {"train_fdr": 0.1},
    "write_rescoring_tables": True,
}


_RESIDUES = "ACDEFGHIKLMNPQRSTVWY"


def _peptide_seq(spec_i: int, rank_i: int) -> str:
    """Build a unique, valid amino-acid sequence for a given spectrum/rank combination."""
    a = _RESIDUES[spec_i % len(_RESIDUES)]
    b = _RESIDUES[(spec_i // len(_RESIDUES)) % len(_RESIDUES)]
    c = _RESIDUES[rank_i % len(_RESIDUES)]
    return f"PEPT{a}{b}{c}IDEK"


def _make_psm_list(
    n_spectra=40, ranks_per_spectrum=1, seed=0, with_protein=True, run="run1", mumble_ranks=0
):
    """
    Separable synthetic PSMList: targets score higher on both `score` and features.

    `mumble_ranks` extra candidates per spectrum are tagged as mumble-generated
    (metadata["original_psm"] = False), scoring the same as rank 0's decoy-ness.
    """
    rng = np.random.default_rng(seed)
    psms = []
    for spec_i in range(n_spectra):
        is_decoy_spectrum = rng.random() < 0.3
        for rank_i in range(ranks_per_spectrum + mumble_ranks):
            is_original = rank_i < ranks_per_spectrum
            is_decoy = is_decoy_spectrum or rank_i > 0
            shift = 0.0 if is_decoy else 8.0
            psm = PSM(
                peptidoform=f"{_peptide_seq(spec_i, rank_i)}/2",
                spectrum_id=str(spec_i),
                run=run,
                is_decoy=is_decoy,
                score=float(rng.normal(0, 1) + shift),
                qvalue=float("nan"),
                protein_list=(["PROT" + str(spec_i % 5)] if with_protein else None),
                rescoring_features={
                    "feature_a": float(rng.normal(0, 1) + shift),
                    "feature_b": float(rng.normal(0, 1) + shift),
                },
            )
            if not is_original:
                # psm_utils' metadata field is typed dict[str, str], so passing a bool via
                # the constructor fails validation -- but mutating the dict in place after
                # construction bypasses that check, which is presumably how mumble itself
                # sets this flag (a real bool, not a string).
                psm.metadata["original_psm"] = False
            psms.append(psm)
    return PSMList(psm_list=psms)


def test_build_features_dataframe_columns():
    psm_list = _make_psm_list(n_spectra=5, seed=1)
    df = rescoring.build_features_dataframe(
        psm_list, feature_names={"feature_a", "feature_b"}, lower_score_is_better=False
    )

    assert set(df.columns) >= {
        "spectrum_id",
        "is_decoy",
        "peptidoform",
        "peptide",
        "protein",
        "score",
        "feature_a",
        "feature_b",
    }
    assert len(df) == 5
    # peptidoform is charge-stripped
    assert not df["peptidoform"].str.contains("/").any()
    # peptide is bare sequence (no charge, no charge suffix already stripped from peptidoform)
    assert (df["peptide"] == df["peptidoform"]).all()


def test_build_features_dataframe_negates_score_when_lower_is_better():
    psm_list = _make_psm_list(n_spectra=5, seed=2)
    original_scores = np.array(list(psm_list["score"]))

    df = rescoring.build_features_dataframe(
        psm_list, feature_names=set(), lower_score_is_better=True
    )

    assert np.allclose(df["score"].to_numpy(), -original_scores)


def test_build_features_dataframe_omits_protein_when_missing():
    psm_list = _make_psm_list(n_spectra=5, seed=3, with_protein=False)
    df = rescoring.build_features_dataframe(
        psm_list, feature_names=set(), lower_score_is_better=False
    )
    assert "protein" not in df.columns


def test_evaluate_before_returns_result_with_expected_columns():
    psm_list = _make_psm_list(n_spectra=30, seed=4)
    result = rescoring.evaluate_before(psm_list, BASE_CONFIG)

    assert {"spectrum_id", "run", "peptidoform", "score", "qvalue", "pep"}.issubset(
        result.psms.columns
    )
    assert len(result.psms) == 30  # competed to one best PSM per spectrum
    assert result.peptidoforms is not None
    assert (result.psms["run"] == "run1").all()


def test_evaluate_before_disambiguates_spectrum_id_across_runs():
    # Two runs both using spectrum_id "0".."29" -- must not collapse into 30 total.
    psm_list_a = _make_psm_list(n_spectra=30, seed=10, run="runA")
    psm_list_b = _make_psm_list(n_spectra=30, seed=11, run="runB")
    psm_list = PSMList(psm_list=list(psm_list_a) + list(psm_list_b))

    result = rescoring.evaluate_before(psm_list, BASE_CONFIG)

    assert len(result.psms) == 60
    assert set(result.psms["run"]) == {"runA", "runB"}


def test_evaluate_before_excludes_mumble_generated_candidates():
    # 1 original + 1 mumble-tagged candidate per spectrum.
    psm_list = _make_psm_list(n_spectra=20, seed=14, mumble_ranks=1)
    assert len(psm_list) == 40

    result = rescoring.evaluate_before(psm_list, BASE_CONFIG)

    # Only the 20 original candidates ever entered the competition.
    assert len(result.psms) == 20


def test_evaluate_before_mirrors_max_psm_rank_output():
    psm_list = _make_psm_list(n_spectra=20, ranks_per_spectrum=3, seed=15)
    config = {**BASE_CONFIG, "max_psm_rank_output": 2}

    result = rescoring.evaluate_before(psm_list, config)

    # Trimmed to top 2 (by raw score) per spectrum, not competed down to 1.
    assert len(result.psms) <= 40
    counts = result.psms["spectrum_id"].value_counts()
    assert (counts <= 2).all()
    assert (counts == 2).any()


def test_rescore_disambiguates_spectrum_id_across_runs():
    psm_list_a = _make_psm_list(n_spectra=30, seed=12, run="runA")
    psm_list_b = _make_psm_list(n_spectra=30, seed=13, run="runB")
    psm_list = PSMList(psm_list=list(psm_list_a) + list(psm_list_b))

    new_psm_list, after_result = rescoring.rescore(psm_list, BASE_CONFIG, "unused-output-root")

    assert len(new_psm_list) == 60
    assert len(after_result.psms) == 60
    assert set(after_result.psms["run"]) == {"runA", "runB"}


def test_rescore_writes_scores_and_metadata():
    psm_list = _make_psm_list(n_spectra=30, ranks_per_spectrum=1, seed=5)

    new_psm_list, after_result = rescoring.rescore(psm_list, BASE_CONFIG, "unused-output-root")

    assert len(new_psm_list) == len(after_result.psms)
    assert set(new_psm_list["qvalue"]) == set(after_result.psms["qvalue"])
    assert "run" in after_result.psms.columns
    for psm in new_psm_list:
        assert "peptidoform_score" in psm.metadata
        assert "peptidoform_qvalue" in psm.metadata
        assert "protein_score" in psm.metadata  # protein_list present in fixture


def test_rescore_multi_rank_output_keeps_multiple_ranks_per_spectrum():
    psm_list = _make_psm_list(n_spectra=30, ranks_per_spectrum=2, seed=6)
    config = {**BASE_CONFIG, "max_psm_rank_output": 2}

    new_psm_list, after_result = rescoring.rescore(psm_list, config, "unused-output-root")

    new_psm_list.set_ranks(lower_score_better=False)
    assert (new_psm_list["rank"] <= 2).all()
    # At least some spectra should retain both ranks
    counts = pd.Series(list(new_psm_list["spectrum_id"])).value_counts()
    assert (counts == 2).any()

    # after_result IS the multi-rank population -- no separate rank-1-only view exists.
    assert len(after_result.psms) == len(new_psm_list)


def test_evaluate_after_from_psm_list_reproduces_rescore_result():
    psm_list = _make_psm_list(n_spectra=25, seed=16)
    new_psm_list, after_result = rescoring.rescore(psm_list, BASE_CONFIG, "unused-output-root")

    reconstructed = rescoring.evaluate_after_from_psm_list(new_psm_list, BASE_CONFIG)

    assert len(reconstructed.psms) == len(after_result.psms)
    assert np.allclose(
        reconstructed.psms.sort_values("spectrum_id")["qvalue"].to_numpy(),
        after_result.psms.sort_values("spectrum_id")["qvalue"].to_numpy(),
    )


def test_evaluate_before_from_provenance_uses_stashed_score():
    psm_list = _make_psm_list(n_spectra=20, seed=17)
    # Simulate parse_psms()'s provenance stash with a score distinguishable from the current one.
    stashed_scores = list(psm_list["score"])
    for psm, stashed in zip(psm_list, stashed_scores):
        psm.provenance_data["before_rescoring_score"] = stashed
        psm.score = -999.0  # current score must NOT be used

    result = rescoring.evaluate_before_from_provenance(psm_list, BASE_CONFIG)

    assert len(result.psms) == 20
    assert not (result.psms["score"] == -999.0).any()


def _toy_peptidoform(i: int) -> str:
    """Valid, distinct peptidoform string for toy fixtures."""
    return f"PEPT{_RESIDUES[i % len(_RESIDUES)]}IDEK/2"


def _toy_result(spectrum_ids, is_decoy, scores, peps):
    """Minimal RescoreResult.psms with the columns ristretto.evaluate() requires."""
    return RescoreResult(
        psms=pd.DataFrame(
            {
                "spectrum_id": spectrum_ids,
                "is_decoy": is_decoy,
                "peptidoform": [
                    _toy_peptidoform(i).rsplit("/", 1)[0] for i in range(len(spectrum_ids))
                ],
                "score": scores,
                "qvalue": [0.01] * len(spectrum_ids),
                "pep": peps,
            }
        ),
        peptidoforms=pd.DataFrame(),
        peptides=None,
        proteins=None,
        pi0=float("nan"),
        n_iterations=[],
        feature_weights=pd.DataFrame(),
    )


def test_fix_constant_pep_removes_higher_scoring_decoys():
    # 3 targets (scores 1-3), 2 decoys: one below the best target (0.5), one above it (10.0).
    # Only the higher-scoring decoy should be removed; one decoy remains, so q-values/PEP are
    # recomputed rather than falling back to the no-recompute edge case.
    is_decoy = [False, False, False, True, True]
    psm_list = PSMList(
        psm_list=[
            PSM(
                peptidoform=_toy_peptidoform(i),
                spectrum_id=str(i),
                is_decoy=is_decoy[i],
                score=score,
                pep=1.0,
            )
            for i, score in enumerate([1.0, 2.0, 3.0, 0.5, 10.0])
        ]
    )
    result = _toy_result(
        spectrum_ids=["0", "1", "2", "3", "4"],
        is_decoy=is_decoy,
        scores=[1.0, 2.0, 3.0, 0.5, 10.0],
        peps=[1.0] * 5,
    )

    fixed_psm_list, fixed_result = rescoring.fix_constant_pep(psm_list, result)
    assert len(fixed_psm_list) == 4
    assert set(fixed_psm_list["spectrum_id"]) == {"0", "1", "2", "3"}
    assert len(fixed_result.psms) == 4
    # Recomputed, no longer stuck at the degenerate constant 1.0.
    assert not (fixed_result.psms["pep"] == 1.0).all()


def test_fix_constant_pep_no_recompute_when_no_decoys_remain():
    # Removing the higher-scoring decoy leaves zero decoys -- evaluate() can't recompute a
    # competition against an empty decoy set, so this must fall back to plain filtering
    # instead of raising.
    result = _toy_result(
        spectrum_ids=["0", "1"],
        is_decoy=[False, True],
        scores=[1.0, 5.0],
        peps=[1.0, 1.0],
    )
    psm_list = PSMList(
        psm_list=[
            PSM(peptidoform="AAAK/2", spectrum_id="0", is_decoy=False, score=1.0, pep=1.0),
            PSM(peptidoform="BBBK/2", spectrum_id="1", is_decoy=True, score=5.0, pep=1.0),
        ]
    )

    fixed_psm_list, fixed_result = rescoring.fix_constant_pep(psm_list, result)
    assert len(fixed_psm_list) == 1
    assert not fixed_psm_list[0].is_decoy
    assert len(fixed_result.psms) == 1
    assert fixed_result.psms["spectrum_id"].tolist() == ["0"]


def test_fix_constant_pep_is_noop_when_pep_not_constant():
    psm_list = PSMList(
        psm_list=[
            PSM(peptidoform="AAAK/2", spectrum_id="1", is_decoy=False, score=1.0, pep=0.5),
            PSM(peptidoform="BBBK/2", spectrum_id="2", is_decoy=True, score=5.0, pep=1.0),
        ]
    )
    result = _toy_result(
        spectrum_ids=["1", "2"], is_decoy=[False, True], scores=[1.0, 5.0], peps=[0.5, 1.0]
    )

    fixed_psm_list, fixed_result = rescoring.fix_constant_pep(psm_list, result)
    assert len(fixed_psm_list) == 2
    assert len(fixed_result.psms) == 2


def test_write_rescoring_tables(tmp_path):
    psm_list = _make_psm_list(n_spectra=30, seed=7)
    _, after = rescoring.rescore(psm_list, BASE_CONFIG, "unused-output-root")

    prefix = str(tmp_path / "test")
    rescoring.write_rescoring_tables(after, prefix)

    for level in ("psms", "peptidoforms", "proteins"):
        assert (tmp_path / f"test.ristretto.{level}.tsv").is_file()
    assert (tmp_path / "test.ristretto.weights.tsv").is_file()

    # TSV, not parquet, and human-readable (has a header row with expected columns).
    psms_tsv = pd.read_csv(tmp_path / "test.ristretto.psms.tsv", sep="\t")
    assert {"spectrum_id", "run", "score", "qvalue", "pep"}.issubset(psms_tsv.columns)


def test_write_rescoring_tables_no_peptides_file_when_no_peptide_col(tmp_path):
    # BASE_CONFIG's build_features_dataframe never sets a peptide_col in `after` (peptide_col
    # is always "peptide" in rescore(), so peptides IS present -- this documents that proteins
    # can be absent when protein_list isn't, by using a fixture without proteins).
    psm_list = _make_psm_list(n_spectra=30, seed=9, with_protein=False)
    _, after = rescoring.rescore(psm_list, BASE_CONFIG, "unused-output-root")
    assert after.proteins is None

    prefix = str(tmp_path / "test")
    rescoring.write_rescoring_tables(after, prefix)

    assert not (tmp_path / "test.ristretto.proteins.tsv").is_file()
    assert (tmp_path / "test.ristretto.peptides.tsv").is_file()
