# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed

- `deeplc` feature generator: `deeplc_retrain` option renamed to `finetune` (redundant naming,
  already scoped under `deeplc`). The old name is still accepted but deprecated, and emits a
  warning.

### Fixed

- Removed stale `ms2_tolerance` option from the default `ms2pip` feature generator configuration,
  schema, and example configs. Fragment tolerance is configured centrally via the top-level
  `tolerance_value`/`tolerance_mode` options; `MS2PIPFeatureGenerator` no longer accepts
  `ms2_tolerance`.
- Documented `model_dir` (`ms2pip` feature generator) in the configuration schema; it already
  worked but was missing from the schema.
- Corrected the configuration schema's default for `im2deep`'s `reference_dataset`, which
  incorrectly claimed `"Meier_unimod.parquet"`. The actual default falls through to IM2Deep's own
  bundled reference dataset.

## [4.0.0] - 2026-07-17

### Added

- New rescoring engine: **ristretto**, a lean, dependency-light (numpy/scikit-learn/pandas only)
  reimplementation of the Percolator/Käll semi-supervised algorithm, purpose-built for
  MS²Rescore.
- New `rescoring` configuration option: `train_fdr` and `model` (`"svm"`, default, or `"lda"`,
  faster but less powerful).
- New top-level `report_fdr` configuration option: FDR threshold used for console-logged
  identification counts, the HTML report's stats/charts, and FlashLFQ output filtering.
  Previously hardcoded at 1% throughout.
- Automatic inference of search-engine score direction (higher-is-better vs. lower-is-better)
  via spectrum-competed target-decoy evaluation, replacing the user-set `lower_score_is_better`
  option. Grouped by run, so multi-file input sharing native spectrum/scan IDs across runs
  doesn't corrupt the inferred direction.
- Rescoring result tables (`<prefix>.psms.tsv`, `.peptidoforms.tsv`, `.peptides.tsv`,
  `.proteins.tsv`, `.weights.tsv`) are now always written as plain TSV.
- `ms2rescore-report` CLI: new `--fdr` option to regenerate a report at a different FDR
  threshold without rerunning rescoring.
- New MS2 feature generator using Rust-based `ms2rescore_rs` for direct spectrum feature
  extraction (intensity ratios, matched ion counts/percentages, hyperscore).
- New basic features: `theoretical_mass`, `experimental_mass`, `mass_error`, `pep_len`.
- `annotate_spectra()`: annotates all PSM spectra once before feature generators run,
  eliminating redundant per-generator spectrum parsing.
- Top-level configuration options `fragmentation_model`, `tolerance_value`, and `tolerance_mode`
  to control centralized fragment ion annotation (defaults: `cidhcd`, `0.02 Da`).
- **Mumble** integration (optional, beta): a new PSM generator for exploring alternative peptide
  identifications via candidate mass-shift modifications (`pip install ms2rescore[mumble]`).
- Intermediate file output on feature-generation or rescoring errors, enabling recovery by
  rerunning with a modified configuration instead of restarting from scratch.
- Feature generators are intelligently skipped when all their features are already present in
  the input PSM file (e.g., on a recovery run).
- Standalone HTML report regeneration from a PSM TSV file alone -- no config or log file
  required; before/after comparisons are reconstructed from the PSM list's provenance data.
- `ParseSpectrumError` exception for spectrum-parsing failures.

### Changed

- MS2 and MS2PIP feature calculation migrated to Rust via `ms2rescore_rs` (~5x speed-up).
- Spectrum files are parsed and annotated once, up front, and shared across all feature
  generators.
- DeepLC upgraded to its v4 API: dataset-wide processing with per-run calibration or finetuning.
  New multitask model leads to much improved performance, even without finetuning.
- IM2Deep upgraded to its v2 API (`im2deep>=2.0.1`): dataset-wide processing with per-run CCS
  calibration using reference peptides.
- Basic feature generator uses fixed charge encoding (charges 1-6) instead of a dynamic
  per-dataset range.
- HTML report generation (in-run and standalone) reconstructs before/after rescoring
  comparisons from the main PSM list's provenance data, rather than relying on separately
  persisted result tables.
- Report/identification-overlap comparisons key on `(run, spectrum_id)` instead of bare
  `spectrum_id`, so multiple input files reusing the same native spectrum IDs no longer collide.
- Multi-run PSM lists are disambiguated during rescoring/competition via a run identifier,
  instead of relying on spectrum ID alone.
- `max_psm_rank_output > 1` now applies consistently across the main output, rescoring tables,
  and report: multiple ranked PSMs per spectrum, with q-values/PEPs computed per-row rather than
  through full spectrum competition. Intended for surfacing ambiguous candidates (e.g. from
  Mumble), not a statistically rigorous FDR-controlled count.
- Main PSM list output renamed `<prefix>.psms.tsv` → `<prefix>.tsv`; the crash-recovery
  intermediate file renamed the same way (`<prefix>.intermediate.tsv`).
- Dependencies upgraded: `deeplc>=4.0.0`, `im2deep>=2.0.1`, `ms2pip>=4.2.0`,
  `ms2rescore_rs>=0.5.0`. Added `pyarrow`.
- numpy 2.0 compatibility.
- Python 3.11 or newer is now required.

### Removed

- [BREAKING] **Mokapot** rescoring engine and dependency removed.
- [BREAKING] **Percolator CLI integration** removed (the separate engine that shelled out to a
  locally-installed `percolator` binary). ristretto is now the only rescoring engine.
- [BREAKING] MaxQuant feature generator removed; functionality consolidated into the MS2
  feature generator.
- [BREAKING] ionmob feature generator removed; replaced by IM2Deep v2.
- [BREAKING] `rescoring_engine` configuration option removed (mokapot/Percolator-specific:
  `fasta_file`, `write_weights`, `write_txt`, `protein_kwargs`), replaced by `rescoring`.
- [BREAKING] Top-level `fasta_file` configuration option and FASTA-based protein inference
  removed (mokapot-specific). Does not affect Mumble's separate
  `psm_generator.mumble.fasta_file` option, which is unrelated and unchanged.
- [BREAKING] `lower_score_is_better` configuration option removed; score direction is always
  auto-inferred now, with no configuration-level override.
- [BREAKING] `write_rescoring_tables` configuration option removed -- rescoring tables are
  unconditionally written.
- [BREAKING] PIN (Percolator) file output removed; the main PSM list TSV already carries all
  rescoring features.
- [BREAKING] Ability to skip rescoring via configuration removed; rescoring always runs.
- [BREAKING] `ms2_tolerance`, `spectrum_path`, and `spectrum_id_pattern` parameters removed from
  `MS2PIPFeatureGenerator`. Fragment mass tolerance is set globally via `tolerance_value` /
  `tolerance_mode`.
- [BREAKING] `spectrum_path`, `spectrum_id_pattern`, `mass_mode`, and `processes` parameters
  removed from `MS2FeatureGenerator`. Spectra are provided via centralized `annotate_spectra()`.
- [BREAKING] `ms2rescore.utils` (public Python API) renamed and split into two internal modules,
  `ms2rescore._utils` and `ms2rescore._ristretto_utils` -- neither is part of the public API.
- `deeplcretrainer` dependency removed (functionality merged into DeepLC v4).
- `tomli` dependency removed (only required for Python <3.11).

### Fixed

- `processes=-1` (ms2rescore default) passed to DeepLC `num_threads`, which requires a positive
  integer or `None`.
- Q-value NaN check in `parse_psms.py` failed when `qvalue` array contained `None` values.
- `BrokenExecutor` not caught during mokapot rescoring, producing unclear crashes on worker
  failure.
- GUI runs never wrote an HTML log file, unlike CLI runs.
- Out-of-memory errors from multiprocessing during spectrum parsing.
