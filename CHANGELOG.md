# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [4.0.0] - 2026-07-17

### Added

- `annotate_spectra()` in `parse_spectra.py`: annotates all PSM spectra once before feature
  generators run, eliminating redundant per-generator spectrum annotation.
- Top-level configuration options `fragmentation_model`, `tolerance_value`, and `tolerance_mode`
  to control centralized fragment ion annotation. Defaults: `cidhcd`, `0.02 Da`.
- New `rescoring` configuration option for the ristretto rescoring engine: `train_fdr` and
  `model` (`"svm"`, default, or `"lda"`, faster but less powerful). Accepts `{}` for defaults,
  or a partial dict -- missing keys fall back to ristretto's own defaults.
- New top-level `report_fdr` configuration option: FDR threshold used for console-logged
  identification counts, the HTML report's stats/charts, and FlashLFQ output filtering.
  Previously hardcoded at 1% throughout.
- `ms2rescore-report` CLI: new `--fdr` option to regenerate a report at a different FDR
  threshold without rerunning rescoring.
- Automatic inference of search-engine score direction (higher-is-better vs. lower-is-better)
  via spectrum-competed target-decoy evaluation, replacing the user-set `lower_score_is_better`
  option. Grouped by run, so multi-file input sharing native spectrum/scan IDs across runs
  doesn't corrupt the inferred direction.
- Rescoring result tables (`<prefix>.psms.tsv`, `.peptidoforms.tsv`, `.peptides.tsv`,
  `.proteins.tsv`, `.weights.tsv`) are now always written as plain TSV, independent of rescoring
  engine internals.
- `ristretto-ms` dependency: a lean, dependency-light (numpy/scikit-learn/pandas)
  reimplementation of the Percolator/Käll semi-supervised rescoring algorithm.
- New `ms2` feature generator backed by `ms2rescore-rs`, providing direct-spectrum features
  including matched-ion counts and percentages, hyperscore, and intensity-ratio features for
  the a, b, c, x, y, and z ion series.
- Configurable IM2Deep CCS reference datasets. CSV, compressed CSV, and Parquet files are
  supported and must contain `peptidoform` and `CCS` columns.
- Mumble integration as an optional dependency (`ms2rescore[mumble]`) for generating
  mass-shift candidate peptidoforms before rescoring.
- Intermediate-file recovery: the current PSM state is written to
  `<prefix>.intermediate.tsv` when feature generation or rescoring fails.
- Feature generators are skipped when all of their output features are already present,
  allowing interrupted runs to resume without recomputing completed features.
- Full per-feature-generator test coverage for the basic, DeepLC, IM2Deep, MS2, and MS²PIP
  feature generators, together with end-to-end rescoring integration tests.
- GUI: rescoring model selector (svm/lda); centralized fragmentation-model and fragment-tolerance
  controls; GUI runs now also write an HTML log file (`<prefix>.log.html`), matching the CLI.

### Changed

- Spectrum annotation is now performed once in `core.py` before all feature generators run.
  MS²PIP and MS2 feature generators reuse `AnnotatedMS2Spectrum` objects attached to each PSM.
- MS²PIP: migrated from `correlate_preloaded` back to the now-unified `correlate()` API. Spectra
  are passed via `psm.spectrum`.
- MS2: migrated from `ms2_features_from_ms2spectra` to `score_ms2_spectra` API. Feature set
  expanded to cover all ion series (a, b, c, x, y, z).
- DeepLC migrated to the v4 functional API, with dataset-wide prediction and selection of the
  best multitask prediction head independently for each run.
- IM2Deep migrated to the v2 functional API, with dataset-wide prediction and per-run linear
  CCS calibration.
- Dependencies bumped: `deeplc>=4.0.0b1`, `im2deep>=2.0.1`, `ms2pip>=4.2.0`,
  `ms2rescore_rs>=0.5.0`, and `ristretto-ms>=0.3.0`. Added `pyarrow>=14`.
- numpy 2.0 compatibility in `charts.py` (`np.trapz` → `np.trapezoid`).
- Rescoring engine replaced: mokapot → ristretto. Rescoring can no longer be skipped -- it
  always runs.
- Main PSM list output renamed `<prefix>.psms.tsv` → `<prefix>.tsv`; the crash-recovery
  intermediate file renamed the same way (`<prefix>.intermediate.tsv`).
- HTML report generation (both in-run and standalone via `ms2rescore-report`) reconstructs
  before/after rescoring comparisons from the main PSM list's provenance data and current state,
  rather than relying on separately persisted result tables.
- Report/identification-overlap comparisons key on `(run, spectrum_id)` instead of bare
  `spectrum_id`, so multiple input files reusing the same native spectrum IDs no longer collide.
- Multi-run PSM lists are disambiguated during rescoring/competition via ristretto's `run_col`,
  instead of relying on `spectrum_id` alone.
- `max_psm_rank_output > 1` now applies consistently across the main output, rescoring tables,
  and report: multiple ranked PSMs per spectrum, with q-values/PEPs computed per-row rather than
  through full spectrum competition. Intended for surfacing ambiguous candidates (e.g. from
  Mumble), not a statistically rigorous FDR-controlled count.
- Protein-level rollups use ristretto's picked-protein competition (Savitski et al. 2015) when
  `id_decoy_pattern` is set.
- Python 3.11 or newer is now required.

### Removed

- [BREAKING] `ms2_tolerance`, `spectrum_path`, and `spectrum_id_pattern` parameters removed from
  `MS2PIPFeatureGenerator`. Fragment mass tolerance is set globally via `tolerance_value` /
  `tolerance_mode` in the top-level configuration (default: `0.02 Da`).
- [BREAKING] `spectrum_path`, `spectrum_id_pattern`, `mass_mode`, and `processes` parameters
  removed from `MS2FeatureGenerator`. Spectra are provided via centralized `annotate_spectra()`.
- [BREAKING] `ionmob` feature generator removed. Ion-mobility prediction and CCS calibration are
  now provided by the IM2Deep v2 feature generator.
- [BREAKING] `maxquant` feature generator removed. Its relevant spectrum-derived features are
  now provided by the `ms2` feature generator.
- [BREAKING] Mokapot rescoring engine and the `mokapot` dependency removed, along with the
  `ms2rescore.rescoring_engines` module.
- [BREAKING] `rescoring_engine` configuration option removed (mokapot-specific: `fasta_file`,
  `write_weights`, `write_txt`, `protein_kwargs`), replaced by `rescoring` (see Added).
- [BREAKING] `fasta_file` configuration option and FASTA-based protein inference removed.
- [BREAKING] `lower_score_is_better` configuration option removed. Score direction is now always
  auto-inferred (see Added) with no config-level override.
- [BREAKING] `write_rescoring_tables` configuration option removed -- rescoring tables are
  unconditionally written now.
- [BREAKING] PIN (Percolator) file output removed for `log_level=debug` -- the main PSM list
  TSV already carries all rescoring features.
- [BREAKING] Ability to skip rescoring via configuration removed. `rescoring: null` is now
  rejected by config validation instead of being silently ignored; rescoring always runs.
- [BREAKING] The public `ms2rescore.utils` module removed. General internal utilities now reside
  in `ms2rescore._utils`, while private ristretto integration helpers reside in
  `ms2rescore._ristretto_utils`; neither module is part of the public API surface.

### Fixed

- MS²PIP features incorrectly computed for multi-rank PSMs (`max_psm_rank_input > 1`): all
  PSMs sharing a spectrum ID received the annotation of the first-seen PSM, producing a bimodal
  `spec_pearson_norm` distribution. Fixed in ms2pip (per-PSM annotation) and reflected in
  ms2rescore via centralized per-PSM `annotate_spectra()`.
- DeepLC RT features incorrectly assigned across PSMs: missing `sort_index()` after q-value sort
  for calibration caused PSMs to receive another PSM's RT predictions.
- Multi-run MS²PIP crash caused by observed-spectrum deduplication using `spectrum_id` alone.
  Spectra are now keyed by `(run, spectrum_id)`, preventing collisions between runs that reuse
  native scan identifiers.
- Missing precursor m/z or search-engine scores produced `NaN` basic features instead of the
  documented `0` values because float `NaN` values were not detected by the existing `None`
  checks.
- DeepLC multitask models always used prediction head 0 instead of selecting the head that best
  matched each run.
- Mumble-generated candidate PSMs leaked into DeepLC and IM2Deep calibration-set selection
  because they inherited the original PSM's score and q-value. Generated candidates are now
  excluded using `get_original_hit_mask()`.
- Configured IM2Deep reference datasets were incorrectly reported as missing because the
  file-existence check was inverted.
- `processes=-1` (ms2rescore default) passed to DeepLC `num_threads`, which requires a positive
  integer or `None`.
- Q-value NaN check in `parse_psms.py` failed when `qvalue` array contained `None` values.
- Fragment mass tolerance fallback defaults in `core.py` incorrectly set to `20.0 ppm` instead
  of `0.02 Da`.
- Empty or partial `rescoring` configuration dictionaries could raise `KeyError` instead of
  falling back to ristretto defaults.
- GUI runs never wrote an HTML log file (`<prefix>.log.html`), unlike CLI runs -- the GUI's
  logging setup only ever attached a plain text-file handler.


## [3.3.0a1] - 2026-04-09

### Added

- New MS2 feature generator using Rust-based `ms2rescore_rs` for direct spectrum feature
  extraction (intensity ratios, matched ion counts/percentages, hyperscore).
- Mumble integration as an optional PSM generator for exploring alternative peptide
  identifications with mass shift modifications (`pip install ms2rescore[mumble]`).
- Intermediate file output (`.intermediate.psms.tsv`) on feature generation or rescoring errors,
  enabling recovery by rerunning with modified configuration.
- Intelligent skipping of feature generators when all their features are already present in the
  PSM file (e.g., from an intermediate recovery run).
- New basic features: `theoretical_mass`, `experimental_mass`, `mass_error`, `pep_len`.
- Standalone report generation from PSM TSV files without requiring full config or log files.
- `ParseSpectrumError` exception for spectrum parsing failures.

### Changed

- Migrated MS2 and MS2PIP feature calculations to Rust via `ms2rescore_rs`, significantly
  improving performance (~5x speed-up).
- Spectrum files are now parsed once and stored as `MS2Spectrum` objects, replacing the previous
  per-feature-generator parsing approach.
- DeepLC integration upgraded to v4 API: dataset-wide processing, fine-tuning enabled by
  default, `SplineTransformerCalibration` for retention time calibration.
- IM2Deep integration upgraded to v2 API: dataset-wide processing with per-run
  `LinearCCSCalibration` using reference peptides.
- MS2PIP integration upgraded to use preloaded spectra and Rust-based feature calculation.
- Basic feature generator now uses fixed charge encoding (charges 1-6) instead of dynamic
  min-max range.
- Report generation CLI now accepts PSM file path with optional `--output` flag.
- Charge-stripping regex pattern consolidated into shared `CHARGE_PATTERN` constant.
- Upgraded dependencies: `deeplc>=4.0.0a2`, `im2deep>=2.0.0a2`, `ms2pip>=4.2.0a0`,
  `ms2rescore_rs>=0.5.0a0`.

### Removed

- MaxQuant feature generator (functionality consolidated into MS2 feature generator).
- ionmob feature generator (replaced by IM2Deep v2).
- `deeplcretrainer` dependency (functionality merged into DeepLC v4).

### Fixed

- Percolator kwargs silently ignored due to parameter name shadowing local variable.
- Unreachable and broken error handlers in Percolator subprocess execution.
- `fdr` parameter ignored in `_log_id_psms_before` (hardcoded to 0.01).
- Out-of-memory errors from multiprocessing in spectrum parsing.
