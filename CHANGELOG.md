# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed

- MS2PIP integration now strips pre-annotated spectra back to raw `MS2Spectrum` before calling
  `ms2pip.correlate()`, ensuring each PSM is annotated with its own peptidoform.

### Fixed

- MS²PIP features incorrectly computed for multi-rank PSMs (`max_psm_rank_input > 1`): all
  PSMs sharing a spectrum ID received the annotation of the first-seen PSM, producing a bimodal
  `spec_pearson_norm` distribution.
- DeepLC RT features incorrectly assigned across PSMs due to missing `sort_index()` after
  sorting by q-value for calibration, causing every PSM to receive another PSM's RT features.

### Breaking changes

- The `ms2_tolerance` parameter of `MS2PIPFeatureGenerator` has been removed. Fragment mass
  tolerance is now set globally via `tolerance_value` / `tolerance_mode` in the top-level
  configuration. The default is `0.02 Da`, matching the previous MS²PIP default.

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
