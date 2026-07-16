# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- `annotate_spectra()` in `parse_spectra.py`: annotates all PSM spectra once before feature
  generators run, eliminating redundant per-generator spectrum annotation.
- Top-level configuration options `fragmentation_model`, `tolerance_value`, and `tolerance_mode`
  to control centralized fragment ion annotation. Defaults: `cidhcd`, `0.02 Da`.
- New `rescoring` configuration option for the ristretto rescoring engine: `train_fdr` and
  `model` (`"svm"`, default, or `"lda"`, faster but less powerful). Accepts `null` to skip
  rescoring, `{}` for defaults, or a partial dict -- missing keys fall back to ristretto's own
  defaults.
- New top-level `report_fdr` configuration option: FDR threshold used for console-logged
  identification counts, the HTML report's stats/charts, and FlashLFQ output filtering.
  Previously hardcoded at 1% throughout.
- `ms2rescore-report` CLI: new `--fdr` option to regenerate a report at a different FDR
  threshold without rerunning rescoring.
- Automatic inference of search-engine score direction (higher-is-better vs. lower-is-better)
  via spectrum-competed target-decoy evaluation, replacing the user-set `lower_score_is_better`
  option.
- Rescoring result tables (`<prefix>.psms.tsv`, `.peptidoforms.tsv`, `.peptides.tsv`,
  `.proteins.tsv`, `.weights.tsv`) are now always written as plain TSV, independent of rescoring
  engine internals.
- `ristretto-ms` dependency: a lean, dependency-light (numpy/scikit-learn/pandas)
  reimplementation of the Percolator/Käll semi-supervised rescoring algorithm.
- GUI: rescoring model selector (svm/lda); GUI runs now also write an HTML log file
  (`<prefix>.log.html`), matching the CLI.

### Changed

- Spectrum annotation is now performed once in `core.py` before all feature generators run.
  MS²PIP and MS2 feature generators reuse `AnnotatedMS2Spectrum` objects attached to each PSM.
- MS²PIP: migrated from `correlate_preloaded` back to the now-unified `correlate()` API. Spectra
  are passed via `psm.spectrum`.
- MS2: migrated from `ms2_features_from_ms2spectra` to `score_ms2_spectra` API. Feature set
  expanded to cover all ion series (a, b, c, x, y, z).
- Dependencies bumped: `ms2pip>=4.2.0b1`, `ms2rescore_rs>=0.5.0b1`. Added `pyarrow>=14`.
- numpy 2.0 compatibility in `charts.py` (`np.trapz` → `np.trapezoid`).
- Rescoring engine replaced: mokapot → ristretto.
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

### Fixed

- MS²PIP features incorrectly computed for multi-rank PSMs (`max_psm_rank_input > 1`): all
  PSMs sharing a spectrum ID received the annotation of the first-seen PSM, producing a bimodal
  `spec_pearson_norm` distribution. Fixed in ms2pip (per-PSM annotation) and reflected in
  ms2rescore via centralized per-PSM `annotate_spectra()`.
- DeepLC RT features incorrectly assigned across PSMs: missing `sort_index()` after q-value sort
  for calibration caused PSMs to receive another PSM's RT predictions.
- `processes=-1` (ms2rescore default) passed to DeepLC `num_threads`, which requires a positive
  integer or `None`.
- Q-value NaN check in `parse_psms.py` failed when `qvalue` array contained `None` values.
- `BrokenExecutor` not caught in mokapot rescoring engine.
- Fragment mass tolerance fallback defaults in `core.py` incorrectly set to `20.0 ppm` instead
  of `0.02 Da`.
- GUI runs never wrote an HTML log file (`<prefix>.log.html`), unlike CLI runs -- the GUI's
  logging setup only ever attached a plain text-file handler.

### Breaking changes

- `ms2_tolerance`, `spectrum_path`, and `spectrum_id_pattern` parameters removed from
  `MS2PIPFeatureGenerator`. Fragment mass tolerance is set globally via `tolerance_value` /
  `tolerance_mode` in the top-level configuration (default: `0.02 Da`).
- `spectrum_path`, `spectrum_id_pattern`, `mass_mode`, and `processes` parameters removed from
  `MS2FeatureGenerator`. Spectra are provided via centralized `annotate_spectra()`.
- Mokapot rescoring engine and the `mokapot` dependency removed, along with the
  `ms2rescore.rescoring_engines` module.
- `rescoring_engine` configuration option removed (mokapot-specific: `fasta_file`,
  `write_weights`, `write_txt`, `protein_kwargs`), replaced by `rescoring` (see Added).
- `fasta_file` configuration option and FASTA-based protein inference removed.
- `lower_score_is_better` configuration option removed (now auto-inferred).
- `write_rescoring_tables` configuration option removed -- rescoring tables are unconditionally
  written now.
- PIN (Percolator) file output removed for disabled rescoring / `log_level=debug` -- the main
  PSM list TSV already carries all rescoring features.
- `ms2rescore.utils` (public) renamed to `ms2rescore._utils` (internal) and merged with the new
  rescoring integration layer -- no longer part of the public API surface.

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
