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

## [3.2.1] - 2026-02-10

### Changed

- Updated installer dependency pins (`pyopenms` 3.5, `psm_utils` 1.5, `ms2rescore-rs` 0.4.3,
  `deeplc` 3.1.13, `im2deep` 1.2.0).

### Removed

- Dropped support for Python 3.10 (required for DeepLC).

### Fixed

- `report`: caught `IndexError` when no confidently identified PSMs are present.
- `report`: fixed empty charts caused by a mismatching Plotly.js version; the version is now
  determined dynamically.
- GUI: fixed missing log file when running MS²Rescore through the GUI.
- GUI: fixed `UnicodeEncodeError`/`None` stdout-stderr issue when running DeepLC through the GUI
  in debug mode.
- PyInstaller: fixed missing dependencies in the Windows installer.
- Fixed a typo in the IM2Deep docstring's default argument.

## [3.2.0] - 2025-10-20

### Added

- Update notifications for both GUI and CLI, with a new `disable_update_check` configuration
  option.
- IM2Deep added to the citations interface.
- If observed m/z values are present in the PSM list, always verify that PSMs match the
  corresponding spectra.
- Validation for `psm_id_pattern` and spectrum ID regex patterns, with clear error messages when
  patterns are invalid or don't match expected spectrum IDs.

### Changed

- FlashLFQ export switched to `psm_utils`'s writer; removed
  `ms2rescore > mokapot > write_flashlfq`, replaced by a new top-level `write_flashlfq` option.
- Rescoring no longer continues if precursor information (m/z, RT, or IM) is missing for a subset
  of PSMs; previously such PSMs were silently removed.
- Improved handling of missing LC-IM-MS fields in spectrum files.
- Clearer logging while parsing precursor information; always raises an error if a PSM cannot be
  found in the spectrum files, including example PSM IDs from both PSM and spectrum files.
- Boolean configuration options now permit `null`/`None` to explicitly represent an unset state.
- Dropped support for Python 3.9 (EOL); test matrix extended to 3.12-3.14 (verified through 3.12).
  Default runtime set to Python 3.12 for Docker and the Windows installer.
- Dependencies upgraded: `deeplc>=3.1`, `im2deep>=0.3.1`.
- Switched to `uv` for package management, with a lockfile for reproducible installs; the version
  single source of truth now lives in `pyproject.toml`.
- Docker base image upgraded to `python:3.12-slim`, pinned via the lockfile.
- Suppressed non-actionable OpenMS data warnings to reduce noise.

### Fixed

- CLI/config precedence: command-line defaults (e.g. `False`) no longer override `True` values
  set in the configuration file.
- Configuration file's `profile` parameter is now respected; previously only the CLI argument was
  honored.
- Corrected parsing of spectrum files with missing precursor information; previously, precursor
  data was parsed even when not required by the configured feature generators.
- Avoided incompatibility between `CUDA_VISIBLE_DEVICES` and XGBoost causing a silent crash on
  Windows when set to `-1`.
- Replaced `Inf`/`np.Inf` usage with `inf` for NumPy 2.0 compatibility.

## [3.1.5] - 2025-04-15

### Fixed

- Fixed compatibility with NumPy 2: usage of `np.inf` instead of `np.Inf`.

## [3.1.4] - 2024-12-04

### Fixed

- Updated IM2Deep import after an API change in v0.3.1.

## [3.1.3] - 2024-10-02

### Fixed

- Pinned DeepLC version to `<3.1`, avoiding a calibration bug.
- Pinned pyOpenMS (upstream dependency of `psm_utils`).
- Fixed incorrect decoy pattern configuration in a documentation example.
- Fixed `UnicodeEncodeError` when running IM2Deep.

## [3.1.2] - 2024-09-18

### Changed

- GUI: improved user experience with descriptions, opening the report upon finishing, etc.
- Updated dependency versions; dropped support for Python 3.8 (required for TensorFlow, nearing
  EOL).

### Fixed

- Fixed ruff linting for the tutorial notebook.
- Fixed a bug where the default TIMS²Rescore configuration items always overwrote user
  configuration.
- Fixed a non-descriptive `IndexError` in Qvality when PEP cannot be calculated.
- Fixed `UnicodeEncodeError` when running DeepLC with transfer learning.

## [3.1.1] - 2024-08-14

### Fixed

- CI: added missing `ionmob` dependency for the Windows installer.
- GUI: correctly parse IM2Deep configuration.
- Limited Percolator processes to 128.

### Changed

- Docker images are now built and pushed to GHCR.

## [3.1.0] - 2024-07-21

### Added

- New `tims2rescore` CLI entrypoint with timsTOF-specific configuration defaults.
- IM2Deep feature generator for ion mobility.
- Support for Bruker raw (`.d`) and miniTDF formats (back-ported to v3.0.3).
- `psm_id_rt_pattern` and `psm_id_im_pattern` options to parse RT/IM from PSM file spectrum
  identifiers.
- `max_psm_rank_input` and `max_psm_rank_output` options to control rescoring and reporting of
  multiple PSMs per spectrum (default 5 and 1).
- Peptide-level q-values and PEPs are now also stored in the output PSM list, under the
  `metadata` field.
- `profile` option to write a cProfile report of the rescoring process.
- Option to pass `train_fdr` to Mokapot (`test_fdr` was already supported).
- Faster spectrum file reading with the `timsrust` and `mzdata` Rust libraries.

### Changed

- CI: switched linting to ruff (from flake8).
- DeepLC: by default use PSMs passing a 1% FDR threshold for calibration, instead of a fixed
  number.
- DeepLC: warn when using `deeplc_retrain` with fewer than 500 calibration PSMs.
- Raise an exception if precursor info (retention time or ion mobility) is missing from both PSM
  and spectrum files.
- Dockerfile: use a Python 3.11 base image instead of Ubuntu.
- Config: changed `write_report` default from `true` to `false`.

### Fixed

- Report: fixed calculation of the number of retained peptides in the identification overlap
  chart.
- DeepLC: fixed an issue with conversion to PEPREC.
- Fixed a bug in `infer_spectrum_path` when the configured path is `None`.
- Fixed a bug that always set `im_required` to `True`, even without IM-based feature generators
  configured.

## [3.0.3] - 2024-04-08

### Added

- Support for passing `model_dir` to MS²PIP in the feature generator.

### Changed

- Missing precursor information (retention time and ion mobility) is now parsed in Rust using
  `mzdata` and `timsrust`, via the new `ms2rescore-rs` package. Adds support for Bruker TDF and
  miniTDF raw files for this feature.

### Fixed

- Fixed writing of the FlashLFQ output file with Mokapot when PSM list `run` fields are empty.
- Write a partial report when confidence before rescoring could not be calculated (e.g. no PSMs
  identified before rescoring).

## [3.0.2] - 2024-03-15

### Fixed

- Removed PSMs with invalid amino acids to avoid errors in feature generators.

## [3.0.1] - 2024-02-22

### Changed

- Improved exception handling and logging for report generation.
- Updated `deeplc_retrainer` version.
- Pinned NumPy version for Python 3.11 (incompatibility between pyGAM, scikit-learn, and
  TensorFlow).

### Fixed

- No longer use an unsupported icon on Linux systems.

## [3.0.0] - 2023-11-30

### Added

- Adopted `psm_utils` to support over 10 PSM file formats from different search engines.
- Support for mzML (MGF still supported).
- HTML output with interactive charts for quality control and result inspection.
- `ionmob` feature generator for ion mobility CCS predictions.
- `mokapot` rescoring engine for more efficient rescoring.
- New graphical user interface with a one-click installer for Windows.
- Documentation on ReadTheDocs.

### Changed

- Major refactor into a modular Python package to facilitate use in other workflows (see
  preprint: https://doi.org/10.26434/chemrxiv-2023-rvr9n).

## [2.1.3] - 2022-08-30

### Fixed

- GUI: fixed parsing of new-style MaxQuant modification names.

## [2.1.2] - 2022-04-19

### Added

- Additional PEAKS search engine features.

### Fixed

- Fixed the MGF TITLE field extraction to be compatible with PEAKS XPro output.
- Removed `shell=True` from the Percolator subprocess call.
- MS²Rescore now returns an error if Percolator fails (`check=True` added to `subprocess.run`).
- Fixed a `DeprecationWarning` from passing a set as a pandas indexer.
- Python 3.9 compatibility for the GUI.

### Changed

- Ignore Deprecation/Future/UserWarnings in GUI mode.
- Removed the `importlib_resources` backport (Python 3.6 compatibility no longer needed).
- GUI can now be started from the terminal with `ms2rescore-gui`, in addition to
  `python -m ms2rescore.gui`.
- Relaxed the scikit-learn dependency to `<2`.

## [2.1.1] - 2022-03-04

### Added

- Refactored parsing of the MaxQuant `Modified sequence` field to support new long-form
  modification labels (e.g. `PEPM(Oxidation (M))K` instead of `PEPM(ox)K`).
- Added a regex pattern option to the config for more flexible PSM ID matching against the MGF
  `TITLE` field.
- Added `MGF TITLE pattern`, `num CPU`, and `feature set` options to the GUI.

### Changed

- Replaced deprecated pandas functions.

### Fixed

- Spaces and hyphens are now supported in input file names (changed Percolator subprocess
  command).
- Default value of `mgf_title_pattern` set to `None` for compatibility with the PEAKS pipeline.

## [2.1.0] - 2021-12-14

### Added

- Support for PEAKS mzIdentML identification files.

### Changed

- Retention time predictor now calibrates per raw file independently, resulting in more accurate
  calibrations.
- PIN pipeline: allow mass-shift modification labels with a `+` sign (e.g. `R.IM[+15.99492]MAR.D`).
- PEAKS and MaxQuant pipelines: precursor charge from the ID file now overwrites `CHARGE` in the
  MGF file.
- Plotting module can only run when Percolator is run as part of MS²Rescore.
- MS²Rescore now raises an error and stops if no decoy PSMs are present.
- Removed the obsolete `mgf` option from `maxquant_to_rescore` configuration.

### Fixed

- MaxQuant pipeline: carbamidomethyl is no longer a default fixed modification (it could not
  previously be overridden); it must now be specified explicitly if applicable.
- Fixed MGF parsing incorrectly removing all lines ending in `0.0` (including lines ending in
  `10.0`) instead of only zero-intensity peaks.

## [2.0.0] - 2021-11-21

### Added

- MS²Rescore now runs on Windows.
- Graphical user interface (GUI), built with Gooey (`pip install ms2rescore[gui]`,
  `python -m ms2rescore.gui`).
- Install, upgrade, and run scripts for Windows users without an existing Python installation.
- MS²Rescore can now output a PDF with result plots.

### Changed

- [BREAKING] Renamed the package from "MS²ReScore" to "ms2rescore" (lowercase s); some class
  names and exceptions changed accordingly, breaking the Python API with previous versions.
- [BREAKING] `feature_sets` option now expects a list of strings (`searchengine`, `rt`,
  `ms2pip`) instead of a predefined feature-set string such as `ms2pip_rt`. By default,
  MS²Rescore uses all three.
- MS²PIP is now called through its Python API instead of the CLI.
