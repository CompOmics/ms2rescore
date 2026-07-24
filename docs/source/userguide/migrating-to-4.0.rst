.. _migrating-to-4.0:

Migrating to 4.0
=================

MS²Rescore 4.0 is a major release: the rescoring engine, feature generation, and reporting
were all substantially reworked. This guide covers what changed from a user's perspective and
what to update in existing configuration files, scripts, and pipelines when upgrading from
3.2.x.

.. admonition:: Summary
   :class: note

   Rescoring now runs on **ristretto** instead of mokapot or Percolator -- both of those
   engines are gone, along with their configuration options. Feature generation is faster
   (Rust-based MS2/MS2PIP calculation, centralized spectrum parsing) and DeepLC/IM2Deep were
   upgraded to their latest major versions. The main PSM output file and rescoring result
   tables were renamed. See below for the concrete changes to make.

Computational performance
--------------------------

Internal benchmarking across several representative datasets shows that 4.0 uses substantially
less CPU time than 3.2.x, mainly thanks to the reworked feature generation step and the faster
rescoring engine. Peak memory usage is higher in 4.0, and scales with the size of the raw spectrum
files. Spectra are now read once into memory instead of separately per worker process, so this is
worth taking into account on memory-constrained machines or with large acquisitions. Identification
rates are unchanged between versions.

Rescoring engine: mokapot/Percolator → ristretto
-------------------------------------------------

MS²Rescore no longer supports mokapot or a standalone Percolator installation as rescoring
engines. Both are removed. Rescoring now always runs through **ristretto**, a lean,
dependency-light reimplementation of the same underlying semi-supervised algorithm.

If your configuration has a ``rescoring_engine`` section, replace it with ``rescoring``:

.. tab:: Before (3.2.x)

  .. code-block:: json

    {
      "ms2rescore": {
        "rescoring_engine": {
          "mokapot": {
            "train_fdr": 0.01,
            "fasta_file": "proteins.fasta",
            "write_weights": true,
            "write_txt": true
          }
        }
      }
    }

.. tab:: After (4.0.0)

  .. code-block:: json

    {
      "ms2rescore": {
        "rescoring": {
          "train_fdr": 0.01,
          "model": "svm"
        }
      }
    }

Notes on that change:

- ``model`` chooses between ``"svm"`` (default, iterative Percolator-style SVM) and ``"lda"``
  (a faster but less powerful single-pass Fisher LDA).
- The top-level ``fasta_file`` option and FASTA-based protein inference are removed entirely --
  there is no direct replacement. Protein-level rollups still work; they use picked-protein
  competition when ``id_decoy_pattern`` is set. (This does not affect Mumble's own, unrelated
  ``psm_generator.mumble.fasta_file`` option, which is still used to validate amino-acid-
  combination candidates.)
- ``write_weights``/``write_txt`` are gone: rescoring result tables (PSM, peptidoform, peptide,
  protein, and feature-weight tables) are now **always** written as plain TSV. There is no
  config option to disable this.
- Rescoring can no longer be skipped. If you previously set ``rescoring_engine`` to an empty
  object (or ``rescoring: null``) to only write engineered features, that path is gone --
  rescoring always runs.
- If you had a locally installed ``percolator`` binary configured via the old Percolator CLI
  integration, that integration no longer exists. Switch to the default ``rescoring`` config
  above.

Score direction is now automatic
---------------------------------

The ``lower_score_is_better`` option is removed. MS²Rescore now automatically infers whether a
higher or lower search-engine score is better, by running a small target-decoy competition
under both hypotheses and picking whichever identifies more PSMs at the configured
``rescoring.train_fdr``. This requires no configuration and handles multi-file input correctly
(spectrum IDs are grouped by run so colliding scan numbers across files don't affect the
result).

If your search engine's score direction is genuinely ambiguous (very few decoys, heavily
quantized scores), double-check the inferred direction in the log output at ``debug`` level.

New: ``report_fdr``
--------------------

The FDR threshold used for console-logged identification counts, the HTML report's charts, and
FlashLFQ output filtering was previously hardcoded to 1%. It's now a configurable top-level
option:

.. code-block:: json

  {
    "ms2rescore": {
      "report_fdr": 0.01
    }
  }

This is independent of ``rescoring.train_fdr``, which only controls ristretto's training.

Output file names changed
--------------------------

+-------------------------------------+-----------------------------------------------------------------------------------------------------------------------+
| 3.2.x                               | 4.0.0                                                                                                                 |
+=====================================+=======================================================================================================================+
| ``<prefix>.psms.tsv`` (main output) | ``<prefix>.tsv``                                                                                                      |
+-------------------------------------+-----------------------------------------------------------------------------------------------------------------------+
| ``<prefix>.intermediate.psms.tsv``  | ``<prefix>.intermediate.tsv``                                                                                         |
+-------------------------------------+-----------------------------------------------------------------------------------------------------------------------+
| *(mokapot weights/txt output)*      | ``<prefix>.psms.tsv``, ``.peptidoforms.tsv``, ``.peptides.tsv``, ``.proteins.tsv``, ``.weights.tsv`` (always written) |
+-------------------------------------+-----------------------------------------------------------------------------------------------------------------------+
| ``<prefix>.psms.pin``               | removed (no PIN output)                                                                                               |
+-------------------------------------+-----------------------------------------------------------------------------------------------------------------------+

If you have downstream scripts or pipelines that read ``<prefix>.psms.tsv`` expecting the main
PSM list, update them to read ``<prefix>.tsv`` instead -- ``<prefix>.psms.tsv`` is now the
PSM-level *rescoring result table* (a leaner view: identifiers and score/q-value/PEP only, no
rescoring features or provenance data).

Feature generators: MaxQuant and ionmob removed
------------------------------------------------

- The MaxQuant feature generator is removed; its functionality is folded into the MS2 feature
  generator, which now works directly from parsed spectra regardless of search engine.
- The ionmob feature generator is removed, replaced by the upgraded IM2Deep v2 integration.
- ``MS2PIPFeatureGenerator`` no longer accepts ``ms2_tolerance``, ``spectrum_path``, or
  ``spectrum_id_pattern`` -- fragment mass tolerance and spectrum access are now centralized at
  the top level via ``tolerance_value``, ``tolerance_mode``, and ``spectrum_path``/
  ``spectrum_id_pattern`` on the main configuration.
- ``MS2FeatureGenerator`` no longer accepts ``spectrum_path``, ``spectrum_id_pattern``,
  ``mass_mode``, or ``processes`` for the same reason.

If you configured either of these generators directly, remove those per-generator parameters
and set ``tolerance_value``/``tolerance_mode``/``fragmentation_model`` at the top level instead.

Python API
----------

``ms2rescore.utils`` (previously public) is renamed and split into two internal modules,
``ms2rescore._utils`` and ``ms2rescore._ristretto_utils``, neither of which is part of the
public API -- if any code imports from ``ms2rescore.utils`` directly, that will break. The
documented, public entry points remain ``ms2rescore.rescore()`` for the full pipeline and
``ms2rescore.rescoring.rescore()`` for calling ristretto directly on a feature DataFrame; see the
:doc:`Python API tutorial </tutorials/in-depth-python-api>` for the updated, step-by-step
version.

MS²Rescore 4.0.0 requires **Python 3.11 or newer**.

New, optional
-------------

Nothing below requires any migration, but may be worth adopting:

- **Mumble** (beta): an optional PSM generator that proposes alternative peptide identifications
  for a PSM's precursor mass shift. Install with ``pip install ms2rescore[mumble]`` and see the
  :doc:`Mumble guide </userguide/mumble>`.
- **GUI**: runs now also produce an HTML log file (``<prefix>.log.html``), matching the CLI.
- **Standalone report regeneration**: ``ms2rescore-report`` can rebuild the HTML report from just
  the main PSM TSV file, with an optional ``--fdr`` override, without needing the original config
  or log files.
