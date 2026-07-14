.. _mumble:

Open modification searching with Mumble
========================================

.. warning::
  Mumble is under active development and should be considered **beta** software. It is
  usable and already improves identification rates on open-modification-search-style data,
  but occasional errors can still occur, especially on unusual input files or uncommon
  modifications. Please report issues on the
  `Mumble GitHub repository <https://github.com/compomics/mumble>`_.

What Mumble does
----------------

`Mumble <https://github.com/compomics/mumble>`_ is a PSM generator that expands each PSM's
peptidoform with candidate `Unimod <http://www.unimod.org/>`_ modifications that explain its
precursor mass shift. This is useful when the PSM file comes from a search that already
allows for a wide precursor mass tolerance or an open/semi-open modification search, but did
not (or could not) resolve the mass shift to a specific modification.

For each input PSM, Mumble:

1. Calculates the mass shift between the observed precursor mass and the theoretical mass of
   the (unmodified) peptidoform.
2. Looks up candidate Unimod modifications - optionally combined, and optionally combined with
   amino acid substitutions - that explain this mass shift within a configurable mass error.
3. Generates one new candidate PSM per matching modification (or combination of modifications),
   each with the modification placed at every allowed site.

MS²Rescore then rescores all candidate PSMs (and, depending on configuration, the original
PSM) together, so that the combination of features and the rescoring model can pick the most
likely explanation for the mass shift.

Mumble is integrated into MS²Rescore as a *PSM generator*, configured through the
``psm_generator`` configuration option, alongside the *feature generators* and *rescoring
engine*. See :ref:`Configuration` for how these fit together.


Installing Mumble
------------------

Mumble is an optional dependency. Install it together with MS²Rescore with:

.. code-block:: bash

    pip install ms2rescore[mumble]

or, in an existing MS²Rescore installation:

.. code-block:: bash

    pip install mumble

If Mumble is not installed and ``mumble`` is present in the ``psm_generator`` configuration,
MS²Rescore raises a clear configuration error rather than failing silently.

.. note::
  Mumble is currently only configurable through the JSON/TOML configuration file or the Python
  API, not through the graphical user interface (GUI) or the command line interface (CLI).


Enabling Mumble
----------------

Add a ``mumble`` entry to the ``psm_generator`` configuration section:

.. tab:: JSON

  .. code-block:: json

    "psm_generator": {
      "mumble": {
        "mass_error": 0.02,
        "combination_length": 1,
        "aa_combinations": 0,
        "exclude_mutations": false
      }
    }

.. tab:: TOML

  .. code-block:: toml

    [ms2rescore.psm_generator.mumble]
    mass_error = 0.02
    combination_length = 1
    aa_combinations = 0
    exclude_mutations = false

Only one PSM generator can be configured at a time, and ``mumble`` is currently the only
available option.


Configuration options
----------------------

- ``mass_error`` (default: ``0.02``): Mass error tolerance in Da used to match the precursor
  mass shift to a candidate modification.
- ``combination_length`` (default: ``1``): Maximum number of modifications to combine per mass
  shift. Lower combination lengths are always included as well. Increasing this quickly
  increases runtime and the number of candidate PSMs.
- ``exclude_mutations`` (default: ``false``): Exclude candidate modifications classified as
  amino acid substitutions.
- ``aa_combinations`` (default: ``0``): Number of amino acid combinations to additionally
  consider as a mass-shift explanation, on top of Unimod modifications. Requires
  ``fasta_file``. Increases runtime combinatorially; keep this low (0-2).
- ``fasta_file`` (default: ``null``): Path to a FASTA file with protein sequences, used to
  validate amino acid combination candidates. Required if ``aa_combinations`` is greater
  than 0.
- ``modification_file`` (default: ``null``): Path to a restriction list of Unimod
  modifications to consider (TSV with ``unimod_id``/``name``, and optionally ``residue``,
  columns). Defaults to a curated subset bundled with Mumble. Ignored if
  ``all_unimod_modifications`` is true.
- ``all_unimod_modifications`` (default: ``false``): Consider all Unimod modifications instead
  of the restricted default/``modification_file`` list. Substantially increases runtime and
  the risk of false-positive candidate modifications.

.. caution::
  MS²Rescore always runs Mumble with ``include_original_psm`` and ``include_decoy_psm`` forced
  to ``true``, regardless of what is configured: both the original PSM and modified decoys are
  needed for correct downstream rescoring and FDR control. Setting these two options manually
  has no effect and can be omitted from the configuration.


How Mumble PSMs are filtered
-----------------------------

Adding candidate PSMs for every matching modification can introduce many low-quality
candidates, especially with a larger ``combination_length`` or ``all_unimod_modifications``
enabled. To limit this, MS²Rescore removes any Mumble-generated candidate PSM whose fraction of
matched fragment ions (``matched_ions_pct``) drops more than 25% below that of the original,
unmodified PSM, *before* rescoring. This keeps the candidate set focused on modifications that
are actually well-supported by the fragmentation spectrum, without requiring extra
configuration.

Mumble-generated candidates are also excluded from the DeepLC and IM2Deep calibration sets:
since a candidate is an unconfirmed hypothesis about the peptide's identity (it inherits the
original PSM's score and q-value verbatim, so it cannot be distinguished from a confirmed hit on
those columns alone), only the original, non-Mumble PSMs are used to calibrate retention time and
ion mobility predictions.


.. caution::
  **Ties at the top rank are broken arbitrarily.** With the default ``max_psm_rank_output: 1``,
  MS²Rescore keeps only the single best-scoring PSM per spectrum. If several Mumble candidates
  (and/or the original PSM) end up with an identical rescoring score for the same spectrum - not
  uncommon, since near-isobaric or poorly discriminated modifications can score alike - which one
  is kept as "the" top hit is **not deterministic**: it depends on the internal sort order of the
  PSM list, not on any chemically meaningful tie-break. MS²Rescore logs a warning when this
  situation can occur (Mumble enabled and ``max_psm_rank_output`` is 1). To inspect these cases
  rather than silently pick one, set ``max_psm_rank_output`` to a value greater than 1 and review
  same-spectrum, same-score groups in the output before drawing conclusions from the winning
  peptidoform.


Tips and known limitations
----------------------------

- Start with the default ``modification_file`` (a curated subset of common modifications).
  Only switch to ``all_unimod_modifications`` if you have a specific reason to search the full
  Unimod database, and be aware that this increases runtime and the risk of
  false-positive candidate modifications.
- Keep ``combination_length`` and ``aa_combinations`` low (1 and 0-2, respectively) unless you
  have a specific reason to expect multiple simultaneous modifications or sequence variants;
  both increase runtime combinatorially.
- Because Mumble is beta software, always inspect a subset of the rescoring results (for
  example with the HTML report, see :ref:`Output files`) to confirm that the added
  modifications are chemically sensible before trusting them in downstream analysis.
- If Mumble raises an unexpected error, first try without ``fasta_file``/``aa_combinations`` and
  with the default ``modification_file`` to isolate whether the issue is related to the search
  space configuration, then report the issue (with a minimal reproducing example if possible) on
  the `Mumble issue tracker <https://github.com/compomics/mumble/issues>`_.
- Mumble caches the parsed Unimod modification list on disk between runs. If you suspect the
  cache is stale (for example after switching ``modification_file``), clear it with
  ``mumble --clear-cache``.


Standalone usage
------------------

Mumble can also be used independently of MS²Rescore, both through its own command line
interface and its Python API. See the
`Mumble README <https://github.com/compomics/mumble#readme>`_ for details.
