.. _mumble:

Open modification searching with Mumble
========================================

.. warning::
  Mumble support is **beta** software. Inspect a subset of the results before using them in
  downstream analysis, and report problems on the
  `MS²Rescore <https://github.com/compomics/ms2rescore/issues>`_ or
  `Mumble <https://github.com/compomics/mumble/issues>`_ issue tracker, together with the
  MS²Rescore version and the ``<output>.log.txt`` file.


What Mumble does
----------------

An open or wide-tolerance search identifies the peptide of a spectrum, but often leaves a
precursor mass shift unexplained. `Mumble <https://github.com/compomics/mumble>`_ proposes
explanations for that shift and MS²Rescore decides between them.

For each PSM, Mumble:

1. computes the mass shift between the observed precursor and the unmodified peptide;
2. looks up every `Unimod <http://www.unimod.org/>`_ modification, combination of
   modifications or amino acid addition that explains the shift within ``mass_error``;
3. writes one candidate PSM per explanation and per allowed residue (the *site*);
4. optionally writes *decoy sites*: the same modification on a residue it cannot occupy.

MS²Rescore computes features for all candidates, rescores them together with the original
PSMs, and, with ``rank_sites``, reports per spectrum which explanation and which site the
spectrum supports, and how confident that call is.

Mumble is a *PSM generator* (``psm_generator`` in the configuration), next to the feature
generators and the rescoring engine (see :ref:`Configuration`). It is configured through the
JSON or TOML configuration file or the Python API, not through the command line or the
graphical user interface.


Installation
------------

Mumble is an optional dependency of MS²Rescore:

.. code-block:: bash

    pip install --pre "ms2rescore[mumble]"

``--pre`` installs the beta release. A CUDA GPU is optional and speeds up DeepLC and IM2Deep.
When ``mumble`` is configured but not installed, MS²Rescore stops with a configuration error.


Quick start
-----------

Input:

- a PSM file from an open or wide-tolerance search, with its decoy PSMs (for example Sage
  with a precursor window of several hundred Da);
- the matching spectrum files (mzML or MGF) in ``spectrum_path``;
- a ``modification_mapping`` from the search engine modification labels to Unimod names.

Configuration for modification localisation:

.. tab:: JSON

  .. code-block:: json

    {
      "ms2rescore": {
        "feature_generators": {
          "basic": {},
          "ms2pip": { "model": "HCD" },
          "ms2": { "add_mod_info": true },
          "deeplc": { "finetune": true },
          "im2deep": {}
        },
        "psm_generator": {
          "mumble": {
            "mass_error": 0.02,
            "combination_length": 1,
            "include_mumble_decoys": true,
            "isotope_errors": [0, 1]
          }
        },
        "rescoring": { "model": "svm", "train_fdr": 0.01 },
        "rank_sites": true,
        "max_psm_rank_output": 10,
        "psm_file_type": "sage",
        "spectrum_path": "spectra/",
        "id_decoy_pattern": "^rev_",
        "modification_mapping": {
          "+79.9663": "Phospho",
          "+15.9949": "Oxidation",
          "+57.0215": "Carbamidomethyl"
        },
        "processes": 8
      }
    }

.. tab:: TOML

  .. code-block:: toml

    [ms2rescore]
    rank_sites = true
    max_psm_rank_output = 10
    psm_file_type = "sage"
    spectrum_path = "spectra/"
    id_decoy_pattern = "^rev_"
    processes = 8

    [ms2rescore.feature_generators]
    basic = {}
    ms2pip = { model = "HCD" }
    ms2 = { add_mod_info = true }
    deeplc = { finetune = true }
    im2deep = {}

    [ms2rescore.psm_generator.mumble]
    mass_error = 0.02
    combination_length = 1
    include_mumble_decoys = true
    isotope_errors = [0, 1]

    [ms2rescore.rescoring]
    model = "svm"
    train_fdr = 0.01

    [ms2rescore.modification_mapping]
    "+79.9663" = "Phospho"
    "+15.9949" = "Oxidation"
    "+57.0215" = "Carbamidomethyl"

Set ``ms2pip.model`` to the instrument and fragmentation of the data (for example
``timsTOF`` for timsTOF data) and ``id_decoy_pattern`` to the decoy prefix of the FASTA file.
Then run:

.. code-block:: bash

    ms2rescore -c config.json -p results.sage.tsv -s spectra/ -o out/results

``max_psm_rank_output`` sets how many candidates per spectrum are written. With ``rank_sites``,
the candidates are written in ranker order, so ``1`` keeps only the reported candidate.


Configuration options
----------------------

Options of the ``mumble`` PSM generator (``psm_generator.mumble``):

=============================  ===========  ===========================================================
Option                         Default      Description
=============================  ===========  ===========================================================
``mass_error``                 ``0.02``     Tolerance in Da to match the mass shift to an explanation.
``combination_length``         ``1``        Maximum number of modifications combined per mass shift;
                                            lower numbers are included.
``isotope_errors``             ``[0]``      Precursor isotope errors (13C spacings, 1.00336 Da) tried
                                            when matching a mass shift. ``[0, 1]`` also finds a
                                            modification when the first 13C peak was selected as
                                            precursor. Stored in ``metadata["isotope_error"]``.
``include_mumble_decoys``      ``false``    Place each single modification also on residues it cannot
                                            occupy, one decoy site per real site. A decoy site with the
                                            same residue masses as a real candidate is not generated.
                                            Required for ``rank_sites``.
``modification_file``          bundled      Restriction list of Unimod modifications (TSV with
                                            ``unimod_id`` or ``name`` and optionally ``residue``).
``all_unimod_modifications``   ``false``    Use all of Unimod instead of the restriction list.
``exclude_mutations``          ``false``    Exclude Unimod entries classified as amino acid
                                            substitution.
``aa_combinations``            ``0``        Number of amino acids that may be added to the peptide as
                                            an explanation. Requires ``fasta_file``.
``fasta_file``                 ``null``     FASTA file to check that an extended peptide occurs in a
                                            protein of the PSM.
=============================  ===========  ===========================================================

MS²Rescore always runs Mumble with ``include_original_psm`` and ``include_decoy_psm`` set to
``true``; setting these options has no effect.

Related MS²Rescore options:

===================================  ==========  ====================================================
Option                               Default     Description
===================================  ==========  ====================================================
``rank_sites``                       ``false``   Rank the candidates of each spectrum and write the
                                                 localisation outcome (see below).
``mumble_matched_ions_threshold``    ``0.5``     Remove Mumble candidates whose fraction of matched
                                                 fragment ions is below this fraction of the
                                                 original PSM's. ``0`` keeps every candidate.
``max_psm_rank_output``              ``1``       Number of candidates per spectrum in the output.
``feature_generators.ms2``           (none)      ``add_mod_info: true`` adds the
                                                 modification-aware features below.
``tolerance_value``                  ``0.02``    Fragment tolerance in Da for annotation, site sets
                                                 and fragment-identical explanations.
===================================  ==========  ====================================================


How it works
------------

1. **Read PSMs.** Modification labels are mapped to Unimod names.
2. **Generate candidates.** Mumble writes the original PSM, every explanation at every allowed
   site and, with ``include_mumble_decoys``, the decoy sites. Decoy peptides from the search
   get candidates as well.
3. **Compute features.** All feature generators run on all candidates, decoy sites included.
   DeepLC and IM2Deep are calibrated on original PSMs only.
4. **Set decoy sites aside.** Decoy sites are not rescored and do not take part in the PSM,
   peptide or protein FDR.
5. **Filter.** Candidates below ``mumble_matched_ions_threshold`` are removed.
6. **Rescore.** The rescoring model scores every candidate; q-values and PEPs are computed
   with one identification per spectrum.
7. **Rank** (``rank_sites``). Per spectrum, in two steps, on all candidates:

   a. *Explanation.* The original PSM and each modification set, each represented by its
      best-scoring candidate, are ranked against each other. When explanations give the same
      fragment masses (for example Phospho and Sulfo on the same residue), the explanation that
      most often wins with a clear gain over the original PSM in the same run is put first.
   b. *Site.* Within each modification set, the candidates are ranked with the decoy sites as
      known negatives.

   Both rankers learn from the differences between candidates of the same spectrum, so
   spectrum-level features (retention time of the spectrum, charge) do not count. Precursor
   mass features and the search engine score are not used.
8. **Localisation statistics.** See below.
9. **Output.** Within each spectrum, the rescoring scores, q-values and PEPs are reassigned in
   ranker order, so the spectrum keeps its FDR status and the best-supported candidate is
   reported. The top ``max_psm_rank_output`` candidates per spectrum are written.


Reading the localisation output
--------------------------------

With ``rank_sites``, every PSM gets the columns below. The localisation columns describe the
spectrum and are the same for all its candidates.

==================================  ============================================================
Column                              Meaning
==================================  ============================================================
``meta:site_rank``                  Rank of the candidate in its spectrum, over all candidates
                                    including decoy sites. 1 is the best.
``meta:site_score``                 Site-ranker log-odds; comparable within a spectrum only.
``meta:localisation_status``        Outcome of the spectrum (see the next table).
``meta:localisation_set``           Sites the observed fragments cannot tell apart from the top
                                    site, as peptidoforms separated by ``;``.
``meta:localisation_alternatives``  Other modifications with the same fragment masses as the top
                                    candidate, for example Sulfo next to Phospho.
``meta:localisation_pep``           Probability that the true site is not in the set.
``meta:modification_pep``           Probability that the modification explains the spectrum no
                                    better than the original PSM.
``meta:localisation_confidence``    (1 - ``modification_pep``) x (1 - ``localisation_pep``).
``meta:localisation_qvalue``        Expected fraction of wrong calls among the target spectra
                                    with at least this confidence.
==================================  ============================================================

Status values:

======================  =====================================================================
Status                  Meaning
======================  =====================================================================
``confident``           One site, confidence at least 0.5.
``ambiguous``           The modification is supported, but the fragments cannot separate two
                        or more sites; all are in ``localisation_set``.
``original``            The original (unmodified) PSM is the best explanation.
``decoy``               A decoy site scores higher than every real site of the modification.
``no_gain``             The modification does not explain the spectrum better than the
                        original PSM (``modification_pep`` at least 0.5).
``low_confidence``      Confidence below 0.5.
``unranked``            The spectrum has no competing candidates.
======================  =====================================================================

How the numbers are computed:

- **Site set.** For the top candidate and each other candidate of the same modification, the
  site-determining ions are the b and y ions (charge 1 to precursor charge minus 1) that differ
  in mass between the two. A candidate joins the set when the spectrum contains no more
  site-determining ions of the top candidate than of that candidate.
- **localisation_pep.** Target-decoy competition: the margin of the site set over the best
  candidate outside it is compared with the margin of the best decoy site over its competitors,
  over all spectra of the run, and converted into a posterior error probability with an
  isotonic fit. Decoy sites with the same fragment masses as a real site do not take part.
- **modification_pep.** The hyperscore gain of the top candidate over the original PSM is
  compared with the same gain in spectra of decoy peptides, which have no true modification.
- **localisation_qvalue.** Running mean of (1 - ``localisation_confidence``) over the target
  spectra sorted by confidence.

Filtering:

- Keep the spectra with status ``confident`` or ``ambiguous`` and a ``localisation_qvalue``
  below the chosen threshold, for example 0.01 or 0.05. The q-value controls the localisation
  error; the PSM ``qvalue`` controls the peptide identification.
- Report an ambiguous spectrum as one identification with several possible sites, not as
  several identifications.
- ``localisation_alternatives`` lists modifications that the fragments cannot separate from the
  reported one. The reported modification is the one that wins most often in the run;
  retention time and precursor mass can still separate such pairs.
- In a spectrum where a decoy site wins, the best returned candidate has ``site_rank`` 2 or
  higher, because decoy sites are not written.

The log reports the share of target spectra per status and the largest feature weights of both
rankers. The weights are written to ``<output>.explanation_weights.tsv`` and
``<output>.site_weights.tsv``: one row per feature, one column per cross-validation fold and
their mean, in the within-spectrum centred and standardised feature space.


Modification-aware features
----------------------------

``add_mod_info`` in the ``ms2`` feature generator adds six features. Unmodified PSMs get 0.

==================================  ===============================================================
Feature                             Meaning
==================================  ===============================================================
``mod_loss_intensity_ratio``        Intensity fraction of backbone fragments that lost a neutral
                                    loss of the modification (for example b5-H3PO4).
``precursor_mod_loss_ratio``        Intensity fraction of the precursor after a neutral loss of the
                                    modification.
``delta_hyperscore_unmod``          Hyperscore of the PSM minus the hyperscore with one
                                    modification removed (a substitution reverted to the residue of
                                    the original PSM); minimum over the modifications.
``mod_site_flank_matched``          Fraction of the four ions next to the site (b and y on both
                                    sides of the modified residue) that is matched; minimum over the
                                    modifications.
``mod_site_flank_intensity_ratio``  Intensity fraction of those ions; minimum over the
                                    modifications.
``mod_mass_error_offset_ppm``       Difference between the median ppm error of fragments that
                                    contain a modified residue and that of the other fragments.
==================================  ===============================================================

The neutral losses of a modification are taken from the PSMs that carry it on an allowed
residue, and are then matched on any residue the modification is placed on, decoy sites
included. Losses of the whole modification and losses that occur without the modification
(for example H2O) are not counted. The features need modification names that can be parsed;
PSMs with numeric mass shifts only get backbone evidence, and the log reports how many.

The ``basic`` feature generator adds ``isotope_error``: the number of 13C spacings between the
selected precursor peak and the monoisotopic mass, taken from Mumble or the search engine. The
precursor mass error features are computed against that peak.


Known limitations
------------------

- Sequence extensions (``aa_combinations``) have no decoy sites, and decoy peptides do not get
  extensions. Their localisation and their PSM-level FDR are therefore not controlled.
- Rare modifications compete with the original PSM and with frequent modifications on equal
  terms; a rare modification is reported when its fragments support it, and is then often
  labelled ``no_gain`` or ``decoy``.
- Mumble only adds mass. A negative shift from a peptide that is too long is not explained.
- With ``rename_to_usi``, decoy sites keep the native spectrum identifier.
- Amino acid substitutions have no decoy sites; only ``modification_pep`` applies to them.
- Runs with one spectrum file give DeepLC and IM2Deep little data for calibration; the log
  warns about this.
- Memory use grows with the number of candidates. Mumble can multiply the number of PSMs
  several times, so memory use can be much higher than for a run without Mumble.


Tips
----

- Start with the default ``modification_file``. ``all_unimod_modifications`` enlarges the
  candidate set and the runtime considerably.
- Keep ``combination_length`` at 1 and ``aa_combinations`` at 0 to 2; both increase the
  runtime combinatorially.
- Mumble caches the parsed modification list. Clear it with ``mumble --clear-cache`` after
  changing ``modification_file``.
- When Mumble fails, run again without ``fasta_file`` and ``aa_combinations`` and with the
  default ``modification_file`` to find the cause.


Standalone usage
------------------

Mumble also runs on its own, through its command line and Python API. See the
`Mumble README <https://github.com/compomics/mumble#readme>`_.
