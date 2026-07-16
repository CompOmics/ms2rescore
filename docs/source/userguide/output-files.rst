############
Output files
############

Depending on the options you choose, the following files will be created. All PSMs, peptides, and
proteins are not yet filtered at any false discovery rate (FDR) level.

Main output files:

+-----------------------------------+----------------------------------------------------------------------------------+
| File                              | Description                                                                      |
+===================================+==================================================================================+
| ``<prefix>.psms.tsv``             | Main output file with rescored PSMs and their new scores                         |
+-----------------------------------+----------------------------------------------------------------------------------+
| ``<prefix>.report.html``          | HTML report with interactive plots showing the results and some quality control  |
|                                   | metrics.                                                                         |
+-----------------------------------+----------------------------------------------------------------------------------+

Log and configuration files:

+--------------------------------------+--------------------------------------------------------------------------------------+
| File                                 | Description                                                                          |
+======================================+======================================================================================+
| ``<prefix>.log.txt``                 | Log file with information about the run                                              |
+--------------------------------------+--------------------------------------------------------------------------------------+
| ``<prefix>.log.html``                | HTML version of the log file                                                         |
+--------------------------------------+--------------------------------------------------------------------------------------+
| ``<prefix>.full-config.json``        | Full configuration file with all the parameters used                                 |
|                                      | as configured in the user-provided configuration file, the command line or graphical |
|                                      | interface, and the default values.                                                   |
+--------------------------------------+--------------------------------------------------------------------------------------+
| ``<prefix>.feature_names.tsv``       | List of the features and their descriptions                                          |
+--------------------------------------+--------------------------------------------------------------------------------------+
| ``<prefix>.intermediate.psms.tsv``   | Created automatically if the process crashes during feature generation or rescoring. |
|                                      | Contains all PSMs with successfully added features up to the crash point. Can be     |
|                                      | used to resume processing with ``-p <prefix>.intermediate.psms.tsv -t tsv``.         |
+--------------------------------------+--------------------------------------------------------------------------------------+

Rescoring result tables (always written): the
post-rescoring score, q-value, and PEP at each identification level, as plain TSV files --
convenient to open directly (e.g. in Excel), unlike the full PSM list output above, which also
includes rescoring features and other provenance data. The protein-level table is only written
if the PSM file provides a ``protein_list``.

Only the post-rescoring ("after") result is written. The pre-rescoring ("before") result isn't
persisted at all: it's fully reconstructable from the main PSM list's provenance data, which is
what the HTML report and ``ms2rescore-report`` use to regenerate it on demand.

If ``max_psm_rank_output`` is set to more than 1, these tables (and the main PSM list output)
hold multiple ranks per spectrum, and their q-values/PEPs are computed treating each rank as an
independent row rather than through proper spectrum competition -- not statistically rigorous
FDR control. ``max_psm_rank_output > 1`` is intended for surfacing ambiguous results (e.g.
multiple candidate peptidoforms per spectrum from Mumble), not for a corrected identification
count.

+-----------------------------------------+--------------------------------------------------------------+
| File                                    | Description                                                  |
+=========================================+==============================================================+
| ``<prefix>.ristretto.psms.tsv``         | PSMs and their rescored score, at PSM-level FDR.             |
+-----------------------------------------+--------------------------------------------------------------+
| ``<prefix>.ristretto.peptidoforms.tsv`` | Peptidoforms and their score, at peptidoform-level FDR.      |
+-----------------------------------------+--------------------------------------------------------------+
| ``<prefix>.ristretto.peptides.tsv``     | Peptides and their score, at peptide-level FDR.              |
+-----------------------------------------+--------------------------------------------------------------+
| ``<prefix>.ristretto.proteins.tsv``     | Proteins and their score, at protein-level FDR.              |
+-----------------------------------------+--------------------------------------------------------------+
| ``<prefix>.ristretto.weights.tsv``      | Feature weights, showing feature usage in the rescoring run. |
+-----------------------------------------+--------------------------------------------------------------+

If rescoring is disabled (``"rescoring": null``) or in DEBUG mode, the following files will also
be written:

+-------------------------------------------------------------+-----------------------------------------------------------+
| File                                                        | Description                                               |
+=============================================================+===========================================================+
| ``<prefix>.pin``                                            | PSMs with all features for rescoring                      |
+-------------------------------------------------------------+-----------------------------------------------------------+
