############
Output files
############

Depending on the options you choose, the following files will be created. All PSMs, peptides, and
proteins are not yet filtered at any false discovery rate (FDR) level and contain both target and
decoy identifications.

Main output files:

+--------------------------+------------------------------------------------------------------------------------------+
| File                     | Description                                                                              |
+==========================+==========================================================================================+
| ``<prefix>.tsv``         | Main output file with rescored PSMs and their new scores                                 |
+--------------------------+------------------------------------------------------------------------------------------+
| ``<prefix>.report.html`` | HTML report with interactive plots showing the results and some quality control metrics. |
+--------------------------+------------------------------------------------------------------------------------------+

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
| ``<prefix>.intermediate.tsv``        | Created automatically if the process crashes during feature generation or rescoring. |
|                                      | Contains all PSMs with successfully added features up to the crash point. Can be     |
|                                      | used to resume processing with ``-p <prefix>.intermediate.tsv -t tsv``.              |
+--------------------------------------+--------------------------------------------------------------------------------------+

Rescoring result tables: the post-rescoring score, q-value, and PEP at each identification level,
as plain TSV files, convenient to open directly (e.g. in Excel). The protein-level table is only
written if the PSM file provides a ``protein_list``.

If ``max_psm_rank_output`` is set to more than 1, these tables (and the main PSM list output)
hold multiple ranks per spectrum, and their q-values/PEPs are computed treating each rank as an
independent row rather than through proper spectrum competition (not statistically rigorous
FDR control). ``max_psm_rank_output > 1`` is intended for surfacing ambiguous results (e.g.
multiple candidate peptidoforms per spectrum from Mumble), not for a corrected identification
count.

+-------------------------------+-------------------------------------------------------------------------------------------------+
| File                          | Description                                                                                     |
+===============================+=================================================================================================+
| ``<prefix>.psms.tsv``         | PSMs and their rescored score, at PSM-level FDR.                                                |
+-------------------------------+-------------------------------------------------------------------------------------------------+
| ``<prefix>.peptidoforms.tsv`` | Peptidoforms (sequence and modifications, no charge) and their score, at peptidoform-level FDR. |
+-------------------------------+-------------------------------------------------------------------------------------------------+
| ``<prefix>.peptides.tsv``     | Peptides (stripped sequence) and their score, at peptide-level FDR.                             |
+-------------------------------+-------------------------------------------------------------------------------------------------+
| ``<prefix>.proteins.tsv``     | Protein groups and their score, at protein group-level FDR.                                     |
+-------------------------------+-------------------------------------------------------------------------------------------------+
| ``<prefix>.weights.tsv``      | Feature weights, showing feature usage in the rescoring run.                                    |
+-------------------------------+-------------------------------------------------------------------------------------------------+
