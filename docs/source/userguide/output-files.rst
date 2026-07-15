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

Rescoring result tables (written if ``write_rescoring_tables`` is enabled, the default),
one pair of before/after files per identification level. The peptide- and protein-level tables are
only written if the PSM file provides the corresponding information (a stripped peptide sequence
and a ``protein_list``, respectively):

+----------------------------------------------------+------------------------------------------------------------------+
| File                                               | Description                                                      |
+====================================================+==================================================================+
| ``<prefix>.ristretto.psms_before.parquet``         | PSMs and their pre-rescoring score, at PSM-level FDR.            |
+----------------------------------------------------+------------------------------------------------------------------+
| ``<prefix>.ristretto.psms_after.parquet``          | PSMs and their new scores at PSM-level FDR.                      |
+----------------------------------------------------+------------------------------------------------------------------+
| ``<prefix>.ristretto.peptidoforms_before.parquet`` | Peptidoforms and their pre-rescoring score, at peptidoform-level |
|                                                    | FDR.                                                             |
+----------------------------------------------------+------------------------------------------------------------------+
| ``<prefix>.ristretto.peptidoforms_after.parquet``  | Peptidoforms and their new scores at peptidoform-level FDR.      |
+----------------------------------------------------+------------------------------------------------------------------+
| ``<prefix>.ristretto.peptides_before.parquet``     | Peptides and their pre-rescoring score, at peptide-level FDR.    |
+----------------------------------------------------+------------------------------------------------------------------+
| ``<prefix>.ristretto.peptides_after.parquet``      | Peptides and their new scores at peptide-level FDR.              |
+----------------------------------------------------+------------------------------------------------------------------+
| ``<prefix>.ristretto.proteins_before.parquet``     | Proteins and their pre-rescoring score, at protein-level FDR.    |
+----------------------------------------------------+------------------------------------------------------------------+
| ``<prefix>.ristretto.proteins_after.parquet``      | Proteins and their new scores at protein-level FDR.              |
+----------------------------------------------------+------------------------------------------------------------------+
| ``<prefix>.ristretto.weights.parquet``             | Feature weights, showing feature usage in the rescoring run.     |
+----------------------------------------------------+------------------------------------------------------------------+

If rescoring is disabled (``"rescoring": null``) or in DEBUG mode, the following files will also
be written:

+-------------------------------------------------------------+-----------------------------------------------------------+
| File                                                        | Description                                               |
+=============================================================+===========================================================+
| ``<prefix>.pin``                                            | PSMs with all features for rescoring                      |
+-------------------------------------------------------------+-----------------------------------------------------------+
