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

Rescoring engine files:

+-------------------------------------------------------------+-------------------------------------------------------------+
| File                                                        | Description                                                 |
+=============================================================+=============================================================+
| ``<prefix>.<mokapot/percolator>.psms.txt``                  | PSMs and their new scores at PSM-level FDR.                 |
+-------------------------------------------------------------+-------------------------------------------------------------+
| ``<prefix>.<mokapot/percolator>.peptides.txt``              | Peptides and their new scores at peptide-level FDR.         |
+-------------------------------------------------------------+-------------------------------------------------------------+
| ``<prefix>.<mokapot/percolator>.proteins.txt``              | Proteins and their new scores at protein-level FDR.         |
+-------------------------------------------------------------+-------------------------------------------------------------+
| ``<prefix>.<mokapot/percolator>.decoy.psms.txt``            | Decoy PSMs and their new scores at PSM-level FDR.           |
+-------------------------------------------------------------+-------------------------------------------------------------+
| ``<prefix>.<mokapot/percolator>.decoy.peptides.txt``        | Decoy peptides and their new scores at peptide-level FDR.   |
+-------------------------------------------------------------+-------------------------------------------------------------+
| ``<prefix>.<mokapot/percolator>.decoy.proteins.txt``        | Decoy proteins and their new scores at protein-level FDR.   |
+-------------------------------------------------------------+-------------------------------------------------------------+
| ``<prefix>.<mokapot/percolator>.weights.txt``               | Feature weights, showing feature usage in the rescoring run |
+-------------------------------------------------------------+-------------------------------------------------------------+

If no rescoring engine is selected, if Percolator was selected, or in DEBUG mode, the following
files will also be written:

+-------------------------------------------------------------+-----------------------------------------------------------+
| File                                                        | Description                                               |
+=============================================================+===========================================================+
| ``<prefix>.pin``                                            | PSMs with all features for rescoring                      |
+-------------------------------------------------------------+-----------------------------------------------------------+
