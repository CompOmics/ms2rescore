************************
Graphical user interface
************************


Installation
============

The MS²Rescore desktop application can be installed on Windows with a
:ref:`one-click installer <Windows installer>`. Alternatively, or on other platforms, follow the
:ref:`Python package installation instructions <Python package>`.


Starting the application
========================

If installed with the one-click installer, simply start MS²Rescore from the start menu or with the
desktop shortcut. Otherwise, start the application from the
:ref:`command line <command line interface>` with the command ``ms2rescore-gui`` or with
``python -m ms2rescore.gui``.


Application overview
====================

The MS²Rescore graphical user interface is divided into three main sections:

1. A side bar with references, window controls, and the current version number.
2. The configuration pane with input file selection, and parameter configuration.
3. The application log pane with the status output.

On the bottom of the window, the application log level can be selected. The log level determines
which messages are shown in the application log pane. On the bottom right, the application can be
started with the "Start" button. The "Stop" button can be used to stop the application at any time
during the execution.

.. figure:: ../_static/img/gui-overview.png
   :width: 100%
   :alt: MS²Rescore graphical user interface

   Overview of the MS²Rescore desktop application.


Configuring MS²Rescore
======================

Input file selection
^^^^^^^^^^^^^^^^^^^^

The main input for MS²Rescore are the PSM file(s) (search engine output) and the spectrum file(s).
See :ref:`Input files` for more information.

One or more PSM files can be selected from the file system with the "Browse files" button under.
To make ensure correct reading of the file, specify the file type with from the drop-down menu.

.. figure:: ../_static/img/gui-example-xtandem-psm-file.png
   :width: 60%
   :alt: PSM file selection

   PSM file selection


.. figure:: ../_static/img/gui-example-xtandem-psm-filetype.png
   :width: 60%
   :alt: PSM file type selection

   PSM file type selection


To select a single spectrum file (mzML or MGF), click the "Browse files" button. To select a
folder with spectrum files, click the "Browse directories" button.

.. figure:: ../_static/img/gui-example-xtandem-spectra.png
   :width: 60%
   :alt: Spectrum file selection

   Spectrum file selection


Number of processes
^^^^^^^^^^^^^^^^^^^

The number of processes can be configured to run the application in parallel. The default is to
use all available CPU cores. The number of processes can be reduced to avoid overloading the
system or to avoid memory issues. A number under 16 is recommended.


Modification mapping
^^^^^^^^^^^^^^^^^^^^

Depending on the search engine, the peptide modification labels will have to be mapped
to labels that can be understood by MS²Rescore. For example, X!Tandem uses mass shift labels, such
as ``+57.02146`` for carbamidomethylation. However, tools such as DeepLC requires the atomic
composition for all modifications. As this cannot be derived from the mass shift (or other labels
that are not known to MS²Rescore), a mapping has to be provided.

.. figure:: ../_static/img/gui-example-xtandem-modifications-before.png
   :width: 70%
   :alt: Modification mapping

   Modification mapping configuration. Click the plus sign to add more rows.


In modification mapping, click the plus sign to add more rows to the table, or click the minus sign
to remove rows. In the first column "Search engine label", enter the modification label as it
appears in the PSM file. In the second column "ProForma label", enter a ProForma-compatible
modification label. More information on accepted labels can be found in :ref:`Parsing modification
labels`.

.. figure:: ../_static/img/gui-example-xtandem-modifications-filled.png
   :width: 70%
   :alt: Modification mapping

   Modification mapping configuration for the X!Tandem example. Mass shift labels from X!Tandem
   are mapped to ProForma UniMod labels.


Fixed modifications
^^^^^^^^^^^^^^^^^^^

If the search engine PSM file does not contain information on which fixed modifications were used,
this must be specified in the MS²Rescore configuration. At the time of writing, only MaxQuant
``msms.txt``` files do not contain this information. For all other search engines, this information
is contained in the PSM file and the following field can be left empty.


Advanced options
^^^^^^^^^^^^^^^^

Most advanced options are only required for specific use cases or with specific search engine PSM
files. All options are listed in the :doc:`userguide/configuration` section of the user guide.

In the X!Tandem example, only the `PSM ID regex pattern` option is required. This option is used
to extract the spectrum ID from the PSM file. The spectrum ID is used to match the PSM to the
spectrum file. See :ref:`Mapping PSMs to spectra` for more information.

.. figure:: ../_static/img/gui-example-xtandem-advanced.png
   :width: 70%
   :alt: Advanced options

   Advanced options


For reference, all parameters for the X!Tandem example are also listed in the example
configuration file on
`GitHub <https://github.com/compomics/ms2rescore/blob/main/examples/xtandem-ms2rescore.toml>`_.


Starting the rescoring process
==============================

After the configuration is complete, click the "Start" button to start the rescoring process.
The application will show the progress in the application log pane. The log level can be changed
before the run to show more or less information.

.. figure:: ../_static/img/gui-example-xtandem-progress.png
   :width: 100%
   :alt: Running application

   Running application with log output


A pop up will appear when the application is finished, or when an error occurred. If an error
has occurred, the error message in the pop up should provide some insight into what went wrong.
If the error message is not clear, please report the issue on the
`GitHub issue tracker <https://github.com/compomics/ms2rescore/issues>`_ or post your question on
the `Discussion forum <https://github.com/compomics/ms2rescore/discussions>`_.

.. figure:: ../_static/img/gui-example-xtandem-finished.png
   :width: 40%
   :alt: Pop up when MS²Rescore is finished

   Pop up when MS²Rescore is finished


Viewing the results
===================

After a successful run, the output files can be found in the directory of the input PSM file, or
in the specified output directory. The most important files are the ``*.ms2rescore.psms.tsv`` file,
which contains all PSMs with their new scores, and the ``*.ms2rescore.report.html`` file, which
contains interactive charts that visualize the results and various quality control metrics. See
:ref:`Output files` for more information.

.. figure:: ../_static/img/gui-example-xtandem-output-files.png
   :width: 100%
   :alt: Output files

   Overview of the output files after rescoring the X!Tandem example.

Double click the ``*.ms2rescore.report.html`` file to open it in the default web browser:

.. figure:: ../_static/img/qc-reports.png
   :width: 100%
   :alt: Rescoring report

   Rescoring QC report with interactive charts.
