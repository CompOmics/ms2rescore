.. _tims2rescore:

TIMS²Rescore
============

Introduction
------------

`TIMS²Rescore` is a specialized version of `MS²Rescore` for timsTOF DDA-PASEF data. This guide
provides an overview of how to use TIMS²Rescore effectively.

Installing TIMS²Rescore
-----------------------

TIMS²Rescore is part of the ``ms2rescore`` package. Check out the :ref:`installation` instructions
to get started.

Usage
-----

To use TIMS²Rescore, follow these steps:

1. Prepare your input files:
    - To boost DDA-PASEF peptide identifications, TIMS²Rescore requires the spectrum files from
      the timsTOF instrument and the PSM files with identifications from a supported search engine.
    - Make sure that the PSM file format comes from a supported search engine or is a standard
      format such as mzIdentML (See
      :external+psm_utils:ref:`supported file formats <supported file formats>`).
    - Spectrum files can directly be passed as ``.d`` or `miniTDF` raw data or can optionally be
      first converted to mzML or MGF. We recommend using the format that was passed to the search
      engine.

2. Run ``tims2rescore``:
    - Open a terminal or command prompt.
    - Navigate to the directory where your input files are located.
    - Execute the following command:

      .. code-block:: bash

          tims2rescore -p <path_to_psm_file> -s <path_to_spectrum_file>

    Replace `<path_to_psm_file>`, `<path_to_tims_file>`, and `<path_to_output_file>` with the
    actual paths to your input and output files.

    .. admonition:: note

        By default, specialized timsTOF models will be used for predictions. Optionally you can
        further configure TIMS²Rescore through a configuration file. For more information, refer
        to the :ref:`configuration` tab in the user guide.

3. Review the results:
    - Once the ``tims2rescore`` process completes, you will find the rescoring results in the
      same directory as the input files.
    - If you want a detailed report of the rescoring performance, you can either give the set
      `write_report` to `True` in the configuration file, use the `--write_report` option in the
      ``tims2rescore`` command line. Alternatively, run the following command after rescoring:

      .. code-block:: bash

          ms2rescore-report <path_to_psm_file>

      Replace `<path_to_psm_file>` with the path to the main PSM list output file (e.g.,
      ``identifications.ms2rescore.tsv``).

Additional options
------------------

`tims2rescore` provides additional options to customize rescoring. You can explore these options
by running the following command:

.. code-block:: bash

    tims2rescore --help


