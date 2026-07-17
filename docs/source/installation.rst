************
Installation
************

Python package
==============

.. image:: https://flat.badgen.net/badge/install%20with/pip/green?icon=pypi
    :alt: Install with pip
    :target: https://pypi.org/project/ms2rescore/

.. image:: https://flat.badgen.net/badge/install%20with/conda/green?icon=conda
    :alt: Install with conda
    :target: https://anaconda.org/bioconda/ms2rescore

MS²Rescore is installable as a Python package on Windows, macOS and Linux.

In a fresh `virtual environment <https://docs.python.org/3/library/venv.html>`_, run::

    pip install ms2rescore


Or, in a fresh `conda environment <https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html>`_, run::

    conda install -c bioconda ms2rescore

Bioconda packages are only available for Linux and macOS.


Optional dependencies
=====================

Some features require additional dependencies that are not installed by default. Install them
with the relevant extra, for example:

.. code-block:: bash

    pip install ms2rescore[mumble]

- ``mumble``: enables the `Mumble <https://github.com/compomics/mumble>`_ PSM generator, which
  proposes candidate modifications for PSMs with an unresolved precursor mass shift. Mumble is
  **beta** software; see :ref:`Open modification searching with Mumble` before enabling it.
- ``idxml``: enables reading OpenMS idXML PSM files.


Windows installer
=================

.. image:: https://flat.badgen.net/badge/install%20for/windows/blue?icon=windows
    :alt: Get for Windows
    :target: https://github.com/compomics/ms2rescore/releases/latest

Download the ``.exe`` file from the
`latest release <https://github.com/compomics/ms2rescore/releases/latest>`_
and go through the installation steps. If Microsoft Defender SmartScreen displays a warning, click
"More info" and then click "Run anyway".


Docker container
================

.. image:: https://flat.badgen.net/badge/pull/biocontainer/blue?icon=docker
    :alt: Pull with Docker
    :target: https://quay.io/repository/biocontainers/ms2rescore

First check the latest version tag on
`biocontainers/ms2rescore/tags <https://quay.io/repository/biocontainers/ms2rescore?tab=tags>`_.
Then pull and run the container with:

.. code-block:: bash

   docker container run -v <working-directory>:/data -w /data quay.io/biocontainers/ms2rescore:<tag> ms2rescore <ms2rescore-arguments>

where ``<working-directory>`` is the absolute path to the directory with your MS²Rescore input
files, ``<tag>`` is the container version tag, and ``<ms2rescore-arguments>`` are the ms2rescore
command line options (see :ref:`Command line interface`).


For development
===============

Clone this repository and use pip to install an editable version:

.. code-block:: bash

   pip install --editable .
