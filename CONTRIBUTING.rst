############
Contributing
############

This document briefly describes how to contribute to
`ms2rescore <https://github.com/compomics/ms2rescore>`_.



Before you begin
################

If you have an idea for a feature, use case to add or an approach for a bugfix,
you are welcome to communicate it with the community by opening a
thread in
`GitHub Discussions <https://github.com/compomics/ms2rescore/discussions>`_
or in `GitHub Issues <https://github.com/compomics/ms2rescore/issues>`_.

Welcome contributions include:

- New features, such as the addition of new feature generators
- Improvements of existing functionality
- Bugfixes



Development setup
#################

Local install
*************

#. Setup Python 3, and preferably create a virtual environment.
#. Clone the `ms2rescore repository <https://github.com/compomics/ms2rescore>`_.
#. Use pip in editable mode to setup the development environment:

.. code-block:: sh

    pip install --editable .[dev,docs]


Pre-commit hooks
****************

Pre-commit hooks ensure that certain checks are performed before making a new commit. For instance,
the ``black`` pre-commit hook is used to format all Python code, and ``jsonschema2md`` is used to
automatically generate Markdown documentation for the configuration file. Setup the pre-commit
hooks with:

.. code-block:: sh

    pre-commit install


Unit tests
**********

Run tests with ``pytest``:

.. code-block:: sh

    pytest ./tests


Documentation
*************

To work on the documentation and get a live preview, install the requirements
and run ``sphinx-autobuild``:

.. code-block:: sh

    pip install .[docs]
    sphinx-autobuild  --watch ./ms2rescore ./docs/source/ ./docs/_build/html/

Then browse to http://localhost:8000 to watch the live preview.


How to contribute
#################

- Fork `ms2rescore <https://github.com/compomics/ms2rescore>`_ on GitHub to
  make your changes.
- Commit and push your changes to your
  `fork <https://help.github.com/articles/pushing-to-a-remote/>`_. Commit messages
  should follow the `Conventional Commits <https://www.conventionalcommits.org/>`_
  format (e.g. ``fix: ...``, ``feat: ...``, ``docs: ...``).
- Ensure that the tests and documentation (both Python docstrings and files in
  ``/docs/source/``) have been updated according to your changes. Python
  docstrings are formatted in the
  `numpydoc style <https://numpydoc.readthedocs.io/en/latest/format.html>`_.
- Update ``CHANGELOG.md`` as part of your pull request: add an entry describing your
  change under the ``[Unreleased]`` section, following
  `Keep a Changelog <https://keepachangelog.com/en/1.0.0/>`_.
- Code should be formatted with `ruff <https://docs.astral.sh/ruff/>`_, using the
  configuration in ``pyproject.toml``. Run ``ruff format .`` and ``ruff check .``
  before submitting your pull request.
- Open a
  `pull request <https://help.github.com/articles/creating-a-pull-request/>`_
  with these changes. You pull request message ideally should include:

    - A description of why the changes should be made.
    - A description of the implementation of the changes.
    - A description of how to test the changes.

- The pull request should pass all the continuous integration tests which are
  automatically run by
  `GitHub Actions <https://github.com/compomics/ms2rescore/actions>`_.



Release workflow
################

- When a new version is ready to be published:

    #. Change the ``__version__`` in ``ms2rescore/__init__.py`` following
       `semantic versioning <https://semver.org/>`_.
    #. Update the changelog (if not already done) in ``CHANGELOG.md`` according to
       `Keep a Changelog <https://keepachangelog.com/en/1.0.0/>`_.
    #. Merge all final changes with the ``main`` branch.
    #. On GitHub, draft a new release with the new version number and the
       changes that are listed in ``CHANGELOG.md``.

- When a new release is published on GitHub, the following GitHub Actions are triggered:

    #. The Python package is build and published to PyPI.
    #. The Windows installer is build with pyInstaller and InnoSetup and published to the GitHub
       release.

- A webhook triggers a new build of the documentation on Read the Docs.

- The Bioconda recipe is automatically updated by the Bioconda bot, and subsequently both the Conda
  Python package and the Docker image are build.
