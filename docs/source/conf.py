"""Configuration file for the Sphinx documentation builder."""

import os
import sys

sys.path.insert(0, os.path.abspath("../../"))

from ms2rescore import __version__  # noqa: E402

# Project information
project = "ms2rescore"
author = "CompOmics"
github_project_url = "https://github.com/compomics/ms2rescore/"
github_doc_root = "https://github.com/compomics/ms2rescore/tree/main/docs/"
release = __version__

# General configuration
extensions = [
    "nbsphinx",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosectionlabel",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    "sphinxarg.ext",
    "sphinx_inline_tabs",
    "sphinx_rtd_theme",
    "myst_parser",
]
source_suffix = [".rst"]
master_doc = "index"
exclude_patterns = ["_build"]

# Options for HTML output
html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
html_css_files = ["css/custom.css"]
html_js_files = ["js/badge.min.js"]

# Autodoc options
autodoc_default_options = {"members": True, "show-inheritance": True}
autodoc_member_order = "bysource"
autodoc_typehints = "description"
autoclass_content = "init"

# Intersphinx options
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "pandas": ("https://pandas.pydata.org/pandas-docs/stable/", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "plotly": ("https://plotly.com/python-api-reference/", None),
    "psm_utils": ("https://psm-utils.readthedocs.io/en/stable/", None),
}

# nbsphinx options
nbsphinx_execute = "never"


def setup(app):
    config = {  # noqa: F841
        "enable_eval_rst": True,
    }
