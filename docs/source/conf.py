"""Sphinx configuration for the cycombinepy documentation site."""

from __future__ import annotations

import os
import sys
from datetime import datetime

# Make the package importable for autodoc when building from a clean checkout.
sys.path.insert(0, os.path.abspath("../../src"))

import cycombinepy  # noqa: E402

project = "cycombinepy"
author = "cycombinepy contributors"
copyright = f"{datetime.now():%Y}, {author}"
release = cycombinepy.__version__
version = ".".join(release.split(".")[:2])

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    "sphinx.ext.viewcode",
    "sphinx_autodoc_typehints",
    "myst_nb",
]

templates_path = ["_templates"]
exclude_patterns: list[str] = ["_build", "Thumbs.db", ".DS_Store"]
# myst_nb registers parsers for .md and .ipynb; .rst is handled by Sphinx itself.
source_suffix = {
    ".rst": "restructuredtext",
    ".md": "myst-nb",
    ".ipynb": "myst-nb",
}

# -- Autodoc / autosummary -------------------------------------------------
autosummary_generate = True
autodoc_default_options = {
    "members": True,
    "undoc-members": False,
    "show-inheritance": True,
}
autodoc_typehints = "description"
napoleon_numpy_docstring = True
napoleon_google_docstring = False

# -- MyST-NB ---------------------------------------------------------------
# Tutorial notebooks are not executed in docs builds. Most render as source-only,
# while selected executed benchmarks store outputs directly in the notebook.
# CI exercises deterministic docs snippets separately.
nb_execution_mode = "off"
myst_enable_extensions = [
    "colon_fence",
    "dollarmath",
    "deflist",
    "smartquotes",
]
myst_heading_anchors = 3

# -- Intersphinx -----------------------------------------------------------
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable", None),
    "pandas": ("https://pandas.pydata.org/docs", None),
    "scipy": ("https://docs.scipy.org/doc/scipy", None),
    "anndata": ("https://anndata.readthedocs.io/en/stable", None),
    "scanpy": ("https://scanpy.readthedocs.io/en/stable", None),
}

# -- HTML output -----------------------------------------------------------
html_theme = "furo"
html_title = f"cyCombinePy {release}"
html_static_path = ["_static"]
html_theme_options = {
    "source_repository": "https://github.com/mdmanurung/cyCombinePy",
    "source_branch": "main",
    "source_directory": "docs/source/",
}

# Don't fail on missing cross-references to optional deps not in intersphinx.
nitpicky = False
