# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys

# -- Path setup --------------------------------------------------------------
# The DRAFT package lives at the repo root (../../DRAFT relative to this file).
# Both the repo root and the modules folder are added so autodoc can import
# the code whether or not __init__.py files are present.
sys.path.insert(0, os.path.abspath('../..'))
sys.path.insert(0, os.path.abspath('../../modules'))

# -- Project information -----------------------------------------------------
project = 'DRAFT'
copyright = '2026, Srini Raghunathan'
author = 'Srini Raghunathan'
release = '1.0'
version = '1.0'

# -- General configuration ---------------------------------------------------
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.napoleon',      # numpy/google style docstrings
    'sphinx.ext.viewcode',      # link to highlighted source
    'sphinx.ext.mathjax',       # LaTeX math in docstrings
    'sphinx.ext.intersphinx',
]

templates_path = ['_templates']
exclude_patterns = []

# Generate autosummary stub pages automatically
autosummary_generate = True

# Autodoc settings
autodoc_member_order = 'bysource'
autodoc_default_options = {
    'members': True,
    'undoc-members': True,
    'show-inheritance': True,
}

# Heavy scientific dependencies are mocked so ReadTheDocs can build the
# docs without installing them. Add/remove entries to match your imports.
autodoc_mock_imports = [
    'numpy',
    'scipy',
    'matplotlib',
    'pylab',
    'healpy',
    'camb',
    'pysm3',
    'astropy',
    'h5py',
    'tqdm',
    'pandas',
    'sklearn',
    'colossus',
    'numba',
]

# Intersphinx mapping to external docs
intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
    'scipy': ('https://docs.scipy.org/doc/scipy/', None),
}

# -- Options for HTML output -------------------------------------------------
html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_theme_options = {
    'collapse_navigation': False,
    'navigation_depth': 3,
}
