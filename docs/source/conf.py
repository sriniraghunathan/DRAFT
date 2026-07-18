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
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',        # To generate autodocs
    'sphinx.ext.mathjax',           # autodoc with maths
    'sphinx.ext.napoleon'           # For auto-doc configuration
]

napoleon_google_docstring = False   # Turn off googledoc strings
napoleon_numpy_docstring = True     # Turn on numpydoc strings
napoleon_use_ivar = True         # For maths symbology

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
