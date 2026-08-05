# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys

# -- Path setup --------------------------------------------------------------
# The code lives at the repository root (../.. from this file) and in
# ../../modules. Both are added to sys.path so that autodoc can import the
# driver and the modules whether or not __init__.py files are present.
# (There is currently no installable package).
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
    'sphinx.ext.autodoc',           # To generate autodocs
    'sphinx.ext.mathjax',           # autodoc with maths
    'sphinx.ext.napoleon',          # For auto-doc configuration
    'sphinx.ext.intersphinx',       # Link numpy/scipy/matplotlib references
    'sphinx.ext.extlinks'           # Shorthand for arXiv links
]

napoleon_google_docstring = False   # Turn off googledoc strings
napoleon_numpy_docstring = True     # Turn on numpydoc strings
napoleon_use_ivar = True            # Use :ivar: for attributes
napoleon_use_rtype = False          # Use inline parenthesis for return type
autodoc_member_order = 'bysource'   # Follow the grouping in each module

intersphinx_mapping = {
    'numpy': ('https://numpy.org/doc/stable/', None),
    'scipy': ('https://docs.scipy.org/doc/scipy/', None),
    'matplotlib': ('https://matplotlib.org/stable/', None),
}

nitpicky = True                     # Fail loudly on unresolved references
nitpick_ignore_regex = [
    ('py:class', r'optional'),      # numpydoc type qualifiers, not classes
    ('py:class', r'array_like'),
    ('py:class', r'ndarray'),
    ('py:class', r'sequence'),
    ('py:class', r"\{?'\w+'\}?"),   # literal-choice sets, e.g. {'qu', 'eb'}
]

# templates_path = ['_templates']   # unused
exclude_patterns = ['build', 'Thumbs.db', '.DS_Store']

extlinks = {'arxiv': ('https://arxiv.org/abs/%s', 'arXiv:%s')}

rst_prolog = """
.. |paper| replace:: :arxiv:`2608.XXXXX`
"""

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_theme_options = {'navigation_depth': 2}
html_copy_source = False
html_show_sourcelink = False
# html_static_path = ['_static']    # unused
