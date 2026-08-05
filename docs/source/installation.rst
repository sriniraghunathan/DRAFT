Installation
============

DRAFT is used by cloning the repository and running the driver from the
repository root. There is no installable package *yet* (see `Planned importable
API`_ below).

Requirements
------------

* Python 3.12+.
* numpy, scipy, matplotlib.

Optional:

* healpy, needed only by :func:`misc.healpix_rotate_coords`, which imports it
  lazily and raises ``ImportError`` if it is missing. The component
  separation pipeline does not use it.

The galactic foreground spectra are read from precomputed pySM3 outputs shipped
under ``data/``, so pySM3 itself is not required. The delensing and
Fisher-matrix stages are performed by two companion codes with their own
requirements, listed in :doc:`modules`.

Getting the code
----------------

.. code-block:: bash

   git clone https://github.com/sriniraghunathan/DRAFT.git
   cd DRAFT
   pip install numpy scipy matplotlib

Running it
----------

All commands are run from the repository root, since the modules are currently
imported by name rather than as an installed package:

.. code-block:: bash

   python3 get_ilc_residuals.py -expname s4wide_202310xx_pbdr_config \
       -include_gal 1 -which_gal_mask 2 -final_comp cmb

See :doc:`quickstart` for the available arguments and :doc:`products` for the
output format.

Verifying the checkout
----------------------

The driver's help output lists every argument and requires no data:

.. code-block:: bash

   python3 get_ilc_residuals.py --help

To check that the modules import, from the repository root:

.. code-block:: python

   import sys
   sys.path.insert(0, 'modules')
   import exp_specs, foregrounds, ilc, misc

Planned importable API
----------------------

A small importable API is planned, initially unstable, e.g. with names
potentially changing initially. The released product files are however
treated as stable. Until then, ``pip install .`` does not work.
``pyproject.toml`` is incomplete. (To be added.)
