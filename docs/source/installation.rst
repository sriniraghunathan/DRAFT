Installation
============

The DRAFT forecasting tool is used by cloning the repository and running the
drivers from the repository root. There is no installable package *yet* (see
`Planned importable API`_ below).

The component separation stage needs only a few standard Python packages. The
forecasting stage additionally needs the two companion codes, which are built
from source. A reader who only wants ILC residuals can stop after
`Component separation`_ or skip the build entirely and use the released
products described in :doc:`products`.

The Python versions below are what the dependencies require. DRAFT itself is
tested on Python 3.12 and 3.13.

Getting the code
----------------

.. code-block:: bash

   git clone https://github.com/sriniraghunathan/DRAFT.git
   cd DRAFT

Component separation
--------------------

* Python 3.9+.
* numpy, scipy.

.. code-block:: bash

   pip install numpy scipy

Optional:

* matplotlib, needed only by the diagnostic plotting helpers
  :func:`get_ilc_residuals.plot_beams` and
  :func:`get_ilc_residuals.plot_noise_spectra`, which import it lazily and are
  reached only when ``-debug`` is set. A default run draws nothing and never
  imports it.
* healpy, needed only by :func:`misc.healpix_rotate_coords`, which imports it
  lazily and raises ``ImportError`` if it is missing. The component
  separation pipeline does not use it.

The galactic foreground spectra are read from precomputed pySM3 outputs shipped
under ``data/``, so pySM3 itself is not required.

Forecasting
-----------

* Python 3.10+, which is the floor set by CAMB.
* camb.
* A C compiler and ``make``, to build CLASS_delens.

`FisherLens <https://github.com/ctrendafilova/FisherLens>`_ is a git submodule
of this repository and brings
`CLASS_delens <https://github.com/selimhotinli/class_delens>`_ as a submodule
of its own, so the checkout has to be recursive:

.. code-block:: bash

   pip install camb
   git submodule update --init --recursive
   cd FisherLens/CLASS_delens
   make class

Use the ``class`` target rather than ``make`` or ``make all``, which also build
the ``classy`` Python wrapper and need Cython.

CAMB is required even though nothing on the CLASS code path, and therefore
nothing in DRAFT, calls it since FisherLens imports it at module scope.

Running it
----------

All commands are run from the repository root since the modules are currently
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

Once the companion codes are built, :func:`delensing.check_setup` reports
whether FisherLens imports, whether CLASS has been compiled and whether the
scratch directory is writable. It collects every problem before raising, so a
fresh checkout reports all of them at once:

.. code-block:: python

   import sys
   sys.path.insert(0, 'modules')
   import delensing
   delensing.check_setup()

Planned importable API
----------------------

A small importable API is planned, initially unstable, e.g. with names
potentially changing initially. The released product files are however
treated as stable. Until then, ``pip install .`` does not work.
``pyproject.toml`` is incomplete. *[To be added.]*

.. TODO: Complete and add.
