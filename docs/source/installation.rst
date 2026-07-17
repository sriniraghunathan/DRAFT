Installation
============

Requirements
------------

* Python >= 3.6
* numpy, scipy, matplotlib
* pyparsing >= 2.0.2

Optional (depending on the analysis):

* `CAMB <https://camb.readthedocs.io/>`_ -- CMB power spectra and
  derivatives for Fisher forecasting.
* `pySM3 <https://pysm3.readthedocs.io/>`_ -- Galactic dust and
  synchrotron simulations.
* healpy -- HEALPix map handling.

From source
-----------

Clone the repository and install with pip:

.. code-block:: bash

   git clone https://github.com/sriniraghunathan/DRAFT.git
   cd DRAFT
   pip install .

For development (editable install):

.. code-block:: bash

   pip install -e .

Verifying the installation
--------------------------

.. code-block:: python

   import DRAFT
   print(DRAFT.__file__)
