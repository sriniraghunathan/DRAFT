DRAFT: Dark Radiation Anisotropy Flowdown Team tool
====================================================

**DRAFT** is the CMB-S4 Dark Radiation Anisotropy Flowdown Team tool for
component separation and cosmological forecasting.

The code:

* Optimally combines data from different frequency bands using noise and
  foreground signals via internal linear combination (ILC), supporting
  standard, constrained, and partial ILC.
* Combines delensed CMB spectra and lensing spectra to forecast
  cosmological parameter constraints using the Fisher formalism.
* Estimates biases in cosmological parameters due to residual foregrounds,
  also using the Fisher formalism.

Foreground modelling
--------------------

* **Extragalactic foregrounds:** Radio, CIB, tSZ, and kSZ power spectra
  from SPT measurements (George et al. 2015, arXiv:1408.3161;
  Reichardt et al. 2020, arXiv:2002.06197). Default polarisation
  fractions: CIB = 2%, Radio = 3%, tSZ/kSZ = 0 (configurable).
* **Galactic foregrounds:** Dust and synchrotron power spectra obtained
  from pySM3 simulations.

.. toctree::
   :maxdepth: 2
   :caption: Getting Started

   installation
   quickstart

.. toctree::
   :maxdepth: 2
   :caption: API Reference

   modules

.. toctree::
   :maxdepth: 1
   :caption: Additional Information

   products
   references

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
