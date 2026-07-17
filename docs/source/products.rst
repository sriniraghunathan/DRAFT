Data Products
=============

The ``products/`` directory of the repository contains pre-computed
results, including:

* **ILC residual curves** for the supported experiment configurations:

  * S4-Wide (Chilean LAT) -- PBDR noise specs (Oct 2023).
  * S4-Ultra deep (South Pole LAT) -- PBDR noise specs (Oct 2023).
  * S4-Wide (Chilean LAT) -- achieved (best possible) performance.

* **Lensing noise curves** for the above configurations.

* **Galactic foreground masks** (Planck Galactic masks with the CMB-S4
  footprint applied).

Use the ``read_ilc_residuals.py`` and ``read_lensing_noise.py`` scripts
to load these products.
