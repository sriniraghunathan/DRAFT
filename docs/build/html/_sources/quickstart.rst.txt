Quickstart
==========

Computing ILC residuals
-----------------------

The main entry point computes ILC residual curves for a given experiment
configuration. For example, for the S4-Wide (Chilean LAT) PBDR
configuration, including Galactic foregrounds:

.. code-block:: bash

   python3 get_ilc_residuals.py \
       -expname s4wide_202310xx_pbdr_config \
       -include_gal 1 \
       -which_gal_mask 2 \
       -final_comp cmb

Key arguments:

``-expname``
    Name of the experiment configuration (frequency bands, beams, and
    noise levels), e.g. ``s4wide_202310xx_pbdr_config``.

``-include_gal``
    Whether to include Galactic foregrounds (``1``) or not (``0``).

``-which_gal_mask``
    Index of the Galactic foreground mask to use.

``-final_comp``
    The component to recover with the ILC, e.g. ``cmb``.

Reading the results
-------------------

* Use ``read_ilc_residuals.py`` to read the ILC residual curves.
* Use ``read_lensing_noise.py`` to read the lensing noise curves.

Supported experiment configurations include:

* **S4-Wide (Chilean LAT)** -- latest PBDR noise specs (Oct 2023).
* **S4-Ultra deep (South Pole LAT)** -- latest PBDR noise specs (Oct 2023).
* **S4-Wide (Chilean LAT) achieved performance** -- best achievable
  performance.
