API reference
=============

**DRAFT** provides an end-to-end forecasting pipeline from experimental
specifications to projected parameter constraints, covering foreground modeling,
sky masking, ILC-based component separation, delensing and the Fisher-matrix
calculation. The pipeline as a whole is described in §3.4 of |paper|, with the
underlying formalism in §3.1 to §3.3.

The component separation stage is implemented in this repository. The delensing
and Fisher-matrix stages are performed by two companion codes and are not yet
integrated here.

.. note::

   The pages below are currently maintained by hand, with each module described
   by its own module docstring. Do not regenerate them with ``sphinx-apidoc``.

Inputs
------

:doc:`exp_specs`
   Experimental specifications for each survey configuration: frequency bands,
   beam widths :math:`\theta_b`, white-noise levels :math:`\Delta_{T,P}` and the
   atmospheric (:math:`1/f`) parameters :math:`\ell_\mathrm{knee}` and
   :math:`\alpha_\mathrm{knee}`. Also covers the Planck configurations whose
   inverse-variance-combined noise enters the forecast (§2.1 of |paper|).

:doc:`foregrounds`
   Galactic dust and synchrotron power spectra obtained from pySM3 simulations,
   and the extragalactic CIB, radio, tSZ and kSZ templates (§3.1 of |paper|).

The fiducial cosmology, the parameter step sizes for the numerical derivatives
and the foreground configuration are read from ``params.ini``.

.. toctree::
   :hidden:
   :maxdepth: 2

   exp_specs
   foregrounds

Component separation
--------------------

:doc:`get_ilc_residuals`
   Driver for this stage. Assembles the per-band beams and noise levels for a
   named configuration, applies the galactic mask, builds the covariance, runs
   the ILC and writes the residual power spectra :math:`C_\ell^\mathrm{res}`
   to a product file.

:doc:`ilc`
   The multi-frequency covariance :math:`\mathbf{C}_\ell` of the noise and
   foreground spectra, the ILC weights and the residual power. Supports the
   minimum-variance, constrained and partial variants (§3.2 of |paper|).

:doc:`misc`
   Supporting calculations: beam transfer functions :math:`B_\ell`, noise power
   spectra :math:`N_\ell` from the experimental specifications, parameter file
   I/O and map-domain utilities.

The data location for the galactic foregrounds and the values of the template
parameters for the extragalactic foregrounds are read from ``params.ini``.

.. toctree::
   :hidden:
   :maxdepth: 2

   get_ilc_residuals
   ilc
   misc

Forecasting
-----------

Not yet part of this repository (to be added):

* Theoretical CMB spectra and iterative delensing, performed by
  `CLASS_delens <https://github.com/selimhotinli/class_delens>`_.
* Fisher-matrix forecast, performed by
  `FisherLens <https://github.com/ctrendafilova/FisherLens>`_.

Both consume the ILC residuals produced by the component separation stage
(§3.3 of |paper|).

Utilities
---------

:doc:`flatsky`
   Flat-sky simulation and power-spectrum estimation, used for small patches in
   the flat-sky approximation. Not part of the forecasting pipeline.

.. toctree::
   :hidden:
   :maxdepth: 2

   flatsky
