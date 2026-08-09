API reference
=============

The **DRAFT** forecasting tool provides an end-to-end pipeline from
experimental specifications to projected parameter constraints, covering
foreground modeling, sky masking, ILC-based component separation, delensing and
the Fisher-matrix calculation. The pipeline as a whole is described in §3.4 of
|paper|, with the underlying formalism in §3.1 to §3.3.

All these pipeline stages are implemented in this repository. The delensing and
the Fisher-matrix calculation are carried out by two companion codes, which
DRAFT drives rather than reimplements. See :doc:`installation` for how to build
them.

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

:doc:`get_fisher_forecasts`
   Driver for this stage. Reads one or more ILC products, treats each survey
   and sky patch as an independent experiment whose Fisher matrices are summed,
   drives the delensing and the Fisher-matrix calculation, and writes the
   forecast together with the projected parameter uncertainties.

:doc:`delensing`
   Iterative delensing. Combines the residual spectra with the effective Planck
   noise and hands the result to `CLASS_delens
   <https://github.com/selimhotinli/class_delens>`_, which reconstructs the
   lensing potential and returns the delensed spectra together with the
   lensing-reconstruction noise :math:`N_\ell^{dd}`.

:doc:`fisher`
   The Fisher matrix :math:`F_{\alpha\beta}` for one survey configuration and
   the separate large-scale contribution from a satellite, reduced to the
   projected uncertainties :math:`\sigma(p_\alpha)`.

Both companion codes are reached through `FisherLens
<https://github.com/ctrendafilova/FisherLens>`_, a git submodule of this
repository which brings CLASS_delens as a submodule of its own
(§3.3 of |paper|).

The fiducial cosmology, the parameter step sizes, the multipole ranges, the
priors and the option flags are read from ``params_fisher.ini``.

.. toctree::
   :hidden:
   :maxdepth: 2

   get_fisher_forecasts
   delensing
   fisher

Utilities
---------

:doc:`flatsky`
   Flat-sky simulation and power-spectrum estimation, used for small patches in
   the flat-sky approximation. Not part of the forecasting pipeline.

.. toctree::
   :hidden:
   :maxdepth: 2

   flatsky
