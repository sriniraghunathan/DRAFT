DRAFT: Dark Radiation Anisotropy Flowdown Team forecasting tool
===============================================================

The **DRAFT** (Dark Radiation Anisotropy Flowdown Team) forecasting tool
provides an end-to-end pipeline for cosmic microwave background (CMB)
experiments, from experimental specifications to projected constraints on
cosmological parameters. It was developed within the Maps-to-Power-Spectra
analysis working group of the CMB-S4 collaboration and is described in §3.4 of
|paper|. This code is however not specific to CMB-S4, but can be applied
directly to other CMB survey designs, and the experiment library already covers
configurations for Planck, the South Pole Telescope, the Simons Observatory,
CMB-HD and AtLAST alongside the CMB-S4 surveys.

The pipeline targets damping-tail science, where precise measurements at
multipoles :math:`\ell \gtrsim 1000` carry the constraining power but
astrophysical foreground emission also becomes significant relative to the
primary signal and where gravitational lensing broadens the acoustic peaks.
Foreground modeling, component separation and delensing are therefore all needed
to extract the available information. A leading example is the effective number
of relativistic species, :math:`N_\mathrm{eff}`, whose sensitivity is driven by
the damping tail and the acoustic peaks.

Pipeline
--------

.. figure:: figures/pipeline.*
   :alt: Schematic overview of the DRAFT forecasting pipeline, from inputs
         through sky preparation, component separation and forecasting to
         outputs.
   :width: 100%

   Schematic overview of the DRAFT forecasting tool and the pipeline it
   implements, from the instrumental and astrophysical inputs to the projected
   sensitivity to cosmological parameters. (Reproduced from |paper|.)

**Inputs.** The fiducial cosmology together with the parameter step sizes for the
numerical derivatives, the Planck information entering through
inverse-variance combination of its noise, the galactic dust and synchrotron
models, taken from `map-based simulations
<https://github.com/CMB-S4/s4mapbasedsims/tree/main/202102_design_tool_run>`_
generated with `PySM 3 <https://so-pysm-models.readthedocs.io/en/latest/>`_,
the extragalactic templates for the thermal and kinematic Sunyaev-Zel'dovich
effects, the cosmic infrared background and radio galaxies, following a
parameterization based on South Pole Telescope measurements, and the
experimental specifications of the survey configurations under consideration.
See :doc:`modules` for the modules that supply each of these (and
Sections 2 and 3 of |paper|).

**Sky preparation.** The survey footprint is combined with a galactic mask to
remove the most contaminated regions of the sky, which component separation can
only partially clean (§3.1 of |paper|).

**Component separation.** The multi-frequency covariance
:math:`\mathbf{C}_\ell` of the noise and foreground spectra is assembled from
the auto- and cross-frequency spectra, and an internal linear combination
optimally combines the frequency bands to extract the required component while
down-weighting modes dominated by foregrounds or instrumental noise. The
standard minimum-variance, constrained and partial variants are supported
*[cross-ILC still to be ported]*, leaving the residual power spectra
:math:`C_\ell^\mathrm{res}` (§3.2 of |paper|).

.. TODO: port cross-ILC

**Forecasting.** The post-ILC noise spectra are passed to
`CLASS_delens <https://github.com/selimhotinli/class_delens>`_, which computes
the theoretical CMB spectra and performs the iterative delensing. The delensed
spectra and the lensing-reconstruction noise are then combined in a
Fisher-matrix forecast performed by
`FisherLens <https://github.com/ctrendafilova/FisherLens>`_ to obtain projected
constraints on the cosmological parameters. Both codes are driven from within
this repository, which treats the surveys and sky patches of a configuration as
independent experiments and sums their Fisher matrices (§3.3 of |paper|). The
same formalism can also be used to estimate the biases induced in those
parameters by residual foregrounds, which is implemented but still experimental.

**Outputs.** The ILC residuals, the lensing-reconstruction noise, the Fisher
matrix :math:`F_{\alpha\beta}` and the projected uncertainties on the
cosmological parameters, in particular :math:`\sigma(N_\mathrm{eff})`.

All three stages are implemented in this repository. The intermediate ILC
residual spectra are released for a range of survey configurations and can
therefore also be used as inputs to other forecasting frameworks, as described
in :doc:`products`.

Where to go next
----------------

:doc:`products` describes the input data and the released data products
together with their file format, which is where to start if you want the
released curves without running anything. :doc:`installation` covers the
requirements and the build, :doc:`quickstart` how to run the two drivers and
the arguments they take, :doc:`modules` the API organized by pipeline stage,
and :doc:`references` the literature underlying the models and methods.

.. toctree::
   :maxdepth: 2
   :caption: Getting started

   installation
   quickstart
   new_experiment

.. toctree::
   :maxdepth: 2
   :caption: API reference

   modules

.. toctree::
   :maxdepth: 1
   :caption: Additional information

   products
   references

Indices
-------

* :ref:`genindex` -- all documented functions and terms.
* :ref:`modindex` -- the module index.
* :ref:`search` -- full-text search.
