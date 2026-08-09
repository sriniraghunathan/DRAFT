Quickstart
==========

The DRAFT forecasting tool is run from the repository root and has two
drivers. ``get_ilc_residuals.py`` computes the ILC residual power spectra for
one experiment configuration and writes them to a product file.
``get_fisher_forecasts.py`` reads one or more of those products and turns them
into a Fisher forecast. See :doc:`installation` for the requirements, which
differ between the two.

First ILC run
-------------

.. code-block:: bash

   python3 get_ilc_residuals.py -expname s4wide_202310xx_pbdr_config

This performs a minimum-variance CMB ILC for the CMB-S4 Wide survey of the
preliminary baseline/conceptual design, without galactic foregrounds, over the
default seven years of observation. The white-noise levels are printed for each
band, the product is written under ``results/`` and its path is printed at the
end.

``-expname`` is the only required argument; everything else has a default.

.. code-block:: bash

   python3 get_ilc_residuals.py --help

lists every argument together with its current default and is the reference for
the exhaustive list. The sections below cover when to depart from the defaults.

Common options
--------------

**Choosing an experiment.** ``-expname`` selects the frequency bands, beams,
white-noise levels and atmospheric noise. The supported names and the grammar
for the ``---patch<N>``, ``---year<Y>`` and ``+advanced_so_*`` suffixes are
documented in :func:`exp_specs.get_exp_specs`. The configurations, for which
results have already been released, are described in :doc:`products`.

**Galactic foregrounds.** ``-include_gal 1`` adds galactic dust and
synchrotron, read from the precomputed simulations under ``data/galactic/``,
and ``-which_gal_mask`` then selects the mask: ``0``, ``1`` and ``2`` for the
Planck GAL070, GAL080 and GAL090 masks intersected with the survey footprint.

.. code-block:: bash

   python3 get_ilc_residuals.py -expname s4wide_202310xx_pbdr_config \
       -include_gal 1 -which_gal_mask 2

Because the galactic simulations are only available for a fixed set of bands,
``-include_gal 1`` requires every band of the configuration to be one of 27,
39, 93, 145, 225 and 278 GHz. A configuration carrying any other band raises a
``ValueError`` naming the offending bands rather than failing later in the
computation.

**What the ILC recovers.** ``-final_comp`` chooses the component preserved with
unit response, one of ``cmb``, ``tsz``, ``y``, ``radio`` or ``cib``. The
general default for standard calculations and forecasts is ``cmb``. Adding
``-null_comp`` performs a constrained ILC that additionally nulls one or more
other components:

.. code-block:: bash

   python3 get_ilc_residuals.py -expname s4wide_202310xx_pbdr_config \
       -final_comp cmb -null_comp tsz

Several components can be nulled at once by listing them,
``-null_comp tsz cib``. Nulling is recorded in the output filename, so
constrained and standard runs of the same configuration do not overwrite one
another.

**Survey duration.** ``-total_obs_time`` sets the observing time in years and
the white-noise levels are rescaled from the baseline of the configuration
accordingly. Experiment names carrying a ``---year<Y>`` suffix already encode
an observing time, so combining the two raises a ``ValueError`` rather than
applying the scaling twice.

**Per-band noise scaling.** ``-noise_scalings_for_bands`` multiplies the
white-noise level of each band by a given factor, one value per band in band
order, which is useful for exploring sensitivity to individual channels.

**Diagnostics and bookkeeping.** ``-debug 1`` prints the intermediate
dictionary keys and progress, and draws the diagnostic plots of the beams and
the noise spectra. A default run draws nothing. ``-interactive_mode 0`` then
forces the non-interactive ``Agg`` backend rather than showing those figures; on
its own, without ``-debug``, it has no effect. ``-save_fg_res_and_weights 0``
omits the per-component residuals and the ILC weights from the product file.
``-paramfile`` points at an alternative parameter file in place of
``params.ini``.

From Python
-----------

The same calculation is available to other scripts through
:func:`get_ilc_residuals.run_ilc`, which takes the command line options as
keyword arguments and returns the product dictionary together with the path it
would be written to. The driver puts ``modules/`` on ``sys.path`` itself, so no
further setup is needed when working from the repository root:

.. code-block:: python

   import get_ilc_residuals

   opdic, opfname = get_ilc_residuals.run_ilc(
       's4wide_202310xx_pbdr_config', include_gal=1, which_gal_mask=2, save=False)

   el = opdic['el']
   cl_residual = opdic['cl_residual']['TT']

Four options are available only this way, with no command line equivalent:
``which_spec_arr`` to solve for spectra other than ``TT`` and ``EE``,
``cl_multiplier_dic`` to rescale a component's spectrum and so perform a
partial ILC, ``freqcalib_fac`` to apply a per-band gain, and ``save`` to return
the result without writing a file.

First forecast
--------------

Forecasting needs the two companion codes built, as described in
:doc:`installation`. It reads ILC products rather than recomputing them, so the
files above or the released ones under ``products/`` are its input:

.. code-block:: bash

   python3 get_fisher_forecasts.py -ilc_fname <wide product> \
       -fsky 0.62 -label s4_wide

The product paths are long, so refer to :doc:`products` for their grammar
and :mod:`get_fisher_forecasts` for a longer version of this command written
out.

Several products can be given at once. They are treated as independent
experiments whose Fisher matrices are summed, which is how the surveys and sky
patches of a configuration are combined into one forecast:

.. code-block:: bash

   python3 get_fisher_forecasts.py \
       -ilc_fname <wide product> <ultradeep product> \
       -fsky 0.59 0.03 -label s4_conceptual -write_cache 1

``-fsky`` takes one value per ``-ilc_fname``, either a number or ``auto`` to use
the ``fsky_val`` recorded in that product. A product built with
``-include_gal 0`` records none, so a number is required for it. ``-label``
names the output folder and the files within it, and is required once more than
one product is given.

Only what does not bear on a number is given on the command line. The
cosmology, the step sizes, the multipole ranges, the delensing mode, the priors
and the option flags are read from ``params_fisher.ini``, so a forecast is
reproducible from a file rather than from a shell history.

**Caching.** Each configuration costs one iterative CLASS_delens run plus two
CLASS_delens runs per varied parameter. ``-write_cache 1`` records the computed
spectra and derivatives, and ``-cache`` reads them back in place of
running CLASS, one value per ``-ilc_fname`` with ``none`` to run CLASS for that
one. A rerun at a different sky fraction, with different priors or with a
parameter held fixed then costs seconds and no CLASS run at all.

**Other options.** ``-write_lensing_noise 1`` writes the
lensing-reconstruction noise curves of each configuration. ``-save 0`` returns
the forecast without writing it, ``-opfname`` overrides the output path built
from ``-label``, ``-overwrite 1`` replaces existing output files rather than
refusing to, ``-verbose 0`` silences the per-step reporting and ``-paramfile``
points at an alternative parameter file in place of ``params_fisher.ini``.

The same calculation is available through
:func:`get_fisher_forecasts.run_forecast`, which returns the forecast together
with the path it would be written to:

.. code-block:: python

   import get_fisher_forecasts

   product, opfname = get_fisher_forecasts.run_forecast(
       ['<wide product>', '<ultradeep product>'], fsky=[0.59, 0.03],
       label='s4_conceptual', save=False)

Three options are available only this way: ``settings`` to pass the parameters
directly in place of reading a file, ``fix_params`` to hold parameters fixed
after the total Fisher matrix is formed and ``paths`` to point at a FisherLens
or CLASS_delens checkout elsewhere.

Reading the results
-------------------

ILC products are written under ``results/`` in a directory tree that records
the experiment, the mask and the spectra combined. The filename grammar and the
contents of the product dictionary are documented in :doc:`products`.

Everything a forecast produces is written together under
``results/forecasts/<label>/``, because a forecast combining several
configurations belongs to no single ILC product. The run above leaves

``<label>_forecast.npy``
   The forecast itself.

``<ilc product name>_cache.npy``
   One per configuration, written with ``-write_cache 1``.

``<ilc product name>_lensing_noise.npy``
   One per configuration, written with ``-write_lensing_noise 1``.

Each is read back by a function in :mod:`get_fisher_forecasts` rather than by
a script:

.. code-block:: python

   import get_fisher_forecasts as gff

   opdir = 'results/forecasts/s4_conceptual'

   ilc      = gff.load_ilc_product('<ilc product>')
   forecast = gff.read_forecast_product(opdir + '/s4_conceptual_forecast.npy')
   cache    = gff.read_forecast_cache('<cache file>')
   curves   = gff.read_lensing_noise_curves('<lensing-noise file>')

   gff.report_forecast(forecast)

:func:`get_fisher_forecasts.read_lensing_noise_curves` also reads the older
lensing-noise files released under ``products/``, converting them to the
current layout and warning that it has done so.

A forecast product records the Fisher matrix, its inverse and its condition
number under ``fisher``, ``covariance`` and ``condition``, each keyed by
spectrum type, ``unlensed``, ``lensed`` and ``delensed``. ``sigmas`` is keyed
by spectrum type and then by parameter name. ``params`` gives the surviving
parameter order, with ``params_supplied``, ``params_fixed``, ``priors``,
``spectrum_types``, ``survey_labels`` and ``satellite_label`` recording how the
total was formed, ``configurations`` holding one entry per input with its
``fsky`` and the files it produced, and ``settings`` the parameters the run
used. :func:`get_fisher_forecasts.report_forecast` prints the projected
uncertainties of one spectrum type, ``delensed`` unless another is given.
