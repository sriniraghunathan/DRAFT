Quickstart
==========

DRAFT is run from the repository root. Each run computes the ILC residual power
spectra for one experiment configuration and writes them to a product file. See
:doc:`installation` for the requirements.

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
``-include_gal 1`` requires a configuration whose bands are 27, 39, 93, 145,
225 and 278 GHz. Any other configuration raises a ``ValueError`` naming the
offending bands rather than failing later in the computation.

**What the ILC recovers.** ``-final_comp`` chooses the component preserved with
unit response, one of ``cmb``, ``tsz``, ``y``, ``radio`` or ``cib``. The
general default for standard calculations and forecasts is ``cmb''. Adding
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

**Diagnostics and bookkeeping.** ``-interactive_mode 0`` suppresses the
diagnostic plots of the beam and noise spectra, which is what you want when
running non-interactively. ``-save_fg_res_and_weights 0`` omits the
per-component residuals and the ILC weights from the product file.
``-debug 1`` prints intermediate dictionary keys and progress. ``-paramfile``
points at an alternative parameter file in place of ``params.ini``.

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

Reading the results
-------------------

Products are written under ``results/`` in a directory tree that records the
experiment, the mask and the spectra combined. The filename grammar and the
contents of the product dictionary are documented in :doc:`products`. The
released products come with reader scripts,
``products/202310xx_PBDR_config/read_ilc_residuals.py`` for the ILC residuals
and ``read_lensing_noise.py`` for the lensing-noise curves.
