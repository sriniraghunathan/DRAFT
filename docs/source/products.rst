Data and products
=================

The repository ships both the input data that the pipeline reads, and the ILC
residual products and lensing noise curves for a range of experimental
configurations that have been pre-computed. Inputs are contained in ``data/``
and the released products can be accessed under ``products/``.

Input data
----------

**Fiducial cosmology.** ``data/output_planck_r_0.0_2015_cosmo_lensedCls.dat``
holds the lensed CMB power spectra used as the fiducial signal. The parameter
values, the step sizes for the numerical derivatives and the foreground
configuration are read from ``params.ini`` at the repository root.

**Galactic spectra.** ``data/galactic/`` holds the dust and synchrotron power
spectra as precomputed pySM3 outputs (based on `these simulations
<https://github.com/CMB-S4/s4mapbasedsims/tree/main/202102_design_tool_run>`_),
evaluated on the masked sky.
`PySM 3 <https://so-pysm-models.readthedocs.io/en/latest/>`_ itself is not run
by this repository.

**Extragalactic templates.** ``data/george_plot_bestfit_line.sav`` contains
the South Pole Telescope best-fit extragalactic foreground power spectra from
George et al. (2015), which are the amplitudes and frequency scalings applied in
:mod:`foregrounds`. Five further template files are shipped, but not used in
the current pipeline:
``data/dl_shaw_tsz_s10_153ghz_norm1_fake12000.txt`` and
``data/cl_tsz_148_batt.dat`` for the thermal Sunyaev-Zel'dovich effect,
``data/dl_ksz_CSFplusPATCHY_13sep2011_norm1_fake12000.txt`` for the kinematic
Sunyaev-Zel'dovich effect, ``data/dl_cib_1halo_norm1_12000.txt`` and
``data/dl_cib_2halo_norm1_12000.txt`` for the cosmic infrared background.

.. TODO: Remove these other files or incorporate them in the pipeline?

**Survey noise levels.**
``data/cmbs4_chile_opt_survey_patch_noise_levels.npy`` holds the per-patch
white-noise levels for the CMB-S4 Chile-only revised configuration, tabulated
in §2.4 of |paper|. All other configurations have their specifications coded
directly in :func:`exp_specs.get_exp_specs`.

**Figures.** ``data/planck_gal_fg_masks_with_cmbs4_footprint.png`` shows the
galactic masks with the CMB-S4 footprint applied, and
``data/s4_wide_specs_pbdr.png`` and ``data/s4_ultradeep_specs_pbdr.png`` show
the instrument and noise specifications of the conceptual design. The masks
themselves are currently not included.

.. TODO: include masks?
.. Earlier: "Planck Galactic masks with the CMB-S4 footprint applied"


Released products
-----------------

``products/`` contains three trees, produced at different times and with
different configurations. Each holds ILC residual power spectra as ``.npy``
files, with lensing-noise curves alongside them where available. Products from
the same experiment can appear in more than one tree, with the later generation
superseding the earlier one where they overlap.

CMB-S4
^^^^^^

.. rubric:: products/202310xx_PBDR_config/

The two-site conceptual design of §2.2 of |paper|, referred to in the
repository as the preliminary baseline design or PBDR. Contains
``s4deepv3r025_202310xx_pbdr_config`` for the CMB-S4 Ultra-deep survey at the
South Pole, ``s4wide_202310xx_pbdr_config`` for the CMB-S4 Wide survey in
Chile, and ``s4wide_achieved_performance_pbdr_202312xx`` for a variant of the
latter with noise depths scaled from achieved performance.

.. TODO: Last part correct above?

For ``s4wide_202310xx_pbdr_config``, the tree also provides a scan over
integration time, ``for1years`` through ``for10years``. In addition,
``lmax_12000/`` repeats the same configurations to a higher maximum multipole,
with the results being identical over the shared multipole range. ``s4wide/``
is the pre-PBDR baseline, which differs from ``s4wide_202310xx_pbdr_config``
only in the 27 and 39 GHz white-noise levels. ``202303xx/`` is an earlier
generation of the same configuration and is superseded by its parent directory.

.. TODO: Remove 202303xx/?

.. rubric:: products/s4_all_chile_config/

.. TODO: Remove the second subtree and simplify?

The Chile-only revised configuration of §2.4 of |paper|. While two subtrees
are currently included, ``report/lmax_6500/`` is the authoritative branch. It
scans the survey year from ``---year1.0`` to ``---year15.0``. ``lmax_6500/``
contains an earlier generation covering the same survey and patch combinations,
but with its noise in polarization incorrectly being the same as in
temperature. The ``20250530_wrong_1_over_f_definitions/`` subdirectory holds a
superseded set computed with different atmospheric-noise definitions (currently
only retained for reference and should not be used).

.. TODO: Remove 20250530_wrong_1_over_f_definitions/?

.. rubric:: products/20220726/

The earliest tree, with its ``s4wide`` and ``s4deepv3r025`` being superseded by
the PBDR tree. It also contains an ``lmax_10000/`` variant, lensing-noise
curves and the plotting scripts that produced ``ilc_residuals.png`` *[with the
latter being moved soon]*.

.. TODO: Perform this move

.. rubric:: Coverage of the paper's configurations

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - Configuration
     - Released products
   * - §2.1 Common specifications, including Planck
     - No standalone products; only used in forecasting
   * - §2.2 Two-site conceptual design
     - ``products/202310xx_PBDR_config/``
   * - §2.3 South-Pole-only alternative
     - none *[yet]*
   * - §2.4 Chile-only revised configuration
     - ``products/s4_all_chile_config/``
   * - §2.5 Cosmic-variance-limited survey
     - none (no ILC performed)

.. TODO: Add South-Pole-only?

Within the revised configuration, the three survey components of §2.4 map onto
the product names as follows: ``lat_wide---patch<N>+advanced_so_goal`` covers
the sky where the CMB-S4 Hybrid wide survey and the SO-like large-aperture
telescope overlap, ``lat_delensing---patch<N>`` is the CMB-S4 Hybrid Delensing
survey, which has no SO combination, and ``advanced_so_goal`` alone covers the
sky observed by the SO-like telescope but by neither CMB-S4 survey. Products
exist for patches 1 to 4, of which the paper uses patches 1 and 2.

.. TODO: check and correct

Three further survey variants are released: ``lat_wide_dc0`` corresponds to
calculations based on the CMB-S4 Data Challenge 0 (DC 0), ``lat_wide_phase2``
is from a contingency-driven expansion phase of the CMB-S4 project that was
considered during its design phase and the configuration of ``lat_roman`` is
related to the Nancy Grace Roman Space Telescope. The supported names are
listed in full in :func:`exp_specs.get_exp_specs`.

Simons Observatory
^^^^^^^^^^^^^^^^^^

``products/20220726/`` holds ILC residuals and lensing-noise curves for the
``sobaseline`` and ``sogoal`` configurations, computed over the same six
frequency bands as the CMB-S4 products.

``products/s4_all_chile_config/`` holds products for the upgraded
``advanced_so_baseline`` and ``advanced_so_goal`` configurations. The
``report/lmax_6500/`` subtree scans the survey year from ``---year1.0`` to
``---year15.0``.

South Pole Telescope
^^^^^^^^^^^^^^^^^^^^

``products/20220726/`` holds ILC residuals and lensing-noise curves for
``spt3g``. These combine three frequency bands and carry the ``nocorrnoise``
token, i.e. they were computed without correlated atmospheric noise between
bands. ``products/s4_all_chile_config/lmax_6500/`` additionally holds
``spt3g_plus_spt3g+_WG2``, which is an SPT-3G and SPT-3G+ combination.

Filename grammar
----------------

Product files are named

.. code-block:: text

   <expname>_ilc_galaxy<0|1>_<bands>_<spectra>[_galmask<N>][_AZ|_CU][_withtszxcib]_lmax<L>_for<M>years.npy

with the following tokens:

``<expname>``
   The experiment configuration, as accepted by ``-expname``. The grammar,
   including the ``---patch<N>``, ``---year<Y>`` and ``+advanced_so_*``
   suffixes, is documented in :func:`exp_specs.get_exp_specs`.

``galaxy0`` / ``galaxy1``
   Galactic foregrounds excluded or included.

``<bands>``
   Band centres in GHz, joined by hyphens, in the order used for the
   covariance and the weights.

``<spectra>``
   The calculated spectra, ``TT-EE`` for every released product.

``galmask<N>``
   Galactic mask, present only when ``galaxy1``. ``0``, ``1`` and ``2`` select
   the Planck GAL070, GAL080 and GAL090 masks intersected with the CMB-S4
   footprint.

``_AZ`` / ``_CU``
   Which set of galactic simulations was used. The driver writes ``_CU`` when
   the simulation folder path contains ``CUmilta`` and ``_AZ`` otherwise, i.e.
   ``_AZ`` is the default.

``_withtszxcib``
   The tSZ-CIB correlation was included. Absent by default.

``lmax<L>``
   Maximum multipole. Note that the array runs from :math:`\ell = 0` to
   :math:`\ell = L - 1`, so ``lmax6500`` gives 6500 entries ending at 6499.

``for<M>years``
   Integration time in years.

Two further tokens appear on some files: ``nocorrnoise`` for SPT configurations
computed without correlated atmospheric noise and ``nulled_<component>`` for
constrained ILC products in which a component was nulled.

Note: Survey duration is encoded in two different ways. In
``products/202310xx_PBDR_config/``, it is the ``for<T>years`` token in the
filename. In ``products/s4_all_chile_config/report/`` it is the ``---year<Y>``
token inside the experiment name, with the ``for<T>years`` token not reflecting
the depth, but the total survey duration from which the ``year<Y>`` noise
levels were derived by rescaling with :math:`\sqrt{T/Y}`.

File format
-----------

Each ``.npy`` file holds a single pickled dictionary with eleven keys.

``el``
   Multipole array, shape ``(lmax,)``, running from 0 to ``lmax - 1``. The array
   index equals the multipole.

``cl_residual``
   Dictionary keyed ``'TT'`` and ``'EE'``, each a 1-D array over ``el`` giving
   the residual power :math:`C_\ell^\mathrm{res}` after component separation.
   These are :math:`C_\ell`, not :math:`D_\ell`, in units of μK².

``weights``
   Dictionary keyed ``'TT'`` and ``'EE'``, each of shape ``(n_bands, lmax)``.
   The ILC weights are dimensionless.

``fg_res_dic``
   Dictionary keyed ``'TT'`` and ``'EE'``, each a dictionary over the
   individual components ``cib``, ``galdust``, ``galsync``, ``noise``,
   ``radio`` and ``tsz``, giving their separate contributions to the residual.

``beam_noise_dic``
   Dictionary keyed ``'T'`` and ``'P'``, mapping band frequency to
   ``[beam_fwhm_arcmin, white_noise_level]``.

``elknee_dic``
   Dictionary keyed ``'T'`` and ``'P'``, mapping band frequency to the
   atmospheric (1/f) noise parameters ``[elknee, alphaknee]``.

``fsky_val``
   Sky fraction of the mask used.

``which_gal_mask``
   Galactic mask index, matching the ``galmask<N>`` token.

``param_dict``
   The parameter dictionary read from ``params.ini``, including the foreground
   configuration and the paths of the galactic simulation files.

``cl_multiplier_dic``
   Component rescaling factors for a partial ILC. Empty in the released
   products, which use the standard minimum-variance form.

``freqcalib_fac``
   Per-band calibration factors. ``None`` in the released products.

Note that ``cl_residual`` is exactly zero below :math:`\ell = 11`, i.e. the
usable range is :math:`11 \le \ell <` ``lmax``.

Reading the products
--------------------

``products/202310xx_PBDR_config/read_ilc_residuals.py`` and
``read_lensing_noise.py`` load a single file and plot it, and
``make_ascii.py`` writes the plain-text exports. All three take the filename
from a variable at the top of the script rather than from the command line.
``products/20220726/`` contains its own plotting scripts, which glob the
current directory and must therefore be run from inside it.
*[A ``load_product()`` function in a new module, superseding these scripts and
applicable to all three trees, is to be added.]*

.. TODO: New module

Loading a product directly requires both ``allow_pickle`` and the ``latin1``
encoding, the latter a compatibility setting for pickles written under
Python 2:

.. code-block:: python

   import numpy as np

   res_dic = np.load(fname, allow_pickle=True, encoding='latin1').item()
   el = res_dic['el']
   cl_residual = res_dic['cl_residual']['TT']

In addition, some plain-text exports accompany a subset of the products under
``products/202310xx_PBDR_config/``. They cover ILC residuals for ``galaxy0``
configurations and some lensing-noise curves.
