Adding an experiment
====================

Experiment configurations are defined in :func:`exp_specs.get_exp_specs`, a
single dispatch on ``expname`` that returns the per-band beams, white-noise
levels and atmospheric-noise parameters. Adding a configuration means adding
one branch to that dispatch, with no other changes in the pipeline necessary.

Required information
--------------------

``specs_dic``
   The specifications, keyed by band centre in GHz. Each value is the
   seven-element list
   ``[beam_fwhm, delta_T, elknee_T, alphaknee_T, delta_P, elknee_P, alphaknee_P]``,
   with the beam full-width at half-maximum in arcmin and the white-noise
   levels in μK arcmin. Setting ``elknee`` to ``-1`` switches off the
   atmospheric (1/f) noise term for that band and spectrum.

``corr_noise``, ``rho``, ``corr_noise_bands``
   Whether atmospheric noise is correlated between bands, the correlation
   coefficient and for each band the list of bands it is correlated with. For
   an uncorrelated configuration, set ``corr_noise = 0``, ``rho = 0.`` and map
   every band to itself alone.

``Nred_dic`` is initialized to an empty dictionary before the dispatch and only
needs setting for configurations that model red noise explicitly, which at
present is the Simons Observatory family.

Worked example
--------------

The Planck branch is the simplest in the file and shows the shape of a
self-contained configuration:

.. code-block:: python

   elif expname.lower() == 'planck':
       #Planck instrumental noise, from Table 1 of Allison et al. 2015 (arXiv:1509.07471).
       no_pol_noise = 1e4  #μK-arcmin
       specs_dic = {
       #freq: [beam_arcmins, white_noise_T, elknee_T, alphaknee_T, whitenoise_P, elknee_P, alphaknee_P]
       30:  [33., 145., -1., 0., no_pol_noise, -1., 0.],
       44:  [23., 149., -1., 0., no_pol_noise, -1., 0.],
       70:  [14., 137., -1., 0., 450.,         -1., 0.],
       100: [10.,  65., -1., 0., 103.,         -1., 0.],
       143: [7.,   43., -1., 0., 81.,          -1., 0.],
       217: [5.,   66., -1., 0., 134.,         -1., 0.],
       353: [5.,  200., -1., 0., 406.,         -1., 0.],
       }

       corr_noise = 0
       rho = 0.
       corr_noise_bands = {}
       for nu in specs_dic:
           corr_noise_bands[nu] = [nu]

Two conventions are visible here. Atmospheric noise is switched off in every
band with ``elknee = -1``, since it is not relevant for a satellite, and the
two channels without quoted polarization sensitivity are given a deliberately
large ``delta_P`` of :math:`10^4` μK arcmin, so that the ILC down-weights them
in polarization rather than excluding them.

Adding your own
---------------

Insert an ``elif`` branch before the final ``else`` that raises for an
unrecognized name, following the pattern above, then run it:

.. code-block:: bash

   python3 get_ilc_residuals.py -expname mytelescope -interactive_mode 0

The white-noise levels are printed per band, so a first check is that the
printed ``Delta T`` and ``Delta P`` are what you intended after any
observing-time rescaling.

Things to watch for
-------------------

**Dispatch is by substring and order matters.** The families are tested in the
order given in :func:`exp_specs.get_exp_specs` and each test is a substring
match, so a name containing ``s4``, ``atlast``, ``advanced_so`` or ``spt``
anywhere is routed to that family before your branch is reached. Choose a name
that contains none of them unless you intend that routing.

**The observing-time token in the output filename is also gated on** ``s4``.
For a name that does not contain it, ``-total_obs_time`` rescales the noise
levels correctly but the filename records no ``for<M>years`` token, so two runs
of the same configuration with different observing times write to the same path
and the second silently overwrites the first. Until this is fixed, either
include ``s4`` in the name or move each result before rerunning.

.. TODO: address the above?

**Galactic foregrounds need the standard bands.** ``-include_gal 1`` requires
the frequency bands of the configuration to be 27, 39, 93, 145, 225 and 278 GHz
because the precomputed simulations exist only for those. Any other set raises
a ``ValueError`` naming the bands that are missing.

**Check that your branch is reachable.** The dispatch is long and several parts
of the code may be unreachable or are overwritten. After adding a branch,
confirm the printed noise levels match your input.

*[Longer term, the specifications are to be moved out of the dispatch into a data
file or per-experiment dictionaries, so that adding a configuration does not
require editing a long chain.]*

.. TODO
