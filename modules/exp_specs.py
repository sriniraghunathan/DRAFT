"""
Instrument specifications for the experiment configurations supported by get_ilc_residuals.py.

The module exposes a single function, ``get_exp_specs``, which maps an experiment name onto the per-band beams, white-noise levels and atmospheric-noise parameters from which the noise power spectra are built in ``misc.get_nl``.

Seven experiment families are covered:

* CMB-S4: the CMB-S4 Wide and CMB-S4 Ultra-deep variants, the preliminary baseline/conceptual design, and the Chile-only revised configurations, which may be combined with an SO-like LAT,
* CMB-HD,
* AtLAST,
* Advanced Simons Observatory,
* Simons Observatory, baseline and goal,
* South Pole Telescope: SPT-SZ, SPTpol and SPT-3G,
* Planck.

``get_exp_specs`` is called by ``get_ilc_residuals.py``.
``data_folder`` points at the repository's ``data`` directory and is used to load the per-patch survey noise levels needed by the Chile-only revised configurations.
"""

import os
import warnings

import numpy as np


# Constants and module state

#data directory, resolved relative to the repo root rather than the current working directory
data_folder = os.path.join( os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'data' )


# Experiment specifications

def get_exp_specs(expname, corr_noise_for_spt=1, remove_atm=0):
    r"""
    Beams, noise levels and atmospheric-noise parameters for one experiment configuration.

    Each frequency band is described by a beam full-width at half-maximum :math:`\theta_b` in arcminutes, a white-noise level :math:`\Delta_X` in :math:`\mathrm{\mu K}\,\mathrm{arcmin}` for temperature and polarization (:math:`X = T, P`), and the low-frequency (1/f)/atmospheric parameters :math:`\ell_\mathrm{knee}` and :math:`\alpha_\mathrm{knee}` entering

    .. math::

        N_\ell^X = \Delta_X^2 B_\ell^{-2} \left[ 1 + \left( \frac{\ell_\mathrm{knee}}{\ell} \right)^{\!\alpha_\mathrm{knee}} \right] ,

    as evaluated by :func:`misc.get_nl`. An :math:`\ell_\mathrm{knee}` of ``-1`` disables the atmospheric term for that band.

    Parameters
    ----------
    expname : str
        Experiment configuration. The recognized configurations, grouped by experiment:

        * **CMB-S4**, selected by ``s4`` appearing anywhere in the name:

          - ``s4wide``, and the detector-scaled variants ``s4wide_scaled_sobaseline``, ``s4wide_scaled_aso``, ``s4wide_single_chlat`` and ``s4wide_single_chlat_plus_aso``
          - ``s4wide_scaled_aso_plus_fulls4scaledsobaseline`` and ``s4wide_single_chlat_plus_aso_plus_fulls4scaledsobaseline``, which add a full SO-baseline-scaled CMB-S4. Only these two names accept that suffix
          - ``s4wide_202310xx_pbdr_config``, ``s4deepv3r025_202310xx_pbdr_config``, ``s4wide_acheived_performance_pbdr_202312xx``, ``s4wide_chlat_el40``
          - ``s4deep``, ``s4deepv3r025``, ``s4deepv3r025_tma``
          - the Chile-only revised configurations, ``s4_all_chile_config_<survey>---patch<N>[---year<Y>]``, optionally combined with an SO-like LAT by appending ``+advanced_so_baseline`` or ``+advanced_so_goal``, itself optionally suffixed ``---year<Z>``.
            ``<survey>`` is one of ``lat_wide``, ``lat_wide_phase2``, ``lat_wide_dc0``, ``lat_roman`` or ``lat_delensing``, and ``<N>`` selects a sky patch, both looked up in ``data/cmbs4_chile_opt_survey_patch_noise_levels.npy``.

        * **CMB-HD**, selected by ``cmbhd`` or ``cmb-hd`` appearing anywhere in the name.
        * **AtLAST**, selected by ``atlast``: ``atlast`` or ``atlast_dummy``.
        * **Advanced Simons Observatory**, selected by ``advanced_so``: ``advanced_so_baseline`` or ``advanced_so_goal``, optionally suffixed ``---year<Y>``.
        * **Simons Observatory**: ``sobaseline`` or ``sogoal``, case insensitive. (``ccat_prime_so`` reaches this family but has no noise levels defined and raises.)
        * **South Pole Telescope**, selected by ``spt``: ``sptsz``, ``sptpol``, ``sptpol_summer``, ``sptpolsummer``, ``sptpolultradeep``, ``sptpolultradeepplus3g``, ``sptpolplusultradeep``, ``sptpolplusultradeepplus3g``, ``sptpolplusultradeepplus3gfull``, ``spt3g``, ``spt3gsummer``, ``spt3g_summer``, ``spt3g_y12``, ``spt3g_ZP``, ``spt3g_TC``, ``spt3g_WG2``, ``spt3g_plus_spt3g+_WG2``.
        * **Planck**: ``planck``, case insensitive.

        The tests are applied in the order above and are substring matches, so a name containing ``s4``, ``atlast``, ``advanced_so`` or ``spt`` anywhere is routed to that family.
        CMB-S4 and CMB-HD share a single dispatch test, so a name matching both is given the CMB-HD specifications.
    corr_noise_for_spt : int, optional
        Whether the SPT bands carry correlated atmospheric noise. Ignored by every other family, whose correlated-noise setting is fixed in the source. Default 1.
    remove_atm : int, optional
        If non-zero, switch off the atmospheric term by setting :math:`\ell_\mathrm{knee}` to ``-1`` and :math:`\alpha_\mathrm{knee}` to zero. Honored only by the Simons Observatory and South Pole Telescope families. Default 0.

    Returns
    -------
    specs_dic : dict
        Per-band specifications, keyed by frequency band centre in GHz.
        Each value is the seven-element list ``[beam_fwhm, delta_T, elknee_T, alphaknee_T, delta_P, elknee_P, alphaknee_P]``, with the beam in arcminutes and both noise levels in :math:`\mathrm{\mu K}\,\mathrm{arcmin}`.
    corr_noise_bands : dict
        For each band, the list of bands its atmospheric noise is correlated with.
        A band mapped to itself alone has no correlated partner, so the corresponding cross-band noise is zero.
    rho : float
        Correlation coefficient between the atmospheric noise of two correlated bands, passed to :func:`misc.get_nl`.
    corr_noise : int
        Whether correlated atmospheric noise is modelled for this configuration.
    Nred_dic : dict
        Red-noise amplitudes in :math:`\mathrm{\mu K}^2\,\mathrm{s}`, keyed by band, each a ``[Nred_T, Nred_P]`` pair in which ``-1`` disables the treatment.
        Empty except for the Simons Observatory family.

    Raises
    ------
    ValueError
        If ``expname`` is not a string, if it matches no family or matches a family but none of its configurations, if the Chile-only or Advanced SO name does not follow the ``---patch<N>`` / ``---year<Y>`` grammar, or if a parsed observing time is not positive.

    Warns
    -----
    UserWarning
        If ``remove_atm`` is set for a CMB-S4, CMB-HD, AtLAST or Advanced SO configuration, none of which implement it. The atmospheric term is retained.

    Notes
    -----
    Polarization noise levels are :math:`\Delta_P = \sqrt{2} \Delta_T` unless the configuration quotes them separately.
    Bands with no usable data in a given configuration are given a large :math:`\Delta_X` rather than being omitted, so that the ILC weights them to zero.

    Observing times are handled by rescaling the reference white-noise level of a configuration from its reference duration :math:`t_\mathrm{ref}` to the requested :math:`t`:

    .. math::

        \Delta_X \rightarrow \Delta_X \sqrt{ t_\mathrm{ref} / t } .

    Configurations that combine two instruments do so by inverse-variance weighting of the white-noise levels band by band:

    .. math::

        \Delta_X = \left( \Delta_{X,1}^{-2} + \Delta_{X,2}^{-2} \right)^{-1/2} .

    The returned dictionaries are built afresh on every call, so callers may mutate them without affecting later calls.
    ``rho`` is returned regardless of ``corr_noise`` and is unused when it is zero, in which case its value differs between families and should not be relied on.
    """

    if not isinstance(expname, str):
        raise ValueError('expname must be a string; got %s' % (type(expname).__name__))

    Nred_dic = {}

    if expname.find('s4') > -1 or expname.find('cmbhd') > -1 or expname.find('cmb-hd') > -1:
        if remove_atm:
            warnings.warn('remove_atm is not implemented for S4/CMB-HD configurations; the atmospheric term is retained', stacklevel=2)
        #if expname == 's4wide' or expname.find('s4wide_scaled_sobaseline') > -1 or expname.find('s4wide_scaled_aso') > -1 or expname.find('s4wide_single_chlat') > -1 or expname.find('s4wide_single_chlat_plus_aso') > -1 or expname.find('s4wide_plus') > -1:
        if expname == 's4wide' or expname.find('s4wide_scaled') > -1 or expname.find('s4wide_single') > -1:
            #20230517 - changing this to match PBDR
            '''
            specs_dic = {
            #freq: [beam_arcmins, white_noise_T, elknee_T, alphaknee_T, whitenoise_P, elknee_P, alphaknee_P]
            27: [7.4, 21.34, 415., 3.5, 30.23, 700, 1.4],
            39: [5.1, 11.67, 391., 3.5, 16.53, 700, 1.4],
            93: [2.2, 1.89, 1932., 3.5, 2.68, 700, 1.4],
            145: [1.4, 2.09, 3917., 3.5, 2.96, 700, 1.4],
            225: [1.0, 6.90, 6740., 3.5, 9.78, 700, 1.4],
            278: [0.9, 16.88, 6792., 3.5, 23.93, 700, 1.4],
            }
            '''
            specs_dic = {
            #freq: [beam_arcmins, white_noise_T, elknee_T, alphaknee_T, whitenoise_P, elknee_P, alphaknee_P]
            27: [7.4, 21.8, 415., 3.5, 30.23, 700, 1.4],
            39: [5.1, 12.4, 391., 3.5, 16.53, 700, 1.4],
            93: [2.2, 2.0, 1932., 3.5, 2.68, 700, 1.4],
            145: [1.4, 2.0, 3917., 3.5, 2.96, 700, 1.4],
            225: [1.0, 6.9, 6740., 3.5, 9.78, 700, 1.4],
            278: [0.9, 16.7, 6792., 3.5, 23.93, 700, 1.4],
            }

            '''
            if remove_atm:
                specs_dic = {
                #freq: [beam_arcmins, white_noise_T, elknee_T, alphaknee_T, whitenoise_P, elknee_P, alphaknee_P]
                27: [7.4, 21.34, 415., 0., 30.23, 700, 0.],
                39: [5.1, 11.67, 391., 0., 16.53, 700, 0.],
                93: [2.2, 1.89, 1932., 0., 2.68, 700, 0.],
                145: [1.4, 2.09, 3917., 0., 2.96, 700, 0.],
                225: [1.0, 6.90, 6740., 0., 9.78, 700, 0.],
                278: [0.9, 16.88, 6792., 0., 23.93, 700, 0.],
                }
            '''

            '''
            #DSR specs
            specs_dic = {
            #freq: [beam_arcmins, white_noise_T, elknee_T, alphaknee_T, whitenoise_P, elknee_P, alphaknee_P]
            #20: [10.0, None, None, None, None, None, None],
            27: [7.4, 21.8, 471., 3.5, 30.8, 700, 1.4],
            39: [5.1, 12.4, 428., 3.5, 17.6, 700, 1.4],
            93: [2.2, 2.0, 2154., 3.5, 2.9, 700, 1.4],
            145: [1.4, 2.0, 4364., 3.5, 2.8, 700, 1.4],
            225: [1.0, 6.9, 7334., 3.5, 9.8, 700, 1.4],
            278: [0.9, 16.7, 7308., 3.5, 23.6, 700, 1.4],
            #225: [1.0, 100., 7334., 3.5, 9.8, 700, 1.4],
            #278: [0.9, 100., 7308., 3.5, 23.6, 700, 1.4],
            }

            if remove_atm:
                specs_dic = {
                #freq: [beam_arcmins, white_noise_T, elknee_T, alphaknee_T, whitenoise_P, elknee_P, alphaknee_P]
                #20: [10.0, None, None, None, None, None, None],
                27: [7.4, 21.8, 471., 0., 30.8, 700, 0.],
                39: [5.1, 12.4, 428., 0., 17.6, 700, 0.],
                93: [2.2, 2.0, 2154., 0., 2.9, 700, 0.],
                145: [1.4, 2.0, 4364., 0., 2.8, 700, 0.],
                225: [1.0, 6.9, 7334., 0., 9.8, 700, 0.],
                278: [0.9, 16.7, 7308., 0., 23.6, 700, 0.],
                }
            '''

            #20220222 - modify S4 noise levels based on S4/SO detector scalings
            if expname == 's4wide_scaled_sobaseline':
                scaling_factors = np.sqrt(2.) * np.sqrt( [3.1, 3.1, 2.6*2., 2.5*2., 2.0*2., 1.8*2.] )
            elif expname == 's4wide_scaled_aso' or expname == 's4wide_scaled_aso_plus_fulls4scaledsobaseline':
                scaling_factors = np.sqrt(2.) * np.sqrt( [3.1, 3.1, 2.6, 2.5, 2.0, 1.8] )
            elif expname == 's4wide_single_chlat':
                scaling_factors = np.sqrt(2.) * np.sqrt( [1., 1., 1., 1., 1., 1.] )
            elif expname == 's4wide_single_chlat_plus_aso' or expname == 's4wide_single_chlat_plus_aso_plus_fulls4scaledsobaseline':
                scaling_factors_1 = np.sqrt(2.) * np.sqrt( [1., 1., 1., 1., 1., 1.] ) #single ch-lat
                scaling_factors_2 = np.sqrt(2.) * np.sqrt( [3.1, 3.1, 2.6, 2.5, 2.0, 1.8] ) #advanced-SO
            #NOTE: 's4wide_plus_fulls4scaledsobaseline' is named below but is unreachable. The guard at the top of this branch admits only 's4wide' exactly, 's4wide_scaled*' and 's4wide_single*'.
            elif expname == 's4wide' or expname == 's4wide_plus_fulls4scaledsobaseline':
                scaling_factors = np.asarray( [1., 1., 1., 1., 1., 1.] )
            else:
                raise ValueError('No detector scaling factors are defined for expname %s' % (expname))

            #specs_dic_ori = copy.deepcopy(specs_dic)  #unused
            if expname == 's4wide_single_chlat_plus_aso' or expname == 's4wide_single_chlat_plus_aso_plus_fulls4scaledsobaseline':
                for nucntr, nu in enumerate( sorted(specs_dic) ):  #pair by ascending frequency, not dict insertion order
                    delta_t_1 = specs_dic[nu][1] * scaling_factors_1[nucntr]
                    delta_t_2 = specs_dic[nu][1] * scaling_factors_2[nucntr]

                    delta_p_1 = specs_dic[nu][4] * scaling_factors_1[nucntr]
                    delta_p_2 = specs_dic[nu][4] * scaling_factors_2[nucntr]

                    delta_t = (1./delta_t_1**2. + 1./delta_t_2**2.)**-0.5
                    delta_p = (1./delta_p_1**2. + 1./delta_p_2**2.)**-0.5

                    specs_dic[nu][1] = delta_t
                    specs_dic[nu][4] = delta_p
            else:
                for nucntr, nu in enumerate( sorted(specs_dic) ):  #pair by ascending frequency, not dict insertion order
                    specs_dic[nu][1] *= scaling_factors[nucntr]
                    specs_dic[nu][4] *= scaling_factors[nucntr]

            #20220222 - modify S4 noise levels based on S4/SO detector scalings

        elif expname.find('s4_all_chile_config_lat') > -1: #20250504
            #https://github.com/sriniraghunathan/cmb_s4_all_chile_optimisation_2025_05/blob/main/get_patches_for_cmbs4.ipynb
            specs_dic = {
            #freq: [beam_arcmins, white_noise_T, elknee_T, alphaknee_T, whitenoise_P, elknee_P, alphaknee_P]
            27: [7.4, 15.8, 415., 3.5, 22.3, 700, 1.4],
            39: [5.1, 8.5, 391., 3.5, 12., 700, 1.4],
            93: [2.2, 1.4, 1932., 3.5, 2., 700, 1.4],
            145: [1.4, 1.3, 3917., 3.5, 1.9, 700, 1.4],
            225: [1.0, 4.7, 6740., 3.5, 6.6, 700, 1.4],
            278: [0.9, 13.7, 6792., 3.5, 19.4, 700, 1.4],
            }

            total_obs_time_default_for_s4_all_chile_config = 10.
            total_obs_time_default_for_advanced_so = 9.

            #s4_all_chile_config_noise_val_dic_fname = '../data/cmbs4_chile_opt_survey_patch_noise_levels.npy'
            s4_all_chile_config_noise_val_dic_fname = '%s/cmbs4_chile_opt_survey_patch_noise_levels.npy' % (data_folder)
            s4_all_chile_config_noise_val_dic = np.load(s4_all_chile_config_noise_val_dic_fname, allow_pickle=True).item()
            ##print(s4_all_chile_config_noise_val_dic.keys())

            tmp_expname = expname.split('+')[0].replace('s4_all_chile_config_', '')
            tmpsplit = tmp_expname.split('---')
            if len(tmpsplit) == 2:
                s4_all_chile_config_survey_keyname, s4_all_chile_config_survey_patchno = tmpsplit
                s4_all_chile_config_survey_yearval = 'year%s' % (total_obs_time_default_for_s4_all_chile_config)
            elif len(tmpsplit) == 3:
                s4_all_chile_config_survey_keyname, s4_all_chile_config_survey_patchno, s4_all_chile_config_survey_yearval = tmpsplit
            else:
                raise ValueError("Expected 'survey---patchN' or 'survey---patchN---yearY' in expname %s" % (expname))
            try:
                s4_all_chile_config_survey_patchno = int( s4_all_chile_config_survey_patchno.replace('patch', '') )
            except ValueError:
                raise ValueError("Expected 'patchN' with N an integer in expname %s" % (expname))
            #if s4_all_chile_config_survey_yearval is not None:  #never None; both branches above assign a string
            s4_all_chile_config_survey_yearval = float( s4_all_chile_config_survey_yearval.replace('year', '') )
            if s4_all_chile_config_survey_yearval <= 0.:
                raise ValueError('Observing time must be positive; got %s in expname %s' % (s4_all_chile_config_survey_yearval, expname))
            s4_all_chile_config_noise_val_dic_curr_survey_noise_dic = s4_all_chile_config_noise_val_dic[s4_all_chile_config_survey_keyname][s4_all_chile_config_survey_patchno]
            mod_nu_dic = {27: 30, 39: 40, 93: 90, 145: 150, 225: 220, 278: 280}

            year_scaling = np.sqrt( total_obs_time_default_for_s4_all_chile_config / s4_all_chile_config_survey_yearval )

            for nu in specs_dic:
                mod_nu = mod_nu_dic[nu]
                noiseval_t = s4_all_chile_config_noise_val_dic_curr_survey_noise_dic[mod_nu]
                noiseval_p = noiseval_t * np.sqrt(2.)

                specs_dic[nu][1] = noiseval_t * year_scaling
                specs_dic[nu][4] = noiseval_p * year_scaling

            ##print(specs_dic)
            #combine with ASO if desired
            if expname.find('+advanced_so') > -1: #20250504
                #https://arxiv.org/pdf/2503.00636
                #Page 44 for 1/f definitions.
                aso_specs_dic = {
                #freq: [beam_arcmins, white_noise_T, elknee_T, alphaknee_T, whitenoise_P, elknee_P, alphaknee_P]
                27: [7.4, 61., 500., 3.5, None, 700, 1.4],
                39: [5.1, 30., 500., 3.5, None, 700, 1.4],
                93: [2.2, 5.3, 2100., 3.5, None, 700, 1.4],
                145: [1.4, 6.6, 3000., 3.5, None, 700, 1.4],
                225: [1.0, 15., 3800., 3.5, None, 700, 1.4],
                278: [0.9, 35., 3800., 3.5, None, 700, 1.4],
                }
                #NOTE: The dict above carries the Advanced SO knees from arXiv:2503.00636, but it is immediately overwritten here, so the SO-like LAT atmosphere is modelled with S4 knees.

                aso_specs_dic = { #same as S4 for 1/f
                27: [7.4, None, 415., 3.5, None, 700, 1.4],
                39: [5.1, None, 391., 3.5, None, 700, 1.4],
                93: [2.2, None, 1932., 3.5, None, 700, 1.4],
                145: [1.4, None, 3917., 3.5, None, 700, 1.4],
                225: [1.0, None, 6740., 3.5, None, 700, 1.4],
                278: [0.9, None, 6792., 3.5, None, 700, 1.4],
                }

                if expname.find('+advanced_so_baseline') > -1:
                    noise_arr_t = np.asarray( [61., 30., 5.3, 6.6, 15., 35.])
                elif expname.find('+advanced_so_goal') > -1:
                    noise_arr_t = np.asarray( [44., 23., 3.8, 4.1, 10., 25.])
                else:
                    raise ValueError("Expected '+advanced_so_baseline' or '+advanced_so_goal' in expname %s" % (expname))
                noise_arr_p = noise_arr_t * np.sqrt(2.)

                tmp_expname = expname.split('+')[1] #aso stuff
                tmpsplit = tmp_expname.split('---')
                if len(tmpsplit) == 1:
                    aso_yearval = total_obs_time_default_for_advanced_so
                elif len(tmpsplit) == 2:
                    aso_yearval = float( tmpsplit[1].replace('year', '') )
                else:
                    raise ValueError('Cannot parse the Advanced SO observing time from expname %s' % (expname))
                if aso_yearval <= 0.:
                    raise ValueError('Observing time must be positive; got %s in expname %s' % (aso_yearval, expname))

                aso_year_scaling = np.sqrt( total_obs_time_default_for_advanced_so / aso_yearval )
                ##print( s4_all_chile_config_survey_yearval, year_scaling, aso_yearval, aso_year_scaling)
                ##sys.exit()

                for nucntr, nu in enumerate( sorted(specs_dic) ):  #pair by ascending frequency, not dict insertion order
                    noiseval_t, noiseval_p = noise_arr_t[nucntr], noise_arr_p[nucntr]

                    aso_specs_dic[nu][1] = noiseval_t * aso_year_scaling
                    aso_specs_dic[nu][4] = noiseval_p * aso_year_scaling

                ##print(aso_specs_dic); ##sys.exit()

                #combine S4 and ASO now
                for nucntr, nu in enumerate( sorted(specs_dic) ):
                    noise_t_1, noise_t_2 = specs_dic[nu][1], aso_specs_dic[nu][1]
                    final_t_noise_level = ( 1./noise_t_1**2. + 1./noise_t_2**2. )**-0.5

                    noise_p_1, noise_p_2 = specs_dic[nu][4], aso_specs_dic[nu][4]
                    final_p_noise_level = ( 1./noise_p_1**2. + 1./noise_p_2**2. )**-0.5

                    specs_dic[nu][1] = final_t_noise_level
                    specs_dic[nu][4] = final_p_noise_level

                    ##print( noise_t_1, noise_t_2, final_t_noise_level )
                    ##print( noise_p_1, noise_p_2, final_p_noise_level )
                    ##sys.exit()

                #NOTE: This corr_noise/corr_noise_bands/rho block is not used. Control falls through to the 'corr_noise = 1' assignment at the end of the s4 branch, which overrides all three for every S4-family expname.
                corr_noise = 0
                if corr_noise:
                    corr_noise_bands = {27: [39], 39: [27], 93: [145], 145: [93], 225: [278], 278: [225], 350: [350]}
                else:
                    corr_noise_bands = {27: [27], 39: [39], 93: [93], 145: [145], 225: [225], 278: [278], 350: [350]}
                rho = 0.9

                ##print(specs_dic); sys.exit()

        elif expname == 's4wide_acheived_performance_pbdr_202312xx': #20231213
            ##PBDR achieved performance. PBDR A.1.3.3
            specs_dic = {
            #freq: [beam_arcmins, white_noise_T, elknee_T, alphaknee_T, whitenoise_P, elknee_P, alphaknee_P]
            27: [7.4, 15.8, 415., 3.5, 22.3, 700, 1.4],
            39: [5.1, 8.5, 391., 3.5, 12., 700, 1.4],
            93: [2.2, 1.4, 1932., 3.5, 2., 700, 1.4],
            145: [1.4, 1.3, 3917., 3.5, 1.9, 700, 1.4],
            225: [1.0, 4.7, 6740., 3.5, 6.6, 700, 1.4],
            278: [0.9, 13.7, 6792., 3.5, 19.4, 700, 1.4],
            }

        elif expname == 's4wide_202310xx_pbdr_config': #20231025
            ##New PBDR (Oct 2023) from https://docs.google.com/spreadsheets/d/10fL76XTzhgP_B_GKsEW4nqNTkRgvp2dh4zYh6Y-G2AE/edit#gid=0
            specs_dic = {
            #freq: [beam_arcmins, white_noise_T, elknee_T, alphaknee_T, whitenoise_P, elknee_P, alphaknee_P]
            27: [7.8, 27.1, 415., 3.5, 30.23 * (27.1/21.8), 700, 1.4],
            39: [5.3, 11.6, 391., 3.5, 16.53 * (11.6/12.4), 700, 1.4],
            93: [2.2, 2.0, 1932., 3.5, 2.68, 700, 1.4],
            145: [1.4, 2.0, 3917., 3.5, 2.96, 700, 1.4],
            225: [1.0, 6.9, 6740., 3.5, 9.78, 700, 1.4],
            278: [0.9, 16.7, 6792., 3.5, 23.93, 700, 1.4],
            }

        elif expname == 's4deepv3r025_202310xx_pbdr_config': #20231025
            ##New PBDR (Oct 2023) from https://docs.google.com/spreadsheets/d/10fL76XTzhgP_B_GKsEW4nqNTkRgvp2dh4zYh6Y-G2AE/edit#gid=0
            specs_dic = {
            #freq: [beam_arcmins, white_noise_T, elknee_T, alphaknee_T, whitenoise_P, elknee_P, alphaknee_P]
            20: [11.4, 11.9, 1200., 4.2, 13.16 * (11.9/9.4), 150., 2.7], #20230517
            27: [9.1, 6.5, 1200., 4.2, 6.5 * (6.5/4.6), 150, 2.7],
            39: [6.2, 3.0, 1200., 4.2, 4.15, 150, 2.7],
            93: [2.5, 0.45, 1200., 4.2, 0.63, 150, 2.6],
            145: [1.6, 0.41, 1900., 4.1, 0.59, 200, 2.2],
            225: [1.1, 1.3, 2100., 4.1, 1.83, 200, 2.2],
            278: [1.0, 3.1, 2100., 3.9, 4.34, 200, 2.2],
            }

        elif expname == 's4wide_chlat_el40':
            #https://cmb-s4.atlassian.net/wiki/spaces/XC/pages/680853505/Neff+forecasts+for+CHLAT+for+different+observing+elevations
            specs_dic = {
            #freq: [beam_arcmins, white_noise_T, elknee_T, alphaknee_T, whitenoise_P, elknee_P, alphaknee_P]
            27: [7.3, 12.96, 415., 3.5, 18.33, 700, 1.4],
            39: [5.5, 9.83, 391., 3.5, 13.90, 700, 1.4],
            93: [2.3, 1.49, 1932., 3.5, 2.11, 700, 1.4],
            145: [1.5, 1.39, 3917., 3.5, 1.97, 700, 1.4],
            225: [1.0, 3.44, 6740., 3.5, 4.87, 700, 1.4],
            278: [0.8, 5.44, 6792., 3.5, 7.69, 700, 1.4],
            }

        elif expname.find('cmbhd') > -1 or expname.find('cmb-hd') > -1:
            specs_dic = {
            #freq: [beam_arcmins, white_noise_T, elknee_T, alphaknee_T, whitenoise_P, elknee_P, alphaknee_P]
            27: [1.4, 6.5, 415., 3.5, 9.2, 700, 1.4],
            39: [1.05, 3.4, 391., 3.5, 4.8, 700, 1.4],
            93: [0.45, 0.73, 1932., 3.5, 1.03, 700, 1.4],
            145: [0.25, 0.79, 3917., 3.5, 1.12, 700, 1.4],
            225: [0.2, 2.0, 6740., 3.5, 2.828, 700, 1.4],
            278: [0.15, 4.6, 6792., 3.5, 6.5, 700, 1.4],
            350: [0.12, 4.6, 6792., 3.5, 6.5, 700, 1.4],
            }

        elif expname == 's4deep':
            '''
            specs_dic = {
            #freq: [beam_arcmins, white_noise_T, elknee_T, alphaknee_T, whitenoise_P, elknee_P, alphaknee_P]
            #20: [10.0, None, None, None, None, None, None],
            27: [7.4, 21.8, 471., 3.5, 30.8, 700, 1.4],
            39: [5.1, 12.4, 428., 3.5, 17.6, 700, 1.4],
            93: [2.2, 0.48, 2154., 3.5, 0.68, 700, 1.4],
            145: [1.4, 0.67, 4364., 3.5, 0.96, 700, 1.4],
            225: [1.0, 4.04, 7334., 3.5, 5.72, 700, 1.4],
            278: [0.9, 6.92, 7308., 3.5, 9.8, 700, 1.4],
            }
            '''
            #different 1/f noise definitions from S4V3R0, I think.
            specs_dic = {
            #freq: [beam_arcmins, white_noise_T, elknee_T, alphaknee_T, whitenoise_P, elknee_P, alphaknee_P]
            #20: [10.0, None, None, None, None, None, None],
            27: [7.4, 21.8, 1200., 4.2, 30.8, 700, 1.4],
            39: [5.1, 12.4, 1200., 4.2, 17.6, 700, 1.4],
            93: [2.2, 0.48, 1200., 4.2, 0.68, 700, 1.4],
            145: [1.4, 0.67, 1900., 4.1, 0.96, 700, 1.4],
            225: [1.0, 4.04, 2100., 4.1, 5.72, 700, 1.4],
            278: [0.9, 6.92, 2100., 3.9, 9.8, 700, 1.4],
            }

            '''
            specs_dic = {
            #freq: [beam_arcmins, white_noise_T, elknee_T, alphaknee_T, whitenoise_P, elknee_P, alphaknee_P]
            #20: [10.0, None, None, None, None, None, None],
            #27: [7.4, 21.8, 471., 3.5, 30.8, 700, 1.4],
            #39: [5.1, 12.4, 428., 3.5, 17.6, 700, 1.4],
            93: [2.2, 0.48, 2154., 0., 0.68, 700, 1.4],
            145: [1.4, 0.67, 4364., 0., 0.96, 700, 1.4],
            225: [1.0, 4.04, 7334., 0., 5.72, 700, 1.4],
            278: [0.9, 6.92, 7308., 0., 9.8, 700, 1.4],
            }
            '''

        elif expname == 's4deepv3r025' or expname == 's4deepv3r025_tma':
            #https://cmb-s4.org/wiki/index.php/Delensing_sensitivity_-_updated_sensitivities,_beams,_TT_noise
            #20230517 - changing this to match PBDR
            '''
            specs_dic = {
            #freq: [beam_arcmins, white_noise_T, elknee_T, alphaknee_T, whitenoise_P, elknee_P, alphaknee_P]
            #20: [11.0, None, None, None, None, None, None],
            27: [8.4, 4.6, 1200., 4.2, 6.5, 150, 2.7],
            39: [5.8, 2.94, 1200., 4.2, 4.15, 150, 2.7],
            93: [2.5, 0.45, 1200., 4.2, 0.63, 150, 2.6],
            145: [1.6, 0.41, 1900., 4.1, 0.59, 200, 2.2],
            225: [1.1, 1.29, 2100., 4.1, 1.83, 200, 2.2],
            278: [1.0, 3.07, 2100., 3.9, 4.34, 200, 2.2],
            }
            '''

            specs_dic = {
            #freq: [beam_arcmins, white_noise_T, elknee_T, alphaknee_T, whitenoise_P, elknee_P, alphaknee_P]
            20: [11.4, 9.4, 1200., 4.2, 13.16, 150., 2.7], #20230517
            27: [8.4, 4.6, 1200., 4.2, 6.5, 150, 2.7],
            39: [5.8, 3.0, 1200., 4.2, 4.15, 150, 2.7],
            93: [2.5, 0.45, 1200., 4.2, 0.63, 150, 2.6],
            145: [1.6, 0.41, 1900., 4.1, 0.59, 200, 2.2],
            225: [1.1, 1.3, 2100., 4.1, 1.83, 200, 2.2],
            278: [1.0, 3.1, 2100., 3.9, 4.34, 200, 2.2],
            }

            if expname == 's4deepv3r025_tma':
                s4deep_CD_dia = 6. #metres
                s4deep_TMA_dia = 5. #metres
                for nu in specs_dic:
                    specs_dic[nu][0] = specs_dic[nu][0] * s4deep_CD_dia/s4deep_TMA_dia
        else:
            raise ValueError('Unrecognized expname: %s' % (expname))

        #freqarr = sorted( specs_dic.keys() )  #unused
        #nc = len( freqarr )  #unused

        corr_noise = 1
        if corr_noise:
            corr_noise_bands = {20: [20], 27: [39], 39: [27], 93: [145], 145: [93], 225: [278], 278: [225], 350: [350]}
        else:
            corr_noise_bands = {20: [20], 27: [27], 39: [39], 93: [93], 145: [145], 225: [225], 278: [278], 350: [350]}
        rho = 0.9

    elif expname.find('atlast') > -1:
        if remove_atm:
            warnings.warn("remove_atm is not implemented for 'atlast' configurations; the atmospheric term is retained", stacklevel=2)
        if expname == 'atlast_dummy': #cmb-hd
            specs_dic = {
            #freq: [beam_arcmins, white_noise_T, elknee_T, alphaknee_T, whitenoise_P, elknee_P, alphaknee_P]
            42: [1.05, 3.4, 391., 3.5, 4.8, 700, 1.4],
            92: [0.45, 0.73, 1932., 3.5, 1.03, 700, 1.4],
            151: [0.25, 0.79, 3917., 3.5, 1.12, 700, 1.4],
            217: [0.2, 2.0, 6740., 3.5, 2.828, 700, 1.4],
            288: [0.15, 4.6, 6792., 3.5, 6.5, 700, 1.4],
            350: [0.12, 4.6, 6792., 3.5, 6.5, 700, 1.4],
            }
        elif expname == 'atlast':
            specs_dic = {
            #freq: [beam_arcsecs, white_noise_T, elknee_T, alphaknee_T, whitenoise_P, elknee_P, alphaknee_P]
            42: [35.34, 2.32, 391., 3.5, 4.8, 700, 1.4],
            92: [16.22, 1.09, 1932., 3.5, 1.03, 700, 1.4],
            151: [9.83, 1.1, 3917., 3.5, 1.12, 700, 1.4],
            217: [6.82, 1.72, 6740., 3.5, 2.828, 700, 1.4],
            288: [5.14, 3.61, 6792., 3.5, 6.5, 700, 1.4],
            350: [4.24, 10.3, 6792., 3.5, 6.5, 700, 1.4],
            403: [3.68, 29.37, 6792., 3.5, 6.5, 700, 1.4],
            654: [2.27, 1397.27, 6792., 3.5, 6.5, 700, 1.4],
            845: [1.76, 29852.36, 6792., 3.5, 6.5, 700, 1.4],
            }

            for nu in specs_dic:
                specs_dic[nu][0] = specs_dic[nu][0] / 60. #arcsecs to arcmins

        else:
            raise ValueError('Unrecognized expname: %s' % (expname))

        rho = 0.9
        corr_noise = 0
        corr_noise_bands = {}
        for nu in specs_dic:
            corr_noise_bands[nu] = [nu]

    elif expname.find('advanced_so') > -1: #20250504
        if remove_atm:
            warnings.warn('remove_atm is not implemented for Advanced SO configurations; the atmospheric term is retained', stacklevel=2)
        #https://arxiv.org/pdf/2503.00636
        #Page 44 for 1/f definitions.
        specs_dic = {
        #freq: [beam_arcmins, white_noise_T, elknee_T, alphaknee_T, whitenoise_P, elknee_P, alphaknee_P]
        27: [7.4, 61., 500., 3.5, None, 700, 1.4],
        39: [5.1, 30., 500., 3.5, None, 700, 1.4],
        93: [2.2, 5.3, 2100., 3.5, None, 700, 1.4],
        145: [1.4, 6.6, 3000., 3.5, None, 700, 1.4],
        225: [1.0, 15., 3800., 3.5, None, 700, 1.4],
        278: [0.9, 35., 3800., 3.5, None, 700, 1.4],
        }
        #NOTE: The dict above carries the Advanced SO knees from arXiv:2503.00636, but it is immediately overwritten here, so the SO-like LAT atmosphere is modelled with S4 knees.

        specs_dic = { #same as S4 for 1/f
        27: [7.4, None, 415., 3.5, None, 700, 1.4],
        39: [5.1, None, 391., 3.5, None, 700, 1.4],
        93: [2.2, None, 1932., 3.5, None, 700, 1.4],
        145: [1.4, None, 3917., 3.5, None, 700, 1.4],
        225: [1.0, None, 6740., 3.5, None, 700, 1.4],
        278: [0.9, None, 6792., 3.5, None, 700, 1.4],
        }

        total_obs_time_default_for_advanced_so = 9.
        tmpsplit = expname.split('---')
        if len(tmpsplit) == 1:
            aso_yearval = total_obs_time_default_for_advanced_so
        elif len(tmpsplit) == 2:
            aso_yearval = float( tmpsplit[1].replace('year', '') )
        else:
            raise ValueError('Cannot parse the Advanced SO observing time from expname %s' % (expname))
        if aso_yearval <= 0.:
            raise ValueError('Observing time must be positive; got %s in expname %s' % (aso_yearval, expname))

        if expname.find('advanced_so_baseline') > -1:
            noise_arr_t = np.asarray( [61., 30., 5.3, 6.6, 15., 35.])
        elif expname.find('advanced_so_goal') > -1:
            noise_arr_t = np.asarray( [44., 23., 3.8, 4.1, 10., 25.])
        else:
            raise ValueError("Expected 'advanced_so_baseline' or 'advanced_so_goal' in expname %s" % (expname))
        noise_arr_p = noise_arr_t * np.sqrt(2.)

        aso_year_scaling = np.sqrt( total_obs_time_default_for_advanced_so / aso_yearval )
        ##print( nu, aso_yearval, aso_year_scaling); ##sys.exit()

        for nucntr, nu in enumerate( sorted(specs_dic) ):  #pair by ascending frequency, not dict insertion order
            noiseval_t, noiseval_p = noise_arr_t[nucntr], noise_arr_p[nucntr]

            specs_dic[nu][1] = noiseval_t * aso_year_scaling
            specs_dic[nu][4] = noiseval_p * aso_year_scaling

        ##print(specs_dic); sys.exit()

        corr_noise = 0
        if corr_noise:
            corr_noise_bands = {27: [39], 39: [27], 93: [145], 145: [93], 225: [278], 278: [225], 350: [350]}
        else:
            corr_noise_bands = {27: [27], 39: [39], 93: [93], 145: [145], 225: [225], 278: [278], 350: [350]}
        rho = 0.9

    #NOTE: 'ccat_prime_so' enters this branch but no white_noise_T_* values are defined for it, because its clause on the 'sobaseline' test below is commented out. It therefore raises ValueError. Update or reimplement.
    elif expname.lower() == 'sobaseline' or expname.lower() == 'sogoal' or expname.lower() == 'ccat_prime_so':
        if not remove_atm:
            '''
            elknee_T_27, elknee_T_39, elknee_T_90, elknee_T_150, elknee_T_225, elknee_T_278 = 415., 391., 1932., 3917., 6740., 6792.
            alphaknee_T_27, alphaknee_T_39, alphaknee_T_90, alphaknee_T_150, alphaknee_T_225, alphaknee_T_278 = 3.5, 3.5, 3.5, 3.5, 3.5, 3.5
            elknee_P_27, elknee_P_39, elknee_P_90, elknee_P_150, elknee_P_225, elknee_P_278 = 700., 700., 700., 700., 700., 700.
            alphaknee_P_27, alphaknee_P_39, alphaknee_P_90, alphaknee_P_150, alphaknee_P_225, alphaknee_P_278 = 1.4, 1.4, 1.4, 1.4, 1.4, 1.4
            '''
            elknee_T_27, elknee_T_39, elknee_T_90, elknee_T_150, elknee_T_225, elknee_T_278 = 1000., 1000., 1000., 1000., 1000., 1000.
            alphaknee_T_27, alphaknee_T_39, alphaknee_T_90, alphaknee_T_150, alphaknee_T_225, alphaknee_T_278 = 3.5, 3.5, 3.5, 3.5, 3.5, 3.5
            elknee_P_27, elknee_P_39, elknee_P_90, elknee_P_150, elknee_P_225, elknee_P_278 = 700., 700., 700., 700., 700., 700.
            alphaknee_P_27, alphaknee_P_39, alphaknee_P_90, alphaknee_P_150, alphaknee_P_225, alphaknee_P_278 = 1.4, 1.4, 1.4, 1.4, 1.4, 1.4

        else:
            elknee_T_27, elknee_T_39, elknee_T_90, elknee_T_150, elknee_T_225, elknee_T_278 = -1., -1., -1., -1., -1., -1.
            alphaknee_T_27, alphaknee_T_39, alphaknee_T_90, alphaknee_T_150, alphaknee_T_225, alphaknee_T_278 = 0.0, 0., 0., 0., 0., 0.
            elknee_P_27, elknee_P_39, elknee_P_90, elknee_P_150, elknee_P_225, elknee_P_278 = -1., -1., -1., -1., -1., -1.
            alphaknee_P_27, alphaknee_P_39, alphaknee_P_90, alphaknee_P_150, alphaknee_P_225, alphaknee_P_278 = 0.0, 0., 0., 0., 0., 0.

        if expname.lower() == 'sobaseline':# or expname.lower() == 'ccat_prime_so':
            white_noise_T_27, white_noise_T_39 = 71., 36.
            white_noise_T_90, white_noise_T_150 = 8., 10.
            white_noise_T_225, white_noise_T_278 = 22., 54.
        elif expname.lower() == 'sogoal':
            white_noise_T_27, white_noise_T_39 = 52., 27.
            white_noise_T_90, white_noise_T_150 = 5.8, 6.3
            white_noise_T_225, white_noise_T_278 = 15., 37.
        else:
            raise ValueError('No white-noise levels are defined for expname %s' % (expname))

        white_noise_P_27, white_noise_P_39 = white_noise_T_27 * np.sqrt(2.), white_noise_T_39 * np.sqrt(2.)
        white_noise_P_90, white_noise_P_150 = white_noise_T_90 * np.sqrt(2.), white_noise_T_150 * np.sqrt(2.)
        white_noise_P_225, white_noise_P_278 = white_noise_T_225 * np.sqrt(2.), white_noise_T_278 * np.sqrt(2.)

        beam_27, beam_39, beam_90, beam_150, beam_225, beam_278 = 7.4, 5.1, 2.2, 1.4, 1.0, 0.9

        specs_dic = {
        27: [beam_27, white_noise_T_27, elknee_T_27, alphaknee_T_27, white_noise_P_27, elknee_P_27, alphaknee_P_27],
        39: [beam_39, white_noise_T_39, elknee_T_39, alphaknee_T_39, white_noise_P_39, elknee_P_39, alphaknee_P_39],
        93: [beam_90, white_noise_T_90, elknee_T_90, alphaknee_T_90, white_noise_P_90, elknee_P_90, alphaknee_P_90],
        145: [beam_150, white_noise_T_150, elknee_T_150, alphaknee_T_150, white_noise_P_150, elknee_P_150, alphaknee_P_150],
        225: [beam_225, white_noise_T_225, elknee_T_225, alphaknee_T_225, white_noise_P_225, elknee_P_225, alphaknee_P_225],
        278: [beam_278, white_noise_T_278, elknee_T_278, alphaknee_T_278, white_noise_P_278, elknee_P_278, alphaknee_P_278],
        }

        #uK^2 seconds.
        Nred_dic[27] = [100., -1.]
        Nred_dic[39] = [39., -1.]
        Nred_dic[93] = [230., -1.]
        Nred_dic[145] = [1500., -1.]
        Nred_dic[225] = [17000., -1.]
        Nred_dic[278] = [31000., -1.]

        #freqarr = sorted( specs_dic.keys() )  #unused
        #nc = len( freqarr )  #unused

        corr_noise = 1
        if corr_noise:
            corr_noise_bands = {20: [20], 27: [39], 39: [27], 93: [145], 145: [93], 225: [278], 278: [225], 350: [350]}
        else:
            corr_noise_bands = {20: [20], 27: [27], 39: [39], 93: [93], 145: [145], 225: [225], 278: [278], 350: [350]}
        rho = 0.9

    elif expname.find('spt') > -1:
        elknee_T_90, elknee_T_150, elknee_T_220 = 1200., 2200., 2300.
        alphaknee_T_90, alphaknee_T_150, alphaknee_T_220 = 3., 4., 4.
        elknee_P_90, elknee_P_150, elknee_P_220 = 300., 300., 300.
        alphaknee_P_90, alphaknee_P_150, alphaknee_P_220 = 1., 1., 1.

        if remove_atm:
            alphaknee_T_90, alphaknee_T_150, alphaknee_T_220 = 0., 0., 0.
            alphaknee_P_90, alphaknee_P_150, alphaknee_P_220 = 0., 0., 0.

        beam_90, beam_150, beam_220 = 1.7, 1.2, 1.0

        if expname == 'sptsz':
            #slide 3 of Lindsey https://pole.uchicago.edu/sptpol/images/Spt_des_clusters_july17_f2f.pdf
            white_noise_T_90, white_noise_T_150, white_noise_T_220 = 40., 17., 80.
            white_noise_P_90, white_noise_P_150, white_noise_P_220 = 1e6, 1e6, 1e6

        elif expname == 'sptpolultradeep':
            #slide 3 of Lindsey https://pole.uchicago.edu/sptpol/images/Spt_des_clusters_july17_f2f.pdf
            white_noise_T_90, white_noise_T_150, white_noise_T_220 = 10., 3.8, 1e6

            white_noise_P_90 = white_noise_T_90 * np.sqrt(2.)
            white_noise_P_150 = white_noise_T_150 * np.sqrt(2.)
            white_noise_P_220 = white_noise_T_220 * np.sqrt(2.)

        elif expname == 'sptpolultradeepplus3g':
            #slide 3 of Lindsey https://pole.uchicago.edu/sptpol/images/Spt_des_clusters_july17_f2f.pdf
            white_noise_T_90, white_noise_T_150, white_noise_T_220 = 10., 3.8, 30.

            white_noise_P_90 = white_noise_T_90 * np.sqrt(2.)
            white_noise_P_150 = white_noise_T_150 * np.sqrt(2.)
            white_noise_P_220 = white_noise_T_220 * np.sqrt(2.)

        elif expname == 'sptpolplusultradeep':
            #slide 3 of Lindsey https://pole.uchicago.edu/sptpol/images/Spt_des_clusters_july17_f2f.pdf
            white_noise_T_90, white_noise_T_220 = 10., 1e6
            white_noise_T_150 = ( (1./3.8**2.) + (1./5**2.) )**-0.5

            white_noise_P_90 = white_noise_T_90 * np.sqrt(2.)
            white_noise_P_150 = white_noise_T_150 * np.sqrt(2.)
            white_noise_P_220 = white_noise_T_220 * np.sqrt(2.)

        elif expname == 'sptpolplusultradeepplus3g':
            #slide 3 of Lindsey https://pole.uchicago.edu/sptpol/images/Spt_des_clusters_july17_f2f.pdf
            white_noise_T_90, white_noise_T_220 = 10., 30.
            white_noise_T_150 = ( (1./3.8**2.) + (1./5**2.) )**-0.5

            white_noise_P_90 = white_noise_T_90 * np.sqrt(2.)
            white_noise_P_150 = white_noise_T_150 * np.sqrt(2.)
            white_noise_P_220 = white_noise_T_220 * np.sqrt(2.)

        elif expname == 'sptpolplusultradeepplus3gfull':
            #slide 3 of Lindsey https://pole.uchicago.edu/sptpol/images/Spt_des_clusters_july17_f2f.pdf
            white_noise_T_90, white_noise_T_220 = 10., 10.
            white_noise_T_150 = ( (1./3.8**2.) + (1./5**2.) )**-0.5

            white_noise_P_90 = white_noise_T_90 * np.sqrt(2.)
            white_noise_P_150 = white_noise_T_150 * np.sqrt(2.)
            white_noise_P_220 = white_noise_T_220 * np.sqrt(2.)

        elif expname == 'sptpolsummer':
            white_noise_T_150 = 28.
            white_noise_T_90 = white_noise_T_150 * np.sqrt(2.)
            white_noise_T_220 = 1e6

            white_noise_P_90 = white_noise_T_90 * np.sqrt(2.)
            white_noise_P_150 = white_noise_T_150 * np.sqrt(2.)
            white_noise_P_220 = white_noise_T_220 * np.sqrt(2.)

        elif expname == 'spt3gsummer':
            white_noise_T_90 = 15.
            white_noise_T_150 = 15.
            white_noise_T_220 = 40.

            white_noise_P_90 = white_noise_T_90 * np.sqrt(2.)
            white_noise_P_150 = white_noise_T_150 * np.sqrt(2.)
            white_noise_P_220 = white_noise_T_220 * np.sqrt(2.)

        elif expname == 'spt3g':
            white_noise_T_90 = 3.0
            white_noise_T_150 = 2.2
            white_noise_T_220 = 8.8

            white_noise_P_90 = np.sqrt(2.) * white_noise_T_90
            white_noise_P_150 = np.sqrt(2.) * white_noise_T_150
            white_noise_P_220 = np.sqrt(2.) * white_noise_T_220

        elif expname == 'spt3g_ZP':
            white_noise_T_90 = 14.3
            white_noise_T_150 = 14.0
            white_noise_T_220 = 48.

            white_noise_P_90 = 25.
            white_noise_P_150 = 17.
            white_noise_P_220 = 56.

        elif expname == 'spt3g_y12':
            white_noise_T_90 = 9.
            white_noise_T_150 = 10.
            white_noise_T_220 = 40.

            white_noise_P_90 = np.sqrt(2.) * white_noise_T_90
            white_noise_P_150 = np.sqrt(2.) * white_noise_T_150
            white_noise_P_220 = np.sqrt(2.) * white_noise_T_220

        elif expname == 'spt3g_summer':
            white_noise_T_90 = 13.
            white_noise_T_150 = 13.0
            white_noise_T_220 = 40.

            white_noise_P_90 = np.sqrt(2.) * white_noise_T_90
            white_noise_P_150 = np.sqrt(2.) * white_noise_T_150
            white_noise_P_220 = np.sqrt(2.) * white_noise_T_220

        elif expname == 'sptpol_summer':
            white_noise_T_90 = 60.
            white_noise_T_150 = 28.0
            white_noise_T_220 = 10000.

            white_noise_P_90 = np.sqrt(2.) * white_noise_T_90
            white_noise_P_150 = np.sqrt(2.) * white_noise_T_150
            white_noise_P_220 = np.sqrt(2.) * white_noise_T_220

        elif expname == 'sptpol':
            white_noise_T_90 = 14.
            white_noise_T_150 = 6.0
            white_noise_T_220 = 10000.

            white_noise_P_90 = np.sqrt(2.) * white_noise_T_90
            white_noise_P_150 = np.sqrt(2.) * white_noise_T_150
            white_noise_P_220 = np.sqrt(2.) * white_noise_T_220

        elif expname == 'spt3g_TC':
            white_noise_T_90 = 7.960
            white_noise_T_150 = 6.330
            white_noise_T_220 = 23.39

            white_noise_P_90 = 11.26
            white_noise_P_150 = 8.950
            white_noise_P_220 = 33.08

        elif expname == 'spt3g_WG2': #202505xx WG2 of CMB-S4
            white_noise_P_90 = 3.5
            white_noise_P_150 = 3.
            white_noise_P_220 = 10.75

            white_noise_T_90 = white_noise_P_90 / np.sqrt(2.)
            white_noise_T_150 = white_noise_P_150 / np.sqrt(2.)
            white_noise_T_220 = white_noise_P_220 / np.sqrt(2.)

        elif expname == 'spt3g_plus_spt3g+_WG2': #202505xx WG2 of CMB-S4
            white_noise_P_90 = 3.5
            white_noise_P_150 = 3.
            white_noise_P_220 = 10.75

            white_noise_T_90 = white_noise_P_90 / np.sqrt(2.)
            white_noise_T_150 = white_noise_P_150 / np.sqrt(2.)
            white_noise_T_220 = white_noise_P_220 / np.sqrt(2.)

            ##print( white_noise_T_90, white_noise_T_150, white_noise_T_220)
            ##print( white_noise_P_90, white_noise_P_150, white_noise_P_220)

            #add spt3g+ to this
            white_noise_P_90_spt3gplus = 1.22
            white_noise_P_150_spt3gplus = 1.07

            white_noise_T_90_spt3gplus = white_noise_P_90_spt3gplus / np.sqrt(2.)
            white_noise_T_150_spt3gplus = white_noise_P_150_spt3gplus / np.sqrt(2.)

            white_noise_T_90 = (1./white_noise_T_90**2. + 1./white_noise_T_90_spt3gplus**2.)**-0.5
            white_noise_T_150 = (1./white_noise_T_150**2. + 1./white_noise_T_150_spt3gplus**2.)**-0.5

            white_noise_P_90 = (1./white_noise_P_90**2. + 1./white_noise_P_90_spt3gplus**2.)**-0.5
            white_noise_P_150 = (1./white_noise_P_150**2. + 1./white_noise_P_150_spt3gplus**2.)**-0.5

            ##print( white_noise_T_90, white_noise_T_150, white_noise_T_220)
            ##print( white_noise_P_90, white_noise_P_150, white_noise_P_220)
            ##sys.exit()

        else:
            raise ValueError('Unrecognized expname: %s' % (expname))

        specs_dic = {
        90: [beam_90, white_noise_T_90, elknee_T_90, alphaknee_T_90, white_noise_P_90, elknee_P_90, alphaknee_P_90],
        150: [beam_150, white_noise_T_150, elknee_T_150, alphaknee_T_150, white_noise_P_150, elknee_P_150, alphaknee_P_150],
        220: [beam_220, white_noise_T_220, elknee_T_220, alphaknee_T_220, white_noise_P_220, elknee_P_220, alphaknee_P_220],
        }

        corr_noise = corr_noise_for_spt
        if corr_noise:
            corr_noise_bands = {90: [90], 150: [220], 220: [150]}
        else:
            corr_noise_bands = {90: [90], 150: [150], 220: [220]}

        rho = 1.0

    elif expname.lower() == 'planck':
        #Planck instrumental noise, from Table 1 of Allison et al. 2015 (arXiv:1509.07471).
        #Low-frequency noise is neglected for Planck, so elknee = -1 in every band, which switches off the atmospheric term in misc.get_nl.
        #No polarization sensitivity is quoted for 30 and 44 GHz. delta_P there is set to 1e4 so that these channels are negligible in the ILC. (Not 1e6 as in the sptsz branch since that overflows to inf below the standard lmax = 6500 with the 33 arcmin beam.)
        no_pol_noise = 1e4  #uK-arcmin
        specs_dic = {
        #freq: [beam_arcmins, white_noise_T, elknee_T, alphaknee_T, whitenoise_P, elknee_P, alphaknee_P]
        30: [33., 145., -1., 0., no_pol_noise, -1., 0.],
        44: [23., 149., -1., 0., no_pol_noise, -1., 0.],
        70: [14., 137., -1., 0., 450., -1., 0.],
        100: [10., 65., -1., 0., 103., -1., 0.],
        143: [7., 43., -1., 0., 81., -1., 0.],
        217: [5., 66., -1., 0., 134., -1., 0.],
        353: [5., 200., -1., 0., 406., -1., 0.],
        }

        corr_noise = 0
        rho = 0.
        corr_noise_bands = {}
        for nu in specs_dic:
            corr_noise_bands[nu] = [nu]

    else:
        raise ValueError('Unrecognized expname: %s' % (expname))

    return specs_dic, corr_noise_bands, rho, corr_noise, Nred_dic
