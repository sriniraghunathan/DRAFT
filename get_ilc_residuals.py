r"""
Compute ILC noise and foreground residuals for a CMB experiment configuration.

This is the driver for the component separation stage of DRAFT.
It collects per-band beams and white-noise levels for a named experiment configuration, builds the analytic multi-frequency covariance of the CMB, foregrounds and instrumental noise, performs the internal linear combination, and writes the resulting residual power spectra to a ``.npy`` product file.
The Fisher forecasting and foreground bias stages of DRAFT consume those product files and are not part of this module.

Run from the repository root::

    python3 get_ilc_residuals.py -expname s4wide_202310xx_pbdr_config -include_gal 1 -which_gal_mask 2 -final_comp cmb

The same calculation is available to other scripts through :func:`run_ilc`, which accepts the command line options as keyword arguments and returns the product dictionary together with the path it would be written to::

    import get_ilc_residuals

    opdic, opfname = get_ilc_residuals.run_ilc('s4wide_202310xx_pbdr_config', include_gal=1, which_gal_mask=2, save=False)

Notes
-----
Survey specifications come from :mod:`exp_specs`, beam and noise curves from :mod:`misc`, and the covariance and ILC weights from :mod:`ilc`, which in turn draws its extragalactic and galactic foreground models from :mod:`foregrounds`.
"""

import argparse
import os
import sys
import warnings

import numpy as np

sys.path.append( os.path.join( os.path.dirname(os.path.abspath(__file__)), 'modules' ) )
import exp_specs
import misc
import ilc


# Constants

COMP_CHOICES = ['cmb', 'tsz', 'y', 'radio', 'cib']
"""
Accepted values of ``final_comp`` and ``null_comp``, matching check of :func:`ilc.residual_power`.
"""

PARAMFILE = os.path.join( os.path.dirname(os.path.abspath(__file__)), 'params.ini' )
"""
Default parameter file.
"""

GAL_SIM_BANDS = [27, 39, 93, 145, 225, 278]  #hardcoded for now; the relevant ``.npy`` file could be read instead.
"""
Bands in GHz covered by the galactic dust and synchrotron simulations under ``data/galactic/``.
"""
 
GAL_MASK_INDICES = [0, 1, 2]  #hardcoded for now, as above.
"""
Galactic mask indices carried by those simulations, indexing ``fsky_arr``.
"""
 
TP_ARR = ['T', 'P']
"""
Temperature and polarization labels, in the order ``Nred_dic`` indexes them.
"""

WHICH_SPEC_ARR = ['TT', 'EE']
"""
Spectra the ILC is solved for.
"""


# Command line

def parse_args(argv=None):
    """
    Parse the command line.

    Parameters
    ----------
    argv : list of str, optional
        Argument list to parse. Default is ``sys.argv[1:]``.

    Returns
    -------
    argparse.Namespace
        Parsed options. The attribute names match the keyword arguments of :func:`run_ilc`.
    """

    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-expname', dest='expname', action='store', help='expname', type=str, required=True)
    parser.add_argument('-total_obs_time', dest='total_obs_time', action='store', help='total_obs_time in years', type=float, default=7.0)
    parser.add_argument('-include_gal', dest='include_gal', action='store', help='include_gal', type=int, default=0)
    parser.add_argument('-which_gal_mask', dest='which_gal_mask', action='store', help='which_gal_mask', type=int, default=2)
    parser.add_argument('-interactive_mode', dest='interactive_mode', action='store', help='interactive_mode', type=int, default=1)  #diagnostic plots only; see plot_noise_spectra
    parser.add_argument('-save_fg_res_and_weights', dest='save_fg_res_and_weights', action='store', help='save_fg_res_and_weights', type=int, default=1)
    parser.add_argument('-s4_so_joint_configs', dest='s4_so_joint_configs', action='store', help='s4_so_joint_configs', type=int, default=0)
    parser.add_argument('-include_fulls4scaledsobaseline', dest='include_fulls4scaledsobaseline', action='store', help='include_fulls4scaledsobaseline', type=int, default=0)

    #20230530 - scale noise levels of bands
    parser.add_argument('-noise_scalings_for_bands', dest='noise_scalings_for_bands', action='store', help='noise_scalings_for_bands', nargs='+', type=float, default=None)

    #20230531 - option to get CMB or y
    parser.add_argument('-final_comp', dest='final_comp', action='store', help='final_comp', type=str, default='cmb', choices=COMP_CHOICES)
    parser.add_argument('-null_comp', dest='null_comp', action='store', help='null_comp', nargs='+', type=str, default=None, choices=COMP_CHOICES)
    parser.add_argument('-paramfile', dest='paramfile', action='store', help='paramfile', type=str, default=PARAMFILE)
    parser.add_argument('-debug', dest='debug', action='store', help='debug', type=int, default=0)  #print intermediate dictionary keys

    return parser.parse_args(argv)


# Experiment specifications

def get_experiment_specs(expname, remove_atm=0, include_fulls4scaledsobaseline=0):
    """
    Look up the beam and noise specifications for an experiment configuration.

    Parameters
    ----------
    expname : str
        Experiment configuration name, passed through to :func:`exp_specs.get_exp_specs`.
        Two names are special cased and combine two surveys: ``'s4deepv3r025_plus_s4wide'`` and any name containing ``'s4wide_single_chlat_plus_2028aso'``.
    remove_atm : int, optional
        When nonzero, drop the atmospheric :math:`1/f` component. Default is 0.
    include_fulls4scaledsobaseline : int, optional
        When nonzero, also look up the scaled SO-baseline specifications so that its noise can be combined in later. Default is 0.

    Returns
    -------
    specs_dic : dict
        Frequency band in GHz mapped to ``[beam_arcmins, Delta_T, l_knee_T, alpha_knee_T, Delta_P, l_knee_P, alpha_knee_P]``.
    corr_noise_bands : dict
        Band mapped to the bands its atmospheric noise is correlated with.
    rho : float
        Atmospheric noise correlation coefficient.
    corr_noise : int
        Whether correlated noise is modelled at all.
    Nred_dic : dict
        Band mapped to red-noise levels, indexed by :data:`TP_ARR`. Empty for every configuration currently defined.
    secondary_specs : dict
        Specifications for any additional survey being combined in, keyed ``'s4wide'``, ``'aso'`` or ``'fulls4scaledsobaseline'``.
        Empty when the configuration involves a single survey.

    Warns
    -----
    UserWarning
        When the CMB-S4 Ultra-deep survey has bands the CMB-S4 Wide survey does not since those bands cannot be inverse-variance combined and keep their Ultra-deep noise.
    """

    secondary_specs = {}

    #S4 specs
    if expname == 's4deepv3r025_plus_s4wide':
        specs_dic, corr_noise_bands, rho, corr_noise, Nred_dic = exp_specs.get_exp_specs('s4deepv3r025', remove_atm=remove_atm)
        specs_dic_s4wide, corr_noise_bands_s4wide, rho_s4wide, corr_noise_s4wide, Nred_dic = exp_specs.get_exp_specs('s4wide')
        secondary_specs['s4wide'] = {'specs_dic': specs_dic_s4wide, 'corr_noise_bands': corr_noise_bands_s4wide, 'rho': rho_s4wide}
        #S4 Ultra-deep has 20 GHz band which S4 Wide does not
        bands_not_in_s4wide = [f for f in sorted(specs_dic.keys()) if f not in specs_dic_s4wide]
        if len(bands_not_in_s4wide) > 0:
            warnings.warn('Bands %s are in CMB-S4 Ultra-deep but not CMB-S4 Wide; they keep their ultra-deep noise (no inverse-variance combination)' % (bands_not_in_s4wide), stacklevel=2)
    #20220225 - 2028 ASO + single-CHLAT - we will take an inverse variance combination of Nyear single CHLAT and N+1 ASO
    elif expname.find('s4wide_single_chlat_plus_2028aso') > -1:
        specs_dic_aso, corr_noise_bands_aso, rho_aso, corr_noise_aso, Nred_dic = exp_specs.get_exp_specs('s4wide_scaled_aso')
        specs_dic, corr_noise_bands, rho, corr_noise, Nred_dic = exp_specs.get_exp_specs('s4wide_single_chlat')
        secondary_specs['aso'] = {'specs_dic': specs_dic_aso}
    else:
        specs_dic, corr_noise_bands, rho, corr_noise, Nred_dic = exp_specs.get_exp_specs(expname, remove_atm=remove_atm)

    #20220223 - include full SO-baseline, if requested
    if include_fulls4scaledsobaseline:
        specs_dic_fulls4scaledsobaseline, corr_noise_bands_fulls4scaledsobaseline, rho_fulls4scaledsobaseline, corr_noise_fulls4scaledsobaseline, Nred_dic_fulls4scaledsobaseline = exp_specs.get_exp_specs('s4wide_scaled_sobaseline')
        secondary_specs['fulls4scaledsobaseline'] = {'specs_dic': specs_dic_fulls4scaledsobaseline}

    return specs_dic, corr_noise_bands, rho, corr_noise, Nred_dic, secondary_specs


def check_gal_sim_inputs(freqarr, which_gal_mask, expname):
    """
    Check that the galactic simulations cover every band and mask requested.

    Parameters
    ----------
    freqarr : list of int
        Frequency bands in GHz.
    which_gal_mask : int
        Galactic mask index.
    expname : str
        Experiment configuration name, used in the error message.

    Raises
    ------
    ValueError
        If any band is absent from :data:`GAL_SIM_BANDS`, or if ``which_gal_mask`` is not in :data:`GAL_MASK_INDICES`.

    Notes
    -----
    Without these checks a missing band or mask raises an opaque ``KeyError`` from inside :func:`ilc.get_analytic_covariance` or when the sky fraction is read.
    """

    #Note: The galactic dust/sync sims under data/galactic/ only cover the six S4-Wide bands. Any other band (20 GHz in s4deepv3r025, for example) crashes in ilc.get_analytic_covariance.
    bands_without_gal_sims = [f for f in freqarr if f not in GAL_SIM_BANDS]
    if len(bands_without_gal_sims) > 0:
        raise ValueError('-include_gal 1 is only supported for bands %s, but expname %s also has %s. Rerun with -include_gal 0.' % (GAL_SIM_BANDS, expname, bands_without_gal_sims))
    if which_gal_mask not in GAL_MASK_INDICES:
        raise ValueError('-which_gal_mask must be one of %s, got %s' % (GAL_MASK_INDICES, which_gal_mask))


def get_total_obs_time_default(expname, total_obs_time):
    r"""
    Return the reference observing time the stored noise levels correspond to.

    Parameters
    ----------
    expname : str
        Experiment configuration name.
    total_obs_time : float
        Requested observing time in years.

    Returns
    -------
    float
        Reference observing time in years, 10 for the all-Chile configurations and 7 otherwise. The white-noise levels are later scaled by

        .. math::

            \sqrt{t_\mathrm{default} / t_\mathrm{obs}}\, .

    Raises
    ------
    ValueError
        If ``total_obs_time`` is not positive, and if ``expname`` carries a ``---year<Y>`` suffix and ``total_obs_time`` disagrees with the reference time. :mod:`exp_specs` already scales the noise using that suffix, so a disagreeing ``total_obs_time`` would apply a second scaling.
    """

    if total_obs_time <= 0:
        raise ValueError('-total_obs_time must be positive, got %g' % (total_obs_time))
    total_obs_time_default = 7. ###10. #years
    if expname.find('s4_all_chile_config_lat_') > -1: #20250504
        total_obs_time_default = 10.  #NOTE: Check if true
    if expname.find('---year') > -1 and total_obs_time_default != total_obs_time:
        raise ValueError('expname %s encodes an observing time via its ---year suffix, which exp_specs already uses to scale the noise levels. -total_obs_time %g would apply a second scaling of sqrt(%g/%g). Rerun with -total_obs_time %g.' % (expname, total_obs_time, total_obs_time_default, total_obs_time, total_obs_time_default))

    return total_obs_time_default


# Output path

def build_output_path(
        expname,
        param_dict,
        freqarr,
        which_spec_arr,
        include_gal,
        which_gal_mask,
        final_comp,
        total_obs_time,
        corr_noise,
        remove_atm,
        noise_scalings_for_bands=None,
        s4_so_joint_configs=0,
        null_comp=None,
        cl_multiplier_dic=None
        ):
    """
    Build the output product filename.

    Parameters
    ----------
    expname : str
        Experiment configuration name.
    param_dict : dict
        Parameter dictionary from :func:`misc.get_param_dict`.
    freqarr : list of int
        Frequency bands in GHz.
    which_spec_arr : list of str
        Spectra being solved for, e.g. ``['TT', 'EE']``.
    include_gal : int
        Whether galactic foregrounds are included.
    which_gal_mask : int
        Galactic mask index.
    final_comp : str
        Component the ILC preserves.
    total_obs_time : float
        Observing time in years.
    corr_noise : int
        Whether correlated atmospheric noise is modelled.
    remove_atm : int
        Whether the atmospheric component was dropped.
    noise_scalings_for_bands : list of float, optional
        Per-band noise scalings. Default is ``None``.
    s4_so_joint_configs : int, optional
        Whether this is a joint S4 and SO configuration. Default is 0.
    null_comp : str or list of str, optional
        Component nulled by a constrained ILC. Default is ``None``.
    cl_multiplier_dic : dict, optional
        Per-component spectrum multipliers for a partial ILC. Default is ``None``.

    Returns
    -------
    str
        Path of the ``.npy`` product file. The containing directory is not created here; see :func:`run_ilc`.
    """

    freqarr_str = '-'.join( np.asarray( freqarr ).astype(str) )
    which_spec_arr_str = '-'.join( np.asarray( which_spec_arr ).astype(str) )
    parent_folder = 'results/20210506_with202102designtoolinputforpySM3sims_sedscalingfordust'
    #20230530
    parent_folder = '%s/202310xx_modified_PBDR_config_for_Neff_paper/' % (parent_folder) #20231025 - modified PBDR config: https://docs.google.com/spreadsheets/d/10fL76XTzhgP_B_GKsEW4nqNTkRgvp2dh4zYh6Y-G2AE/edit#gid=0

    if expname.find('s4_all_chile_config_lat_') > -1 or expname.find('advanced_so') > -1:
        #parent_folder = 'results/s4_all_chile_config'
        parent_folder = 'results/s4_all_chile_config/report/'

    if final_comp != 'cmb':
        parent_folder = '%s/%s/' % (parent_folder, final_comp)

    if noise_scalings_for_bands is not None and len(np.unique(noise_scalings_for_bands)) > 1: #20230530 - scale noise levels of bands
        parent_folder = '%s/noise_scalings/' % (parent_folder)
        noise_scalings_for_bands_str = '-'.join([str(n) for n in noise_scalings_for_bands])
        noise_scalings_for_bands_str = '_noisescalings%s' % (noise_scalings_for_bands_str)

    if s4_so_joint_configs:
        parent_folder = '%s/s4_so_joint_configs/' % (parent_folder)
    if null_comp is not None:
        null_comp_str = 'nulled_%s' % ('-'.join(null_comp))
        parent_folder = '%s/%s/' % (parent_folder, null_comp_str)

    if param_dict['lmax'] > 5002:
        parent_folder = '%s/lmax_%s/' % (parent_folder, param_dict['lmax'])
    parent_folder = '%s/%s/' % (parent_folder, expname)

    opfname = '%s/%s_ilc_galaxy%s_%s_%s.npy' % (parent_folder, expname, include_gal, freqarr_str, which_spec_arr_str)
    if null_comp is not None:
        opfname = '%s_%s.npy' % (opfname.replace('.npy', ''), null_comp_str)

    if not corr_noise:
        opfname = opfname.replace('.npy', '_nocorrnoise.npy')

    if expname.find('s4') > -1 or expname.find('cmbhd') > -1:
        opfname = opfname.replace(parent_folder, '%s/planck_mask/TT-EE/baseline/' % (parent_folder))

    if include_gal:
        opfname = opfname.replace('.npy', '_galmask%s.npy' % (which_gal_mask))

    if remove_atm:
        opfname = opfname.replace('.npy', '_noatmnoise.npy')

    if cl_multiplier_dic is not None:
        if len(cl_multiplier_dic) > 1:
            cl_multiplier_str = 'pilcmultfacs'
            for kkk in cl_multiplier_dic:
                cl_multiplier_str = '%s-%s%s' % (cl_multiplier_str, kkk, cl_multiplier_dic[kkk])
            cl_multiplier_str = cl_multiplier_str.strip('-')
            opfname = opfname.replace('.npy', '_%s.npy' % (cl_multiplier_str))

    if include_gal:
        cl_gal_folder = param_dict['cl_gal_folder']
        if cl_gal_folder.find('CUmilta') > -1:
            opfname = opfname.replace('.npy', '_CU.npy')
        else:
            opfname = opfname.replace('.npy', '_AZ.npy')

    try:
        param_dict['cl_gal_dic_sync_fname_forced']
        opfname = opfname.replace('.npy', '_forcingsynctoCU.npy')
    except KeyError:
        pass

    if param_dict['lmax'] != 5000:
        opfname = opfname.replace('.npy', '_lmax%s.npy' % (param_dict['lmax']))

    if noise_scalings_for_bands is not None and len(np.unique(noise_scalings_for_bands)) > 1: #20230530 - scale noise levels of bands
        opfname = opfname.replace('.npy', '%s.npy' % (noise_scalings_for_bands_str))

    if expname.lower().find('s4') > -1:#total_obs_time_default != total_obs_time:
        opfname = opfname.replace('.npy', '_for%gyears.npy' % (total_obs_time))

    return opfname


# Beams and noise

def get_beam_and_noise_dic(expname, specs_dic, freqarr, total_obs_time, total_obs_time_default, noise_scalings_for_bands=None, specs_aso=None, specs_fulls4scaledsobaseline=None):
    r"""
    Collect per-band beams, white-noise levels and :math:`1/f` parameters.

    Parameters
    ----------
    expname : str
        Experiment configuration name.
    specs_dic : dict
        Band specifications from :func:`get_experiment_specs`.
    freqarr : list of int
        Frequency bands in GHz.
    total_obs_time : float
        Requested observing time in years.
    total_obs_time_default : float
        Reference observing time from :func:`get_total_obs_time_default`.
    noise_scalings_for_bands : list of float, optional
        Multiplicative white-noise scaling per band, in the order of ``freqarr``. Default is ``None``.
    specs_aso : dict, optional
        Advanced SO specifications, inverse-variance combined for the ``s4wide_single_chlat_plus_2028aso`` configurations. Default is ``None``.
    specs_fulls4scaledsobaseline : dict, optional
        Scaled SO-baseline specifications to combine in. Default is ``None``.

    Returns
    -------
    beam_noise_dic : dict
        ``'T'`` and ``'P'`` mapped to band mapped to ``[beam_arcmins, Delta_X]``.
    elknee_dic : dict
        ``'T'`` and ``'P'`` mapped to band mapped to ``[l_knee, alpha_knee]``.
    noisearr_T : list of float
        White-noise levels :math:`\Delta_T` in μK arcmin.
    noisearr_P : list of float
        White-noise levels :math:`\Delta_P` in μK arcmin.

    Raises
    ------
    ValueError
        If ``noise_scalings_for_bands`` has a different length from ``freqarr``.

    Notes
    -----
    White-noise levels are scaled by :math:`\sqrt{t_\mathrm{default} / t_\mathrm{obs}}` and any additional survey is combined in as

    .. math::

        \Delta_X = \left( \Delta_{X,1}^{-2} + \Delta_{X,2}^{-2} \right)^{-1/2}\, .
    """

    #beam and noise arr
    beamarr = []
    noisearr_T, elkneearr_T, alphakneearr_T = [], [], []
    noisearr_P, elkneearr_P, alphakneearr_P = [], [], []
    for fcntr, freq in enumerate( freqarr ):
        beam_arcmins, white_noise_T, elknee_T, alphaknee_T, white_noise_P, elknee_P, alphaknee_P = specs_dic[freq]

        #20230530 - noise scaling of bands
        if noise_scalings_for_bands is not None:
            if len(noise_scalings_for_bands) != len(freqarr):
                raise ValueError('-noise_scalings_for_bands has %d entries but expname %s has %d bands %s' % (len(noise_scalings_for_bands), expname, len(freqarr), freqarr))
            white_noise_T = white_noise_T * noise_scalings_for_bands[fcntr]
            white_noise_P = white_noise_P * noise_scalings_for_bands[fcntr]

        #noise scaling based on total_obs_time
        ###print(total_obs_time, total_obs_time_default)
        noise_scaling_fac = (total_obs_time_default / total_obs_time)**0.5
        ###print(noise_scaling_fac); sys.exit()
        white_noise_T = white_noise_T * noise_scaling_fac
        white_noise_P = white_noise_P * noise_scaling_fac

        #add N+1 year aso via inverse variance combination now
        if specs_aso is not None:
            #the above noise numbers for N year single CHLAT
            #now we will add N+1 year ASO
            total_obs_time_2028_aso = total_obs_time + 1
            white_noise_T_aso, white_noise_P_aso = specs_aso['specs_dic'][freq][1], specs_aso['specs_dic'][freq][4]
            noise_scaling_fac_2028_aso = (total_obs_time_default / (total_obs_time_2028_aso))**0.5
            white_noise_T_2028_aso = white_noise_T_aso * noise_scaling_fac_2028_aso
            white_noise_P_2028_aso = white_noise_P_aso * noise_scaling_fac_2028_aso

            white_noise_T = ( (1./white_noise_T**2.) + (1./white_noise_T_2028_aso**2.) )**-0.5
            white_noise_P = ( (1./white_noise_P**2.) + (1./white_noise_P_2028_aso**2.) )**-0.5

        #20220223 - include full SO-baseline, if requested
        if specs_fulls4scaledsobaseline is not None:
            beam_arcmins_2, white_noise_T_2, elknee_T_2, alphaknee_T_2, white_noise_P_2, elknee_P_2, alphaknee_P_2 = specs_fulls4scaledsobaseline['specs_dic'][freq]
            white_noise_T_2 = white_noise_T_2 * np.sqrt(7./5.)
            white_noise_P_2 = white_noise_P_2 * np.sqrt(7./5.)
            white_noise_T = (1./white_noise_T**2. + 1./white_noise_T_2**2.)**-0.5
            white_noise_P = (1./white_noise_P**2. + 1./white_noise_P_2**2.)**-0.5

        beamarr.append(beam_arcmins)
        noisearr_T.append(white_noise_T)
        noisearr_P.append(white_noise_P)
        elkneearr_T.append(elknee_T)
        elkneearr_P.append(elknee_P)
        alphakneearr_T.append(alphaknee_T)
        alphakneearr_P.append(alphaknee_P)

    #collect beam and noise into a dic; elknee and alpha into a dic
    beam_noise_dic = {}
    elknee_dic = {}
    for TP in TP_ARR:
        beam_noise_dic[TP] = {}
        elknee_dic[TP] = {}
        if TP == 'T':
            noisearr, elkneearr, alphakneearr = noisearr_T, elkneearr_T, alphakneearr_T
        elif TP == 'P':
            noisearr, elkneearr, alphakneearr = noisearr_P, elkneearr_P, alphakneearr_P

        for (freq, beam, noise, elknee, alphaknee) in zip(freqarr, beamarr, noisearr, elkneearr, alphakneearr):
            beam_noise_dic[TP][freq] = [beam, noise]
            elknee_dic[TP][freq] = [elknee, alphaknee]

    #these used to sit commented beside the 'Delta T'/'Delta P' prints, now in run_ilc
    #print(elkneearr_T)
    #print(beamarr)

    return beam_noise_dic, elknee_dic, noisearr_T, noisearr_P


def get_noise_spectra(
        freqarr,
        el,
        param_dict,
        beam_noise_dic,
        elknee_dic,
        corr_noise_bands,
        rho,
        Nred_dic,
        specs_s4wide=None,
        total_obs_time=None,
        total_obs_time_default=None
        ):
    r"""
    Build beam deconvolved noise spectra for every band pair.

    Parameters
    ----------
    freqarr : list of int
        Frequency bands in GHz.
    el : array_like
        Multipoles.
    param_dict : dict
        Parameter dictionary. ``lmin`` is used to zero the low multipoles.
    beam_noise_dic : dict
        Beams and white-noise levels from :func:`get_beam_and_noise_dic`.
    elknee_dic : dict
        Knee multipoles and slopes from :func:`get_beam_and_noise_dic`.
    corr_noise_bands : dict
        Band mapped to the bands its atmospheric noise is correlated with.
    rho : float
        Atmospheric noise correlation coefficient.
    Nred_dic : dict
        Band mapped to red-noise levels, indexed by :data:`TP_ARR`.
    specs_s4wide : dict, optional
        S4-Wide specifications to inverse-variance combine with, holding keys ``'specs_dic'``, ``'corr_noise_bands'`` and ``'rho'``.
        If None (default), no combination is performed. Default is ``None``.
    total_obs_time : float, optional
        Observing time in years. Used only by the disabled rescaling block below. Default is ``None``.
    total_obs_time_default : float, optional
        Reference observing time in years. Used only by the disabled rescaling block below. Default is ``None``.

    Returns
    -------
    dict
        ``'T'`` and ``'P'`` mapped to ``(freq1, freq2)`` mapped to :math:`N_\ell`. Cross-spectra between uncorrelated bands are zero.

    Notes
    -----
    Where ``specs_s4wide`` is given, the two surveys are combined as

    .. math::

        N_\ell = \left( N_\ell^{-1} + N_{\ell,\mathrm{wide}}^{-1} \right)^{-1} .

    Bands present in the CMB-S4 Ultra-deep survey but absent from the CMB-S4 Wide survey have no wide-survey counterpart and keep their Ultra-deep :math:`N_\ell`.
    """

    #get beam deconvolved noise nls
    nl_dic = {}
    for TPcntr, TP in enumerate( TP_ARR ):
        nl_dic[TP] = {}
        for freq1 in freqarr:
            beamval1, noiseval1 = beam_noise_dic[TP][freq1]
            elknee1, alphaknee1 = elknee_dic[TP][freq1]
            for freq2 in freqarr:
                beamval2, noiseval2 = beam_noise_dic[TP][freq2]
                elknee2, alphaknee2 = elknee_dic[TP][freq2]

                ##elknee1, elknee2 = -1, -1.
                ##alphaknee1, alphaknee2 = -1., -1.

                Nred1, Nred2 = -1., -1.
                if freq1 in Nred_dic:
                    Nred1 = Nred_dic[freq1][TPcntr]
                if freq2 in Nred_dic:
                    Nred2 = Nred_dic[freq2][TPcntr]

                if freq1 == freq2:
                    nl = misc.get_nl(noiseval1, el, beamval1, elknee=elknee1, alphaknee=alphaknee1, Nred1=Nred1)
                else:
                    if freq2 in corr_noise_bands[freq1]:
                        nl = misc.get_nl(noiseval1, el, beamval1, elknee=elknee1, alphaknee=alphaknee1, beamval2=beamval2, noiseval2=noiseval2, elknee2=elknee2, alphaknee2=alphaknee2, rho=rho, Nred1=Nred1, Nred2=Nred2)
                    else:
                        nl = np.zeros( len(el) )
                nl[el <= param_dict['lmin']] = 0.
                ##nl[nl == 0.] = np.min(nl[nl != 0.])/1e3

                if specs_s4wide is not None and np.sum(nl) > 0. and freq1 in specs_s4wide['specs_dic'] and freq2 in specs_s4wide['specs_dic']: #20210428: inv. var. combination of S4-Wide and S4-Ultra deep nl
                    specs_dic_s4wide = specs_s4wide['specs_dic']
                    beam_arcmins_s4wide_f1, white_noise_T_s4wide_f1, elknee_T_s4wide_f1, alphaknee_T_s4wide_f1, white_noise_P_s4wide_f1, elknee_P_s4wide_f1, alphaknee_P_s4wide_f1 = np.copy(specs_dic_s4wide[freq1])
                    beam_arcmins_s4wide_f2, white_noise_T_s4wide_f2, elknee_T_s4wide_f2, alphaknee_T_s4wide_f2, white_noise_P_s4wide_f2, elknee_P_s4wide_f2, alphaknee_P_s4wide_f2 = np.copy(specs_dic_s4wide[freq2])
                    if TP == 'T':
                        white_noise_s4wide_f1, elknee_s4wide_f1, alphaknee_s4wide_f1 = white_noise_T_s4wide_f1, elknee_T_s4wide_f1, alphaknee_T_s4wide_f1
                        white_noise_s4wide_f2, elknee_s4wide_f2, alphaknee_s4wide_f2 = white_noise_T_s4wide_f2, elknee_T_s4wide_f2, alphaknee_T_s4wide_f2
                    elif TP == 'P':
                        white_noise_s4wide_f1, elknee_s4wide_f1, alphaknee_s4wide_f1 = white_noise_P_s4wide_f1, elknee_P_s4wide_f1, alphaknee_P_s4wide_f1
                        white_noise_s4wide_f2, elknee_s4wide_f2, alphaknee_s4wide_f2 = white_noise_P_s4wide_f2, elknee_P_s4wide_f2, alphaknee_P_s4wide_f2

                    if (0): #noise scaling based on total_obs_time
                        noise_scaling_fac = (total_obs_time_default / total_obs_time)**0.5
                        white_noise_s4wide_f1 = white_noise_s4wide_f1 * noise_scaling_fac
                        white_noise_s4wide_f2 = white_noise_s4wide_f2 * noise_scaling_fac

                    #print('s4wide', white_noise_s4wide_f1, elknee_s4wide_f1, alphaknee_s4wide_f1, white_noise_s4wide_f2, elknee_s4wide_f2, alphaknee_s4wide_f2)
                    #print('s4deep', noiseval1, elknee1, alphaknee1, noiseval2, elknee2, alphaknee2)
                    if freq1 == freq2:
                        nl_s4_wide = misc.get_nl(white_noise_s4wide_f1, el, beam_arcmins_s4wide_f1, elknee=elknee_s4wide_f1, alphaknee=alphaknee_s4wide_f1)
                    else:
                        if freq2 in specs_s4wide['corr_noise_bands'][freq1]:
                            nl_s4_wide = misc.get_nl(white_noise_s4wide_f1, el, beam_arcmins_s4wide_f1, elknee=elknee_s4wide_f1, alphaknee=alphaknee_s4wide_f1, beamval2=beam_arcmins_s4wide_f2, noiseval2=white_noise_s4wide_f2, elknee2=elknee_s4wide_f2, alphaknee2=alphaknee_s4wide_f2, rho=specs_s4wide['rho'])
                        else:
                            nl_s4_wide = np.zeros( len(el) )
                    nl_s4_wide[el <= param_dict['lmin']] = 0.

                    #perform inverse variance combination now
                    nl_s4_ultradeep = np.copy(nl)
                    nl = 1. / ( (1./nl) + (1./nl_s4_wide) )

                    if (0):#freq1 == 93:# and freq2 == 145:
                        import matplotlib.pyplot as plt  #not imported at module scope; see setup_matplotlib
                        plt.loglog(el, nl_s4_ultradeep, color='black', label='S4-Ultra deep'); plt.loglog(el, nl_s4_wide, color='red', label='S4-Wide');
                        plt.loglog(el, nl, color='darkgreen', label='S4');
                        plt.title('%s: (%s, %s)' % (TP, freq1, freq2)); plt.legend(loc=3); plt.show()

                nl_dic[TP][(freq1, freq2)] = nl

    return nl_dic


# Diagnostic plots

def setup_matplotlib(interactive_mode=1):
    """
    Select the matplotlib backend and set the plotting defaults.

    Parameters
    ----------
    interactive_mode : int, optional
        When zero, force the non-interactive ``Agg`` backend. Default is 1.

    Notes
    -----
    Called by :func:`plot_beams` and :func:`plot_noise_spectra`, and by nothing else.
    The rest of this module never draws anything, so matplotlib is imported here rather than at module scope and a run that does not plot never imports it at all.

    ``mpl.use`` runs before :mod:`matplotlib.pyplot` is imported, which is what makes the backend selection take effect on the first call.
    """

    import matplotlib as mpl
    if not interactive_mode:
        mpl.use('Agg')
    mpl.rcParams['figure.dpi'] = 150
    mpl.rcParams['font.family'] = 'serif'
    mpl.rcParams["figure.facecolor"] = 'white'


def plot_beams(freqarr, bl_dic, interactive_mode=1):
    r"""
    Plot the beam transfer functions.

    Parameters
    ----------
    freqarr : list of int
        Frequency bands in GHz.
    bl_dic : dict
        Band mapped to :math:`B_\ell`, from :func:`misc.get_beam_dic`.
    interactive_mode : int, optional
        Passed to :func:`setup_matplotlib`. Default is 1.

    Notes
    -----
    Diagnostic only; nothing in this module calls it.
    """

    setup_matplotlib(interactive_mode)
    import matplotlib.pyplot as plt
    for freq in freqarr:
        plt.plot(bl_dic[freq], label=freq)
    plt.legend(loc=1)


def plot_noise_spectra(expname, freqarr, el, nl_dic, bl_dic, beam_noise_dic, use_dls=True, beam_decon=True, interactive_mode=1):
    r"""
    Plot the temperature noise spectra band by band.

    Parameters
    ----------
    expname : str
        Experiment configuration name, used for the plot title.
    freqarr : list of int
        Frequency bands in GHz.
    el : array_like
        Multipoles.
    nl_dic : dict
        Noise spectra from :func:`get_noise_spectra`.
    bl_dic : dict
        Band mapped to :math:`B_\ell`.
    beam_noise_dic : dict
        Beams and white-noise levels, printed alongside the measured level.
    use_dls : bool, optional
        Plot :math:`\ell(\ell+1) N_\ell / 2\pi` rather than :math:`N_\ell`. Default is ``True``.
    beam_decon : bool, optional
        Undo the beam weighting before plotting. Default is ``True``.
    interactive_mode : int, optional
        Passed to :func:`setup_matplotlib`. The figure is only shown when this is nonzero. Default is 1.

    Notes
    -----
    Diagnostic only; nothing in this module calls it.
    The white-noise level quoted in the legend is the median of :math:`N_\ell` over :math:`3000 \le \ell \le 5000`.
    """

    setup_matplotlib(interactive_mode)
    import matplotlib.pyplot as plt
    color_arr = ['navy', 'blue', 'teal', 'darkgreen', 'olive', 'goldenrod', 'orangered', 'darkred', 'maroon']
    plt.subplot(111, yscale='log')
    for fcntr, freq in enumerate( freqarr ):
        currnl = nl_dic['T'][(freq, freq)]*bl_dic[freq]**2.
        #currnl = nl_dic['P'][(freq, freq)]*bl_dic[freq]**2.
        tmpinds = np.where( (el >= 3000) & (el <= 5000) )[0]
        meannl = np.median(currnl[tmpinds])
        noise_uk_arcmin = np.sqrt(meannl)/np.radians(1./60.)
        print(noise_uk_arcmin, beam_noise_dic['T'][freq])

        if beam_decon:
            currnl = currnl / bl_dic[freq]**2.

        if use_dls:
            dl_fac = el * (el+1)/2/np.pi
            plt.plot(el, dl_fac * currnl, label='%s GHz ' % (freq), ls='-', color=color_arr[fcntr % len(color_arr)])
        else:
            plt.plot(el, currnl, label=r'%s GHz (%.2f $\mu$K-arcmin)' % (freq, noise_uk_arcmin), ls='-', color=color_arr[fcntr % len(color_arr)])
        #plt.plot(nl_dic_actual['T'][(freq, freq)], lw=2., color=colordic[freq])
    plt.legend(loc=1)
    if not use_dls:
        plt.xlim(0, 5000); plt.ylim(1e-8, 1.)
        plt.ylabel(r'N$_{\ell}$ [$\mu$K$^{2}$]', fontsize=14)
    else:
        plt.ylabel(r'$\ell(\ell+1)/(2\pi)$ N$_{\ell}$ [$\mu$K$^{2}$]', fontsize=14)
        plt.xlim(0, 5000); plt.ylim(.1, 1e5)
    plt.xlabel(r'Multipole $\ell$', fontsize=14)
    expname_str = expname.replace('spt3g_', 'SPT-3G: ').replace('summer', 'Summer').replace('_', r'\_')
    ##expname_str = 'S4-Wide'
    plt.title(r'%s' % (expname_str), fontsize=14)
    ##plt.savefig('s4_wide_nl.png', dpi= 200.); sys.exit()
    #was 'plt.show(); sys.exit()'. sys.exit() in a function would kill a library caller's process, so the caller decides whether to stop.
    if interactive_mode:
        plt.show()


# Covariance and ILC

def get_ignore_fg(param_dict, final_comp):
    """
    Build the list of components to leave out of the covariance.

    Parameters
    ----------
    param_dict : dict
        Parameter dictionary, optionally holding ``ignore_fg``.
    final_comp : str
        Component the ILC preserves, which need not enter the covariance.

    Returns
    -------
    list of str
        Components to ignore.

    Notes
    -----
    :func:`misc.get_param_dict` returns scalars, so an ``ignore_fg`` set in ``params.ini`` arrives as a string rather than a list.
    It is coerced and copied here so that the list stored in ``param_dict``, which is written to the product file, is not modified.
    """

    ignore_fg = param_dict.get('ignore_fg', [])
    if isinstance(ignore_fg, str):
        ignore_fg = [ignore_fg]
    ignore_fg = list(ignore_fg)

    ignore_fg.append(final_comp.lower()) #the required component need not go into the covariance matrix.
    ignore_fg.append('tsz_cib')

    return ignore_fg


def get_covariances(
        param_dict,
        freqarr,
        el,
        nl_dic,
        bl_dic,
        ignore_fg,
        include_gal,
        which_spec_arr,
        reduce_cib_power=None,
        cl_multiplier_dic=None
        ):
    r"""
    Build the analytic signal, foreground and noise covariance.

    Parameters
    ----------
    param_dict : dict
        Parameter dictionary from :func:`misc.get_param_dict`.
    freqarr : list of int
        Frequency bands in GHz.
    el : array_like
        Multipoles.
    nl_dic : dict
        Noise spectra from :func:`get_noise_spectra`.
    bl_dic : dict
        Band mapped to :math:`B_\ell`.
    ignore_fg : list of str
        Components to leave out, from :func:`get_ignore_fg`.
    include_gal : int
        Whether to include galactic dust and synchrotron.
    which_spec_arr : list of str
        Spectra to build, e.g. ``['TT', 'EE']``.
    reduce_cib_power : float, optional
        CIB power reduction at 150 GHz after point source removal. Default is ``None``.
    cl_multiplier_dic : dict, optional
        Per-component spectrum multipliers for a partial ILC. Default is ``None``.

    Returns
    -------
    el : ndarray
        Multipoles as returned by :func:`ilc.get_analytic_covariance`, which need not be the input array.
    cl_dic : dict
        Spectrum label mapped to the band-band covariance.
    fg_cl_dic : dict
        Spectrum label mapped to the per-component spectra.
    """

    cl_dic = {}
    fg_cl_dic = {}
    for which_spec in which_spec_arr:
        if which_spec == 'TT':
            el, cl_dic[which_spec], fg_cl_dic[which_spec] = ilc.get_analytic_covariance(param_dict, freqarr, el=el, nl_dic=nl_dic['T'], ignore_fg=ignore_fg, include_gal=include_gal, bl_dic=bl_dic, reduce_cib_power=reduce_cib_power, cl_multiplier_dic=cl_multiplier_dic)
        else:
            el, cl_dic[which_spec], fg_cl_dic[which_spec] = ilc.get_analytic_covariance(param_dict, freqarr, el=el, nl_dic=nl_dic['P'], ignore_fg=ignore_fg, which_spec=which_spec, pol_frac_per_cent_dust=param_dict['pol_frac_per_cent_dust'], pol_frac_per_cent_radio=param_dict['pol_frac_per_cent_radio'], pol_frac_per_cent_tsz=param_dict['pol_frac_per_cent_tsz'], pol_frac_per_cent_ksz=param_dict['pol_frac_per_cent_ksz'], include_gal=include_gal, bl_dic=bl_dic, reduce_cib_power=reduce_cib_power)
        ###print(which_spec, cl_dic[which_spec])

    #op_dic_for_fg_noise = {}
    #op_dic_for_fg_noise['fg_cl_dic'] = fg_cl_dic
    #op_dic_for_fg_noise['nl_dic'] = nl_dic
    ##np.save('data/cmbs4_nl_and_fgcl_dict.npy', op_dic_for_fg_noise); ##sys.exit()

    return el, cl_dic, fg_cl_dic


def get_residual_power(param_dict, freqarr, el, cl_dic, nc, which_spec_arr, final_comp='cmb', freqcalib_fac=None, null_comp=None):
    r"""
    Run the ILC and return the residual power and weights.

    Parameters
    ----------
    param_dict : dict
        Parameter dictionary from :func:`misc.get_param_dict`.
    freqarr : list of int
        Frequency bands in GHz.
    el : array_like
        Multipoles.
    cl_dic : dict
        Band-band covariance from :func:`get_covariances`.
    nc : int
        Number of bands.
    which_spec_arr : list of str
        Spectra to solve for.
    final_comp : str, optional
        Component the ILC preserves. Default is ``'cmb'``.
    freqcalib_fac : array_like, optional
        Multiplicative gain of each band, in the same order as ``freqarr``, by which the assumed response is scaled.
        Default is ``None``, which is equivalent to unity in every band.
    null_comp : str or list of str, optional
        Component to null with a constrained ILC. None implies that ``TT`` and ``EE`` are solved jointly. Default is ``None``.

    Returns
    -------
    cl_residual : dict
        Spectrum label mapped to the residual :math:`C_\ell` after the ILC.
    weights_dic : dict
        Spectrum label mapped to the per-band ILC weights.
    """

    cl_residual, weights_dic = {}, {}
    if null_comp is None:
        cl_residual_arr, weights_arr = ilc.residual_power(param_dict, freqarr, el, cl_dic, final_comp=final_comp, freqcalib_fac=freqcalib_fac, return_weights=1)
        single_spec = len([ws for ws in which_spec_arr if ws in ('TT', 'EE')]) == 1
        for which_spec in which_spec_arr:
            if single_spec:
                cl_res = cl_residual_arr[0]
                weights = weights_arr[:nc, 0]
            elif which_spec == 'TT':
                cl_res = cl_residual_arr[0]
                weights = weights_arr[:nc, 0]
            elif which_spec == 'EE':
                cl_res = cl_residual_arr[1]
                weights = weights_arr[nc:, 1]
            elif which_spec == 'TE':
                cl_res = cl_residual_arr[2]
                weights = np.asarray( [weights_arr[nc:, 0], weights_arr[:nc, 1]] )
            cl_residual[which_spec], weights_dic[which_spec] = cl_res, weights
    else:
        cl_dic_TT, cl_dic_EE = {}, {}
        cl_dic_TT['TT'] = cl_dic['TT']
        cl_dic_EE['EE'] = cl_dic['EE']
        cl_residual_TT_arr, weights_TT_arr = ilc.residual_power(param_dict, freqarr, el, cl_dic_TT, final_comp=final_comp, freqcalib_fac=freqcalib_fac, return_weights=1, null_comp=null_comp)
        cl_residual_EE_arr, weights_EE_arr = ilc.residual_power(param_dict, freqarr, el, cl_dic_EE, final_comp=final_comp, freqcalib_fac=freqcalib_fac, return_weights=1, null_comp=null_comp)
        cl_residual['TT'], weights_dic['TT'] = cl_residual_TT_arr[0], weights_TT_arr[:nc, 0]
        cl_residual['EE'], weights_dic['EE'] = cl_residual_EE_arr[0], weights_EE_arr[:nc, 0]

    return cl_residual, weights_dic


def get_fg_residuals(freqarr, el, fg_cl_dic, weights_dic, include_gal, which_spec_arr=None, debug=0):
    r"""
    Propagate each foreground component through the ILC weights.

    Parameters
    ----------
    freqarr : list of int
        Frequency bands in GHz.
    el : array_like
        Multipoles.
    fg_cl_dic : dict
        Per-component spectra from :func:`get_covariances`.
    weights_dic : dict
        ILC weights from :func:`get_residual_power`.
    include_gal : int
        Whether galactic dust and synchrotron are included in the component list.
    which_spec_arr : list of str, optional
        Spectra to propagate. ``'TE'`` is skipped. Default is ``None``, which uses :data:`WHICH_SPEC_ARR`.
    debug : int, optional
        When nonzero, print progress every 2500 multipoles. Default is 0.

    Returns
    -------
    dict
        Spectrum label mapped to component name mapped to the residual power, as

        .. math::

            C_\ell^\mathrm{res} = \mathbf{w}_\ell^\mathsf{T}
                                  \mathbf{C}_\ell \mathbf{w}_\ell\, .

    Warns
    -----
    UserWarning
        If ``which_spec_arr`` contains ``'TE'``, whose weights are two dimensional per multipole.
    """

    fg_res_dic = {}
    #signal_arr = ['galdust', 'galsync']
    if include_gal:
        signal_arr = ['tsz', 'cib', 'radio', 'galdust', 'galsync', 'noise']#, 'tsz-cib']
    else:
        signal_arr = ['tsz', 'cib', 'radio', 'noise']#, 'tsz-cib']
    #was hardcoded to ['TT', 'EE']
    if which_spec_arr is None:
        which_spec_arr = WHICH_SPEC_ARR
    if 'TE' in which_spec_arr:
        warnings.warn("'TE' is skipped", stacklevel=2)
    for which_spec in [ws for ws in which_spec_arr if ws in ('TT', 'EE')]:
        fg_res_dic[which_spec] = {}
        for elcnt, currel in enumerate(el):
            if debug and (elcnt % 2500) == 0:
                print(which_spec, elcnt)
                print('\n')
            for s in signal_arr:
                #print(s)
                if s == 'galdust':
                    #curr_cl_dic = cl_gal_dust_dic[which_spec]
                    curr_cl_dic = fg_cl_dic[which_spec]['galdust']
                elif s == 'galsync':
                    #curr_cl_dic = cl_gal_sync_dic[which_spec]
                    curr_cl_dic = fg_cl_dic[which_spec]['galsync']
                elif s == 'ksz':
                    #curr_cl_dic = cl_ksz_dic[which_spec]
                    curr_cl_dic = fg_cl_dic[which_spec]['ksz']
                elif s == 'cmb':
                    #curr_cl_dic = cl_cmb_dic[which_spec]
                    curr_cl_dic = fg_cl_dic[which_spec]['cmb']
                elif s == 'tsz':
                    #curr_cl_dic = cl_tsz_dic[which_spec]
                    curr_cl_dic = fg_cl_dic[which_spec]['tsz']
                elif s == 'cib':
                    #curr_cl_dic = cl_dust_dic[which_spec]
                    curr_cl_dic = fg_cl_dic[which_spec]['cib']
                elif s == 'radio':
                    #curr_cl_dic = cl_radio_dic[which_spec]
                    curr_cl_dic = fg_cl_dic[which_spec]['radio']
                elif s == 'tsz-cib':
                    #curr_cl_dic = cl_tsz_cib_dic[which_spec]
                    curr_cl_dic = fg_cl_dic[which_spec]['tsz_cib']
                elif s == 'noise':
                    curr_cl_dic = fg_cl_dic[which_spec]['noise']

                tmp_cl_dic = {which_spec: curr_cl_dic}
                clmat = ilc.create_clmat(freqarr, elcnt, tmp_cl_dic)
                currw_ilc = weights_dic[which_spec][:, elcnt]

                curr_res_ilc = float(np.dot(currw_ilc, np.dot(clmat, currw_ilc)))
                if s not in fg_res_dic[which_spec]:
                    fg_res_dic[which_spec][s] = []
                fg_res_dic[which_spec][s].append( curr_res_ilc )

        for s in signal_arr:
            fg_res_dic[which_spec][s] = np.asarray(fg_res_dic[which_spec][s])

    return fg_res_dic


# Output

def get_fsky_val(param_dict):
    """
    Read the unmasked sky fraction from the galactic simulation products.

    Parameters
    ----------
    param_dict : dict
        Parameter dictionary holding ``cl_gal_dic_dust_fname``, optionally ``cl_gal_folder``, and ``which_gal_mask``.

    Returns
    -------
    float
        Sky fraction for the selected galactic mask.
    """

    cl_gal_dic_dust_fname = param_dict['cl_gal_dic_dust_fname']
    try:
        cl_gal_folder = param_dict['cl_gal_folder']
        cl_gal_dic_dust_fname = '%s/%s' % (cl_gal_folder, cl_gal_dic_dust_fname)
    except KeyError:
        pass
    galdustsims_cl = np.load(cl_gal_dic_dust_fname, allow_pickle=1, encoding='latin1').item()
    #used to be computed inline inside an 'if include_gal:' block that held its own 'if not include_gal' branch, so the 0.68 fallback was unreachable. The include_gal guard now sits at the call site in build_output_dic.
    #if not include_gal:
    #    fsky_val = 0.68
    #else:
    #    fsky_val = galdustsims_cl['fsky_arr'][param_dict['which_gal_mask']]

    return galdustsims_cl['fsky_arr'][param_dict['which_gal_mask']]


def build_output_dic(
        el,
        cl_residual,
        fg_res_dic,
        weights_dic,
        beam_noise_dic,
        elknee_dic,
        param_dict,
        which_gal_mask,
        include_gal,
        save_fg_res_and_weights=1,
        freqcalib_fac=None,
        cl_multiplier_dic=None,
        noise_scalings_for_bands=None
        ):
    """
    Assemble the dictionary provided and/or written to the product file.

    Parameters
    ----------
    el : array_like
        Multipoles.
    cl_residual : dict
        Residual power from :func:`get_residual_power`.
    fg_res_dic : dict
        Per-component residuals from :func:`get_fg_residuals`.
    weights_dic : dict
        ILC weights from :func:`get_residual_power`.
    beam_noise_dic : dict
        Beams and white-noise levels.
    elknee_dic : dict
        Knee multipoles and slopes.
    param_dict : dict
        Parameter dictionary, stored so products are self describing.
    which_gal_mask : int
        Galactic mask index.
    include_gal : int
        Whether galactic foregrounds were included.
    save_fg_res_and_weights : int, optional
        Whether to store the per-component residuals and the weights. Default is 1.
    freqcalib_fac : array_like, optional
        Multiplicative gain of each band, in the same order as ``freqarr``, by which the assumed response is scaled.
        Default is ``None``, which is equivalent to unity in every band.
    cl_multiplier_dic : dict, optional
        Per-component spectrum multipliers. Default is ``None``.
    noise_scalings_for_bands : list of float, optional
        Per-band noise scalings. Default is ``None``.

    Returns
    -------
    dict
        Contents of the product file.

    Warns
    -----
    UserWarning
        When ``noise_scalings_for_bands`` is given, since that forces ``save_fg_res_and_weights`` to zero and so overrides an explicit request.
    """

    #freq0, lmax = param_dict['freq0'], param_dict['lmax']  #unused

    #save residual files
    if noise_scalings_for_bands is not None:
        if save_fg_res_and_weights:
            warnings.warn('noise_scalings_for_bands is set, so save_fg_res_and_weights is forced to 0 and neither fg_res_dic nor weights will be stored', stacklevel=2)
        save_fg_res_and_weights = 0

    opdic = {}
    opdic['el'] = el
    if cl_multiplier_dic is not None:
        opdic['cl_multiplier_dic'] = cl_multiplier_dic
    opdic['cl_residual'] = cl_residual
    if save_fg_res_and_weights:
        opdic['fg_res_dic'] = fg_res_dic
        opdic['weights'] = weights_dic
    opdic['freqcalib_fac'] = freqcalib_fac
    opdic['param_dict'] = param_dict
    if include_gal:
        opdic['fsky_val'] = get_fsky_val(param_dict)
    opdic['which_gal_mask'] = which_gal_mask
    opdic['beam_noise_dic'] = beam_noise_dic
    opdic['elknee_dic'] = elknee_dic

    return opdic


# Driver

def run_ilc(
        expname,
        total_obs_time=7.0,
        include_gal=0,
        which_gal_mask=2,
        save_fg_res_and_weights=1,
        s4_so_joint_configs=0,
        include_fulls4scaledsobaseline=0,
        noise_scalings_for_bands=None,
        final_comp='cmb',
        debug=0,
        paramfile=PARAMFILE,
        which_spec_arr=None,
        null_comp=None,
        cl_multiplier_dic=None,
        freqcalib_fac=None,
        save=True
        ):
    """
    Compute the ILC residuals for one experiment configuration.

    Parameters
    ----------
    expname : str
        Experiment configuration name, passed to :func:`get_experiment_specs`.
    total_obs_time : float, optional
        Observing time in years. Default is 7.0.
    include_gal : int, optional
        Include galactic dust and synchrotron. Default is 0.
    which_gal_mask : int, optional
        Galactic mask index. Forced to -1 when ``include_gal`` is zero. Default is 2.
    save_fg_res_and_weights : int, optional
        Store the per-component residuals and the ILC weights. Default is 1.
    s4_so_joint_configs : int, optional
        Mark this as a joint S4 and SO configuration in the output path. Default is 0.
    include_fulls4scaledsobaseline : int, optional
        Combine in the scaled SO-baseline noise. Default is 0.
    noise_scalings_for_bands : list of float, optional
        Multiplicative white-noise scaling per band. Default is ``None``.
    final_comp : str, optional
        Component the ILC preserves, one of :data:`COMP_CHOICES`. Default is ``'cmb'``.
    debug : int, optional
        Print intermediate dictionary keys and progress. Default is 0.
    paramfile : str, optional
        Parameter file to read. Default is :data:`PARAMFILE`.
    which_spec_arr : list of str, optional
        Spectra to solve for. Default is :data:`WHICH_SPEC_ARR`.
    null_comp : str or list of str, optional
        Component or components to null with a constrained ILC, each one of :data:`COMP_CHOICES`. Default is ``None``.
    cl_multiplier_dic : dict, optional
        Per-component spectrum multipliers for a partial ILC. Default is ``None``, which is treated as an empty dictionary.
    freqcalib_fac : array_like, optional
        Multiplicative gain of each band, in the same order as ``freqarr``, by which the assumed response is scaled.
        Default is ``None``, which is equivalent to unity in every band.
    save : bool, optional
        Write the result to the product file. Default is ``True``.

    Returns
    -------
    opdic : dict
        Contents of the product file, from :func:`build_output_dic`.
    opfname : str
        Path the result was or would have been written to.

    Raises
    ------
    ValueError
        If ``final_comp`` is not in :data:`COMP_CHOICES`, if a band or mask lacks galactic simulations, if ``expname`` and ``total_obs_time`` imply a double observing-time scaling, or if ``noise_scalings_for_bands`` has the wrong length.
    """

    if final_comp not in COMP_CHOICES:
        raise ValueError('-final_comp must be one of %s, got %s' % (COMP_CHOICES, final_comp))
    if null_comp is not None and np.ndim(null_comp) == 0:
        null_comp = [null_comp]
    if which_spec_arr is None:
        which_spec_arr = list(WHICH_SPEC_ARR)
    if cl_multiplier_dic is None:
        #cl multipler - multiply a given spectra by some amount to perform partial ILC. similar to https://arxiv.org/abs/2102.05033
        cl_multiplier_dic = {}

    # read and store param dict
    param_dict = misc.get_param_dict(paramfile)

    if not os.path.exists(param_dict['data_folder']):
        #param_dict['data_folder'] = '/data/spt/sri-data48/git/DRAFT/data/'
        param_dict['data_folder'] = os.path.join( os.path.dirname(os.path.abspath(__file__)), 'data' )
    if 'cl_gal_folder' in param_dict and not os.path.exists(param_dict['cl_gal_folder']):
        param_dict['cl_gal_folder'] = os.path.join( os.path.dirname(os.path.abspath(__file__)), param_dict['cl_gal_folder'] )

    el = np.arange(param_dict['lmax'])

    #20220112 - moved to argparse
    ###include_gal = param_dict['include_gal'] ##1
    param_dict['include_gal'] = include_gal
    if not include_gal:
        which_gal_mask = -1
    param_dict['which_gal_mask'] = which_gal_mask

    try:
        remove_atm = param_dict['remove_atm']
    except KeyError:
        remove_atm = 0

    specs_dic, corr_noise_bands, rho, corr_noise, Nred_dic, secondary_specs = get_experiment_specs(expname, remove_atm=remove_atm, include_fulls4scaledsobaseline=include_fulls4scaledsobaseline)

    freqarr = sorted(specs_dic.keys())
    nc = len(freqarr)

    if include_gal:
        check_gal_sim_inputs(freqarr, which_gal_mask, expname)

    reduce_cib_power = None
    if expname.find('cmbhd') > -1:
        reduce_cib_power = 17. #150 GHz power reduction after removing sources above 0.04 mJy

    total_obs_time_default = get_total_obs_time_default(expname, total_obs_time)

    opfname = build_output_path(expname, param_dict, freqarr, which_spec_arr, include_gal, which_gal_mask, final_comp, total_obs_time, corr_noise, remove_atm, noise_scalings_for_bands=noise_scalings_for_bands, s4_so_joint_configs=s4_so_joint_configs, null_comp=null_comp, cl_multiplier_dic=cl_multiplier_dic)
    if save:
        opfolder = '/'.join(opfname.split('/')[:-1])
        #if not os.path.exists(opfolder): os.system('mkdir -p %s' % (opfolder))
        os.makedirs(opfolder, exist_ok=True)

    beam_noise_dic, elknee_dic, noisearr_T, noisearr_P = get_beam_and_noise_dic(expname, specs_dic, freqarr, total_obs_time, total_obs_time_default, noise_scalings_for_bands=noise_scalings_for_bands, specs_aso=secondary_specs.get('aso'), specs_fulls4scaledsobaseline=secondary_specs.get('fulls4scaledsobaseline'))

    print('Delta T =', noisearr_T)
    print('Delta P =', noisearr_P)
    print('')

    nl_dic = get_noise_spectra(freqarr, el, param_dict, beam_noise_dic, elknee_dic, corr_noise_bands, rho, Nred_dic, specs_s4wide=secondary_specs.get('s4wide'), total_obs_time=total_obs_time, total_obs_time_default=total_obs_time_default)
    if debug:
        print(nl_dic['T'].keys())
    ##sys.exit()

    #get beams
    bl_dic = misc.get_beam_dic(freqarr, beam_noise_dic['T'], param_dict['lmax'])
    if debug:
        print(bl_dic.keys())

    #get the CMB, noise, and foreground covriance
    ignore_fg = get_ignore_fg(param_dict, final_comp)
    #freqarr = [145]
    #param_dict['which_gal_mask'] = 0
    if debug:
        print(ignore_fg)

    el, cl_dic, fg_cl_dic = get_covariances(param_dict, freqarr, el, nl_dic, bl_dic, ignore_fg, include_gal, which_spec_arr, reduce_cib_power=reduce_cib_power, cl_multiplier_dic=cl_multiplier_dic)
    if debug:
        print(cl_dic.keys())
        print(el)

    #get the residual power now
    #null_comp = None
    cl_residual, weights_dic = get_residual_power(param_dict, freqarr, el, cl_dic, nc, which_spec_arr, final_comp=final_comp, freqcalib_fac=freqcalib_fac, null_comp=null_comp)
    if debug:
        print(cl_residual.keys())

    fg_res_dic = get_fg_residuals(freqarr, el, fg_cl_dic, weights_dic, include_gal, which_spec_arr=which_spec_arr, debug=debug)
    if debug:
        print(fg_res_dic.keys())

    opdic = build_output_dic(el, cl_residual, fg_res_dic, weights_dic, beam_noise_dic, elknee_dic, param_dict, which_gal_mask, include_gal, save_fg_res_and_weights=save_fg_res_and_weights, freqcalib_fac=freqcalib_fac, cl_multiplier_dic=cl_multiplier_dic, noise_scalings_for_bands=noise_scalings_for_bands)
    #opdic['nl_dic'] = nl_dic

    if save:
        np.save(opfname, opdic)
    print(opfname)

    return opdic, opfname


def main(argv=None):
    """
    Parse the command line and run the ILC.

    Parameters
    ----------
    argv : list of str, optional
        Argument list. Default is ``sys.argv[1:]``.
    """

    args = parse_args(argv)
    #parsed values used to be pushed into module globals with exec()
    #args_keys = args.__dict__
    #for kargs in args_keys:
    #    param_value = args_keys[kargs]
    #    if isinstance(param_value, str):
    #        cmd = '%s = "%s"' % (kargs, param_value)
    #    else:
    #        cmd = '%s = %s' % (kargs, param_value)
    #    exec(cmd)
    run_ilc(args.expname, total_obs_time=args.total_obs_time, include_gal=args.include_gal, which_gal_mask=args.which_gal_mask, save_fg_res_and_weights=args.save_fg_res_and_weights, s4_so_joint_configs=args.s4_so_joint_configs, include_fulls4scaledsobaseline=args.include_fulls4scaledsobaseline, noise_scalings_for_bands=args.noise_scalings_for_bands, final_comp=args.final_comp, null_comp=args.null_comp, paramfile=args.paramfile, debug=args.debug)
    print('\nDone.')


if __name__ == '__main__':
    main()
