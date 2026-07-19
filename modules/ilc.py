import numpy as np, sys, os, scipy as sc, foregrounds as fg, misc, re

################################################################################################################
def get_analytic_covariance(param_dict, 
    freqarr, 
    el = None, 
    nl_dic = None, 
    bl_dic = None, 
    ignore_fg = [], 
    which_spec = 'TT', 
    pol_frac_per_cent_dust = 0.02, 
    pol_frac_per_cent_radio = 0.03, 
    pol_frac_per_cent_tsz = 0., 
    pol_frac_per_cent_ksz = 0., 
    include_gal = 0, 
    cib_corr_coeffs = None, 
    null_highfreq_radio = 1, 
    reduce_radio_power_150 = None, 
    reduce_tsz_power = None, 
    reduce_cib_power = None, 
    remove_cib_decorr = 0, 
    cl_multiplier_dic = None, 
    return_fg_spectra = True, 
    force_cl_dic = None,
    ):

    """
    Get signal + noise covariance for ILC.
    Supports MV-ILC, cILC, partial ILC, etc..

    Parameters
    ----------
    param_dict: dict
        Dictionary containing the params
    freqarr : array
        array of frequency bands for which we need the covariance.
    el : array
        Multipoles over which the covariance must be defined.
        Default is None in which case el = np.arange( len(cl_cmb) ) derived below.
    nl_dic : dict
        Dictionary containing the noise covariance in the bands.
        Default is None.
    bl_dic : dict
        Dictionary containing the beams. Only used for galactic foreground file since cl_gal have S4 beams.
        Default is None.
    ignore_fg : list
        List of signals that must be excluded from the covariance. 
    which_spec : str
        spectra name TT/EE/TE.
    pol_frac_per_cent_dust : float
        Pol fraction for CIB.
        Default is 0.02.
    pol_frac_per_cent_radio : float
        Pol fraction for Radio galaxies.
        Default is 0.03.
    pol_frac_per_cent_tsz : float
        Pol fraction for tSZ. Default is zero.
    pol_frac_per_cent_ksz : float
        Pol fraction for kSZ. Default is zero.
    include_gal : bool
        Include galactic foregrounds.
        Default is None.
    cib_corr_coeffs : dict
        Cross-correlation coefficients across bands for CIB.
        Default is None.
    null_highfreq_radio : bool
        Null radio beyond 230 GHz.
        Deault is True.
    reduce_radio_power_150 : float
        Reduce radio power by some number.
        Default is None.
    reduce_tsz_power : float
        Reduce tSZ power by some number.
        Default is None.
    reduce_cib_power : float
        Reduce CIB power by some number.
        Default is None.
    remove_cib_decorr : bool
        Remove CIB decorrelations.
        Default is False.
    cl_multiplier_dic : dict 
        Dictionary containing factors by which a certain signal must be suppressed in the covariance.
        For partial ILC.
        Default is None.
    return_fg_spectra : bool
        Return the spectra for each foreground signal along with the total covariance.
        Default is True.
    force_cl_dic : dict
        Supply own foreground Cl dict.
        Default is None.

    Returns
    -------
    el : array
        Multipoles over which the covariance is defined.
    cl_dic : dict
        Total covariance in each band as a dictionary.        
    fg_cl_dic: dict
        Spectra for each foreground signal in each band as a dictionary.
        Only returned when return_fg_spectra is True.
    """

    #ignore_fg = foreground terms that must be ignored
    debug=False
    possible_ignore_fg = ['cmb', 'tsz', 'y', 'ksz', 'radio', 'dust', 'noise', 'tsz_cib']
    if len(ignore_fg)>0:
        if 'cmb' in ignore_fg: ignore_fg.append('ksz')
        if not all( [ currfg in possible_ignore_fg for currfg in ignore_fg] ):
            print( '\n\t Alert: Elements of ignore_fg should be one of the following: %s\n\n' %(np.array2string(np.asarray(possible_ignore_fg))) )
            sys.exit()
        ignore_fg = np.unique(ignore_fg)


    el_, cl_cmb = fg.get_foreground_power_spt('CMB', freq1 = param_dict['freq0'], freq2 = param_dict['freq0'])
    if el is None: el = np.copy(el_)    
    cl_cmb = np.interp(el, el_, cl_cmb)
    el_, cl_ksz = fg.get_foreground_power_spt('kSZ', freq1 = param_dict['freq0'], freq2 = param_dict['freq0'])
    cl_ksz = np.interp(el, el_, cl_ksz)

    if which_spec == 'EE':
        cl_ksz = cl_ksz * pol_frac_per_cent_ksz**2.
    if which_spec == 'TE':
        cl_ksz = cl_ksz * 0.

    if cl_multiplier_dic is not None:
        if 'cmb' in cl_multiplier_dic:
            cl_cmb = np.copy(cl_cmb) * cl_multiplier_dic['cmb']
        if 'ksz' in cl_multiplier_dic:
            cl_ksz = np.copy(cl_ksz) * cl_multiplier_dic['ksz']

    #20220428 - force cl if force_cl_dic is supplied
    if force_cl_dic is None: force_cl_dic = {}
    if 'cmb' in force_cl_dic:
        cl_cmb = np.interp(el, np.arange(len(force_cl_dic['cmb'])), force_cl_dic['cmb'])
    if 'ksz' in force_cl_dic:
        cl_ksz = np.interp(el, np.arange(len(force_cl_dic['ksz'])), force_cl_dic['ksz'])

    cl_dic = {}
    cl_ori = np.zeros(len(el))

    if return_fg_spectra:
        fg_cl_dic = {}

    for freq1 in freqarr:
        for freq2 in freqarr:

            #get tsz
            el_, cl_tsz = fg.get_cl_tsz(freq1, freq2, freq0 = param_dict['freq0'], fg_model = param_dict['fg_model'], reduce_tsz_power = reduce_tsz_power)
            if which_spec == 'EE':
                cl_tsz = cl_tsz * pol_frac_per_cent_tsz**2.
            elif which_spec == 'TE':
                cl_tsz = cl_tsz * 0.
            cl_tsz = np.interp(el, el_, cl_tsz)

            #get radio
            el_, cl_radio = fg.get_cl_radio(freq1, freq2, freq0 = param_dict['freq0'], fg_model = param_dict['fg_model'], spec_index_rg = param_dict['spec_index_rg'], null_highfreq_radio = null_highfreq_radio, reduce_radio_power_150 = reduce_radio_power_150)
            if which_spec == 'EE':
                cl_radio = cl_radio * pol_frac_per_cent_radio**2.
            elif which_spec == 'TE':
                cl_radio = cl_radio * 0.
            cl_radio = np.interp(el, el_, cl_radio)

            #get CIB
            el_,  cl_dg_po, cl_dg_clus = fg.get_cl_dust(freq1, freq2, freq0 = param_dict['freq0'], fg_model = param_dict['fg_model'], spec_index_dg_po = param_dict['spec_index_dg_po'], spec_index_dg_clus = param_dict['spec_index_dg_clus'], Tcib = param_dict['Tcib'], reduce_cib_power = reduce_cib_power)
            cl_dust = cl_dg_po + cl_dg_clus
            if remove_cib_decorr:
                el_,  cl_dg_po1, cl_dg_clus1 = fg.get_cl_dust(freq1, freq1, freq0 = param_dict['freq0'], fg_model = param_dict['fg_model'], spec_index_dg_po = param_dict['spec_index_dg_po'], spec_index_dg_clus = param_dict['spec_index_dg_clus'], Tcib = param_dict['Tcib'], reduce_cib_power = reduce_cib_power)
                el_,  cl_dg_po2, cl_dg_clus2 = fg.get_cl_dust(freq2, freq2, freq0 = param_dict['freq0'], fg_model = param_dict['fg_model'], spec_index_dg_po = param_dict['spec_index_dg_po'], spec_index_dg_clus = param_dict['spec_index_dg_clus'], Tcib = param_dict['Tcib'], reduce_cib_power = reduce_cib_power)
                cl_dust1 = cl_dg_po1 + cl_dg_clus1
                cl_dust2 = cl_dg_po2 + cl_dg_clus2
                cl_dust = np.sqrt( cl_dust1 * cl_dust2 )
                cib_corr_coeffs = None

            cl_dust[np.isnan(cl_dust)] = 0.
            if which_spec == 'EE':
                cl_dust = cl_dust * pol_frac_per_cent_dust**2.
            elif which_spec == 'TE':
                cl_dust = cl_dust * 0.

            #include CIB decorrelation if available
            if cib_corr_coeffs is not None:
                if freq1 == freq2:
                    corr_coeff = 1.
                else:
                    if (freq1, freq2) in cib_corr_coeffs:
                        corr_coeff = cib_corr_coeffs[(freq1, freq2)]
                    elif (freq2, freq1) in cib_corr_coeffs:
                        corr_coeff = cib_corr_coeffs[(freq2, freq1)]
                cl_dust *= corr_coeff
            cl_dust = np.interp(el, el_, cl_dust)

            #get tSZ x CIB
            el_, cl_tsz_cib = fg.get_cl_tsz_cib(freq1, freq2, freq0 = param_dict['freq0'], fg_model = param_dict['fg_model'], spec_index_dg_po = param_dict['spec_index_dg_po'], spec_index_dg_clus = param_dict['spec_index_dg_clus'], Tcib = param_dict['Tcib'], reduce_tsz_power = reduce_tsz_power)
            if which_spec == 'EE' or which_spec == 'TE':
                cl_tsz_cib = cl_tsz_cib * 0.
            cl_tsz_cib = np.interp(el, el_, cl_tsz_cib)
            
            #galaxy
            if include_gal:# and not pol: #get galactic dust and sync
                el_, cl_gal_dust = fg.get_cl_galactic(param_dict, 'dust', freq1, freq2, which_spec, el = el, bl_dic = bl_dic)
                el_, cl_gal_sync = fg.get_cl_galactic(param_dict, 'sync', freq1, freq2, which_spec, el = el, bl_dic = bl_dic)

            cl = np.copy( cl_ori )

            #20220428 - force cl if force_cl_dic is supplied
            if 'y' in force_cl_dic:
                cl_y_force = force_cl_dic['y']
                tsz_fac_freq1 = compton_y_to_delta_Tcmb(freq1)
                tsz_fac_freq2 = compton_y_to_delta_Tcmb(freq2)
                tsz_fac = tsz_fac_freq1 * tsz_fac_freq2
                cl_tsz_force = cl_y_force * tsz_fac
                cl_tsz = np.interp(el, np.arange(len(cl_tsz_force)), cl_tsz_force)
            if 'tsz' in force_cl_dic:
                cl_tsz_force = force_cl_dic['tsz'][(freq1, freq2)]
                cl_tsz = np.interp(el, np.arange(len(cl_tsz_force)), cl_tsz_force)
            if 'cib' in force_cl_dic:
                if (freq1, freq2) in force_cl_dic['cib']:
                    cl_dust_force = force_cl_dic['cib'][(freq1, freq2)]
                else:
                    cl_dust_force = force_cl_dic['cib'][(freq2, freq1)]
                cl_dust = np.interp(el, np.arange(len(cl_dust_force)), cl_dust_force)
            if 'tsz_cib' in force_cl_dic or 'cib_tsz' in force_cl_dic:
                if 'tsz_cib' in force_cl_dic:
                    cl_tsz_cib_force = force_cl_dic['tsz_cib'][(freq1, freq2)]
                else:
                    cl_tsz_cib_force = force_cl_dic['cib_tsz'][(freq1, freq2)]
                cl_tsz_cib = np.interp(el, np.arange(len(cl_tsz_cib_force)), cl_tsz_cib_force)
            if 'rad' in force_cl_dic or 'radio' in force_cl_dic:
                if 'rad' in force_cl_dic:
                    cl_radio_force = force_cl_dic['rad'][(freq1, freq2)]
                else:
                    cl_radio_force = force_cl_dic['radio'][(freq1, freq2)]
                cl_radio = np.interp(el, np.arange(len(cl_radio_force)), cl_radio_force)

            if cl_multiplier_dic is not None:
                if 'tsz' in cl_multiplier_dic:
                    cl_tsz = np.copy(cl_tsz) * cl_multiplier_dic['tsz']
                if 'radio' in cl_multiplier_dic:
                    cl_radio = np.copy(cl_radio) * cl_multiplier_dic['radio']
                if 'dust' in cl_multiplier_dic:
                    cl_dust = np.copy(cl_dust) * cl_multiplier_dic['dust']
                if 'tsz_cib' in cl_multiplier_dic:
                    cl_tsz_cib = np.copy(cl_tsz_cib) * cl_multiplier_dic['tsz_cib']
                if include_gal:
                    if 'gal_dust' in cl_multiplier_dic:
                        cl_gal_dust = np.copy(cl_gal_dust) * cl_multiplier_dic['gal_dust']
                    if 'gal_sync' in cl_multiplier_dic:
                        cl_gal_sync = np.copy(cl_gal_sync) * cl_multiplier_dic['gal_sync']

            if 'cmb' not in ignore_fg:
                if len(cl_cmb)<len(cl):
                    cl_cmb = np.interp(el, np.arange(len(cl_cmb)), cl_cmb)
                cl = cl + np.copy(cl_cmb[el])
            if 'ksz' not in ignore_fg:
                if len(cl_ksz)<len(cl):
                    cl_ksz = np.interp(el, np.arange(len(cl_ksz)), cl_ksz)
                cl = cl + cl_ksz[el]
            if 'tsz' not in ignore_fg and 'y' not in ignore_fg:
                if len(cl_tsz)<len(cl):
                    cl_tsz = np.interp(el, np.arange(len(cl_tsz)), cl_tsz)
                cl = cl + cl_tsz[el]
            if 'radio' not in ignore_fg:
                if len(cl_radio)<len(cl):
                    cl_radio = np.interp(el, np.arange(len(cl_radio)), cl_radio)
                cl = cl + cl_radio[el]
            if 'dust' not in ignore_fg:
                if len(cl_dust)<len(cl):
                    cl_dust = np.interp(el, np.arange(len(cl_dust)), cl_dust)
                cl = cl + cl_dust[el]

            #20220503 - add tszxcib if either tsz or cib is included.
            add_cl_tsz_cib = True
            #if ('dust' in ignore_fg and 'tsz' in ignore_fg) or 'tsz_cib' in ignore_fg or 'cib_tsz' in ignore_fg:
            if 'tsz_cib' in ignore_fg or 'cib_tsz' in ignore_fg or 'cib_y' in ignore_fg or 'y_cib' in ignore_fg:
                add_cl_tsz_cib = False 
            if add_cl_tsz_cib: #'dust' not in ignore_fg and 'tsz' not in ignore_fg and 'tsz_cib' not in ignore_fg:
                if len(cl_tsz_cib)<len(cl):
                    cl_tsz_cib = np.interp(el, np.arange(len(cl_tsz_cib)), cl_tsz_cib)
                cl = cl + cl_tsz_cib[el]

            if include_gal:# and not pol: #get galactic dust and sync

                cl = cl + cl_gal_dust
                cl = cl + cl_gal_sync

            #make sure cl start from el=0 rather than el=10 which is the default for SPT G15 results
            lmin = min(el)
            cl = np.concatenate( (np.zeros(lmin), cl) )

            #noise auto power spectrum
            if nl_dic is not None:
                if (freq1, freq2) in nl_dic:
                    nl = nl_dic[(freq1, freq2)]
                else:
                    nl = nl_dic[freq1]
                    if freq1 != freq2: 
                        nl = np.copy(nl) * 0.

                if len(cl) > len(nl):
                    cl = cl[:len(nl)]
                elif len(cl) < len(nl):
                    nl = nl[:len(cl)]

                #remove very large numbers because of beam deconvolution                
                ini_nl = np.median(nl[:100])
                end_nl = np.median(nl[-100:])
                if end_nl>ini_nl: #this implies beam deconvolution has made end nl pretty large
                    max_nl_value = 5e4 #some large number
                    #having end_nl pretty large introduces covariance inversion issues                        
                    badinds = np.where(nl>=max_nl_value)[0]
                    nl[badinds] = max_nl_value
                    #print(ini_nl, end_nl)

                el = np.arange(len(cl))

            else:
                nl = np.zeros(len(cl))

            if cl_multiplier_dic is not None:
                if 'noise' in cl_multiplier_dic:
                    nl = np.copy(nl) * cl_multiplier_dic['noise']

            if 'noise' not in ignore_fg:
                if which_spec != 'TE': cl = cl + np.copy(nl)

            if return_fg_spectra:
                if 'cmb' not in fg_cl_dic: fg_cl_dic['cmb'] = {}
                fg_cl_dic['cmb'][(freq1, freq2)] = fg_cl_dic['cmb'][(freq2, freq1)] = cl_cmb
                if 'ksz' not in fg_cl_dic: fg_cl_dic['ksz'] = {}
                fg_cl_dic['ksz'][(freq1, freq2)] = fg_cl_dic['ksz'][(freq2, freq1)] = cl_ksz
                if 'tsz' not in fg_cl_dic: fg_cl_dic['tsz'] = {}
                fg_cl_dic['tsz'][(freq1, freq2)] = fg_cl_dic['tsz'][(freq2, freq1)] = cl_tsz
                if 'radio' not in fg_cl_dic: fg_cl_dic['radio'] = {}
                fg_cl_dic['radio'][(freq1, freq2)] = fg_cl_dic['radio'][(freq2, freq1)] = cl_radio
                if 'cib' not in fg_cl_dic: fg_cl_dic['cib'] = {}
                fg_cl_dic['cib'][(freq1, freq2)] = fg_cl_dic['cib'][(freq2, freq1)] = cl_dust
                if 'tsz_cib' not in fg_cl_dic: fg_cl_dic['tsz_cib'] = {}
                fg_cl_dic['tsz_cib'][(freq1, freq2)] = cl_tsz_cib
                if 'noise' not in fg_cl_dic: fg_cl_dic['noise'] = {}
                fg_cl_dic['noise'][(freq1, freq2)] = fg_cl_dic['noise'][(freq2, freq1)] = nl
                if include_gal:
                    if 'galdust' not in fg_cl_dic: 
                        fg_cl_dic['galdust'] = {}
                        fg_cl_dic['galsync'] = {}
                    fg_cl_dic['galdust'][(freq1, freq2)] = fg_cl_dic['galdust'][(freq2, freq1)] = cl_gal_dust
                    fg_cl_dic['galsync'][(freq1, freq2)] = fg_cl_dic['galsync'][(freq2, freq1)] = cl_gal_sync

            cl[np.isnan(cl)] = 0.
            cl[np.isinf(cl)] = 0.

            ##########################################################################################
            #20200516 - adjusting Nl when beam is too large (for 30/40 GHz bands)
            adjust_for_large_beams = False                 
            if adjust_for_large_beams:
                beam_tol_for_ilc = 1000. #some large number
                bl = bl_dic[freq1]
                if 'effective' in bl_dic:
                    bl_eff = bl_dic['effective']
                else:
                    bl_eff = bl_dic[145]
                bad_inds = np.where( (bl_eff/bl > beam_tol_for_ilc) )[0]
                print(freq1, freq2, bad_inds)
                cl[bad_inds] = 0.
            ##########################################################################################

            cl_dic[(freq1, freq2)] = cl
    if return_fg_spectra:
        return el, cl_dic, fg_cl_dic
    else:
        return el, cl_dic  

def get_acap(freqarr, final_comp = 'cmb', freqcalib_fac = None, nspecs = 1):

    """
    get frequency dependence of a sky signal.

    Parameters
    ----------
    freqarr : array
        array of frequency bands for which we need the covariance.

    final_comp : str
        It can be
        'cmb' for CMB spectra
        'tsz' or 'y' for thermal SZ
        'radio' for radio galaxies (see get_radio_freq_dep() function)
        'cib' or 'cibpo' for Poisson component of the CIB.
        'cibclus' for clustered component of the CIB.
        'misc_cib_tcibxx_betayy' for misc CIB with T_d=xx and beta=yy.
        'misc_radio_alphaxx' for misc radio with spectral index=xx.
        In the above case, the code will get the freq dep of radio as
        get_radio_freq_dep(...., alpha = xxx)

    freqcalib_fac: array
        array containing calibration factors / mis-matches between different bands.
        default is None which is equivalent to [1, 1, 1, ... 1] in all bands.

    nspecs : int
        tells if we are performing ILC for T alone or T/E/B together.
        default is 1. For only one map component.

    Returns
    -------
    acap : array
        freq. dependene of the respective sky component.
        for example: CMB will be  [1., 1., ...., 1.] in all bands.
    """
    nc = len(freqarr)

    if freqcalib_fac is None: freqcalib_fac = np.ones(nc)

    if final_comp.lower() == 'cmb':
        freqscale_fac = np.ones(nc)

    elif final_comp.lower() == 'tsz' or final_comp.lower() == 'y':

        freqscale_fac = []
        for freq in sorted( freqarr ):
            freqscale_fac.append( compton_y_to_delta_Tcmb(freq) )

        freqscale_fac = np.asarray( freqscale_fac )

    elif final_comp.lower() == 'tsz_cib' or final_comp.lower() == 'cib_tsz':
        pass

    elif final_comp.lower() == 'cib' or final_comp.lower() == 'cibpo':
        freqscale_fac = []
        for freq in sorted( freqarr ):
            freqscale_fac.append( get_cib_freq_dep(freq) )

        freqscale_fac = np.asarray( freqscale_fac )
        freqscale_fac /= np.max(freqscale_fac)
        
    elif final_comp.lower() == 'cibclus':
        freqscale_fac = []
        for freq in sorted( freqarr ):
            freqscale_fac.append( get_cib_freq_dep(freq, beta = 2.505) )

        freqscale_fac = np.asarray( freqscale_fac )
        freqscale_fac /= np.max(freqscale_fac)

    elif final_comp.lower().find('misc_cib')>-1:
        #default values
        misc_tcib = 20.
        misc_beta = 1.505
        tcib_tmp =  re.findall('tcib\d*\.?\d+', final_comp.lower())
        if len(tcib_tmp)>0:
            tcib_tmp = tcib_tmp[0]
            misc_tcib = float(tcib_tmp.replace('tcib', ''))

        beta_tmp =  re.findall('beta\d*\.?\d+', final_comp.lower())
        if len(beta_tmp)>0:
            beta_tmp = beta_tmp[0]
            misc_beta = float(beta_tmp.replace('beta', ''))

        #freqarr = [30, 44, 70, 100, 150, 217, 353, 545]
        freqscale_fac = []
        for freq in sorted( freqarr ):
            freqscale_fac.append( get_cib_freq_dep(freq, Tcib = misc_tcib, beta = misc_beta) )

        freqscale_fac = np.asarray( freqscale_fac )
        freqscale_fac /= np.max(freqscale_fac)

    elif final_comp.lower() == 'radio':
        freqscale_fac = []
        for freq in sorted( freqarr ):
            freqscale_fac.append( get_radio_freq_dep(freq, spec_index_rg = spec_index_rg) )

        freqscale_fac = np.asarray( freqscale_fac )

    acap = np.zeros(nc) + (freqscale_fac * freqcalib_fac) #assuming CMB is the same and calibrations factors are same for all channel

    if nspecs>1:
        acap_full = np.zeros( (nspecs, len(acap) * nspecs) )
        acap_full[0,:len(acap)] = acap
        if final_comp.lower() == 'cmb':
            acap_full[1,len(acap):] = acap
        else: #polarisation weights are zero for other foregrounds
            acap_full[1,len(acap):] = 0.

        acap_full = np.mat(acap_full).T #should be nspecs*nc x nspecs
        acap = acap_full
    else:
        acap = np.mat(acap).T #should be nspecs*nc x nspecs
    
    return acap

def get_cib_freq_dep(nu, 
    Tcib = 20., 
    Tcmb = 2.7255, 
    beta = 1.505,
    ):
    nu *= 1e9
    bnu1 = fg.get_BnuT(nu, temp = Tcib)
    dbdt = fg.get_dB_dT(nu)
    value = (nu**beta) * bnu1 / dbdt

    return value

def get_radio_freq_dep(nu, 
    nu0 = 150., 
    spec_index_rg = -0.76, 
    null_highfreq_radio = True,
    highfreq_radio_threhsold = 230,
    ):
    nu *= 1e9
    nu0 *= 1e9
    
    nr = fg.get_dB_dT(nu0)
    dr = fg.get_dB_dT(nu)
    epsilon_nu1_nu0 = nr/dr
    scaling = (nu/nu0)**spec_index_rg
    value = epsilon_nu1_nu0 * scaling

    if null_highfreq_radio and (nu>highfreq_radio_threhsold):
        value = 0.

    return value

def get_teb_spec_combination(cl_dic):
    """
    uses cl_dict to determine if we are using ILC jointly for T/E/B.

    Parameters
    ----------
    cl_dic: dict
        dictionary containing (signal+noise) auto- and cross- spectra of different freq. channels at all \ell.
        Keys must be TT, EE, TE, etc.

    Returns
    -------
    nspecs : int
        tells if we are performing ILC for T alone or T/E/B together.
        default is 1. For only one map component.

    specs : list
        creates ['TT', 'EE', 'TE', ... etc.] based on cl_dict that is supplied.
        For example:
        ['TT'] = ILC for T-only
        ['EE'] = ILC for E-only
        ['TT', 'EE'] = ILC for T and E separately.
        ['TT', 'EE', 'TE'] = ILC for T and E jointly.
    """

    specs = sorted( list( cl_dic.keys() ) )

    if specs == ['TT'] or specs == ['EE']: #only TT is supplied
        nspecs = 1
    elif specs == sorted(['TT', 'EE']) or specs == sorted(['TT', 'EE', 'TE']): #TT/EE/TE are supplied
        nspecs = 2
    elif specs == sorted(['TT', 'EE', 'BB']) or specs == sorted(['TT', 'EE', 'BB', 'TE', 'TB', 'EB']): #TT/EE/BB are supplied
        nspecs = 3
    else:
        nspecs = 1
        specs = None

    return nspecs, specs

def corr_from_cov(covmat):
    """
    Get correlation matrix from covariance matrix.

    Parameters
    ----------
    covmat: array
        Covariance matrix.

    Returns
    -------
    corrmat: array
        Correlation matrix.
    """

    diags = np.diag(covmat)
    corrmat = np.zeros_like(covmat)
    for i in range(covmat.shape[0]):
        for j in range(covmat.shape[0]):
            corrmat[i, j] = covmat[i, j] / np.sqrt( diags[i] *  diags[j] )
    #corrmat = covmat / np.outer(np.sqrt(diags), np.sqrt(diags))
    return corrmat

def create_clmat(freqarr, elcnt, cl_dic):
    """
    Get the inverse covariance matrix at a specific multipole moment \ell.
    
    Parameters
    ----------
    freqarr: list
        Frequency array
    elcnt: int
        \ell index.
    cl_dic: dict
        dictionary containing (signal+noise) auto- and cross- spectra of different freq. channels at all \ell.
        Keys must be TT, EE, TE, etc.

    Returns
    -------
    clmat: array
        Covariance matrix at the \ell index.
    """
    nc = len(freqarr)
    nspecs, pspec_arr = get_teb_spec_combination(cl_dic)
    clmat = np.zeros( (nspecs * nc, nspecs * nc) )

    for pspecind, pspec in enumerate( pspec_arr ):
        curr_cl_dic = cl_dic[pspec]
        if nspecs == 1: #clmat for TT or EE or BB
            for ncnt1, freq1 in enumerate(freqarr):
                for ncnt2, freq2 in enumerate(freqarr):
                    j, i = ncnt2, ncnt1
                    clmat[j, i] = curr_cl_dic[(freq1, freq2)][elcnt]
        else: #joint or separate TT/EE constraints #fix me: include BB for joint constraints.
            if pspec == 'TT':
                for ncnt1, freq1 in enumerate(freqarr):
                    for ncnt2, freq2 in enumerate(freqarr):
                        j, i = ncnt2, ncnt1
                        clmat[j, i] = curr_cl_dic[(freq1, freq2)][elcnt]
            elif pspec == 'EE':
                for ncnt1, freq1 in enumerate(freqarr):
                    for ncnt2, freq2 in enumerate(freqarr):
                        j, i = ncnt2 + nc, ncnt1 + nc
                        clmat[j, i] = curr_cl_dic[(freq1, freq2)][elcnt]
            elif pspec == 'TE':
                for ncnt1, freq1 in enumerate(freqarr):
                    for ncnt2, freq2 in enumerate(freqarr):
                        j, i = ncnt2 + nc, ncnt1
                        clmat[j, i] = curr_cl_dic[(freq1, freq2)][elcnt]
                        clmat[i, j] = curr_cl_dic[(freq1, freq2)][elcnt]
        
    return clmat

def get_clinv(freqarr, 
    elcnt, 
    cl_dic, 
    return_clmat = False,
    ):
    
    """
    Get the inverse covariance matrix at a specific multipole moment \ell.
    
    Parameters
    ----------
    freqarr: list
        Frequency array
    elcnt: int
        \ell index.
    cl_dic: dict
        dictionary containing (signal+noise) auto- and cross- spectra of different freq. channels at all \ell.
        Keys must be TT, EE, TE, etc.
    return_clmat: bool
        If True, return covariance at a specific ell. (For debugging pruposes)
        Default is False

    Returns
    -------
    clinv: array
        Inverse covariance matrix at the given \ell index.
    clmat: array
        Covariance matrix at the \ell index.
        Only returned if return_clmat is set to True.
    """

    clmat = np.mat( create_clmat(freqarr, elcnt, cl_dic) )
    clinv = np.linalg.pinv(clmat)

    if return_clmat:
        return clinv, clmat
    else:
        return clinv

def residual_power(param_dict, 
    freqarr, 
    el, 
    cl_dic, 
    final_comp = 'cmb', 
    null_comp = None, 
    spec_index_rg = -0.76, 
    freqcalib_fac = None, 
    lmin = 10, 
    return_weights = True, 
    ):

    """
    Get the residual ILC power.
    
    Parameters
    ----------
    freqarr: list
        Frequency array
    el: int
        \Multipole moment ell.
    cl_dic: dict
        dictionary containing (signal+noise) auto- and cross- spectra of different freq. channels at all \ell.
        Keys must be TT, EE, TE, etc.
    final_comp: str
        Name of the signal that is being minimised.
        Default is 'cmb'.
        Options are ['cmb', 'tsz', 'y', 'ksz', 'radio', 'dust', 'noise', 'tsz_cib']
    null_comp: str
        Name of the signal that is being nulled.
        Default is None.
        Options are ['cmb', 'tsz', 'y', 'ksz', 'radio', 'dust', 'noise', 'tsz_cib']
    spec_index_rg: float
        Default radio spectral index set to -0.76.
    freqcalib_fac: array
        array containing calibration factors / mis-matches between different bands.
        default is None which is equivalent to [1, 1, 1, ... 1] in all bands.
    lmin: int
        Minimum multipole moment to use.
        Default is 10.
    return_weights: bool
        If True return ILC weights along with the residuals


    Returns
    -------
    cl_residual: array
        ILC residuals defined over the desired multiples.
    weightsarr: array
        ILC weights defined over the desired multiples for all the frequency bands.
        Only returned if return_weights is set to True.
    """
    
    assert final_comp in ['cmb', 'tsz', 'y', 'ksz', 'radio', 'dust', 'noise', 'tsz_cib']
    assert null_comp in [None, 'cmb', 'tsz', 'y', 'ksz', 'radio', 'dust', 'noise', 'tsz_cib']
    if freqcalib_fac is not None:
        assert len(freqcalib_fac) == len(freqarr)

    nspecs, pspec_arr = get_teb_spec_combination(cl_dic) #20200527 - teb
    acap = get_acap(freqarr, final_comp = final_comp, freqcalib_fac = freqcalib_fac, nspecs = nspecs)

    if null_comp is not None:
        total_comp_to_null = 0
        if np.ndim(null_comp) == 0:
            bcap = get_acap(freqarr, final_comp = null_comp, freqcalib_fac = freqcalib_fac, nspecs = nspecs)
            total_comp_to_null += 1
        else:
            bcap = None
            for curr_null_comp in null_comp:
                curr_bcap = get_acap(freqarr, final_comp = curr_null_comp, freqcalib_fac = freqcalib_fac, nspecs = nspecs)
                if bcap is None:
                    bcap = np.copy( curr_bcap )
                else:
                    bcap = np.column_stack( (bcap, curr_bcap) )
                total_comp_to_null += 1
            bcap = np.mat(bcap)

    nc = len(freqarr)
    weightsarr = np.zeros( (nspecs * nc, nspecs, len( el ) ) )
    cl_residual = np.zeros( (3, len(el)) )

    cl_residual_tmp = []

    for elcnt, currel in enumerate(el):
        if currel <= lmin: continue ## or el>=lmax: continue
        #clinv = get_clinv( freqarr, elcnt, cl_dic )
        clinv, clmat = get_clinv( freqarr, elcnt, cl_dic, return_clmat = True )

        if null_comp is None:
            nr = np.dot(clinv, acap)
            dr = np.dot( acap.T, np.dot(clinv, acap) )
            drinv = np.linalg.pinv(dr)
            weight = np.dot(nr, drinv)
        else:

            G = np.column_stack( (acap, bcap) )
            G = np.mat(G)

            total_comps = G.shape[1]
            ncap = np.zeros( total_comps )#total_comp_to_null + 1 )
            ncap[0] = 1.
            ncap = np.mat( ncap ).T

            nr = np.dot(clinv, G)
            dr = np.dot( G.T, np.dot(clinv, G) )
            drinv = np.dot( np.linalg.pinv(dr), ncap)
            weight = np.dot(nr, drinv)

        #ILC residuals
        if nspecs>1:
            cl_residual_tt, cl_residual_ee, cl_residual_te = drinv[0,0], drinv[1,1], drinv[0,1]
            cl_residual[:, elcnt] = cl_residual_tt, cl_residual_ee, cl_residual_te
        else:
            cl_residual[0, elcnt] = drinv[0]


        weightsarr[:, :, elcnt] = weight

        cl_residual_tmp.append( drinv )

    weightsarr = np.asarray( weightsarr )
    cl_residual = np.asarray( cl_residual )

    cl_residual[np.isinf(cl_residual)] = 0.
    cl_residual[np.isnan(cl_residual)] = 0.
    
    if return_weights:
        return cl_residual, weightsarr
    else:
        return cl_residual

################################################################################################################

def coth(x):
    """
    Coth function for tSZ frequency response

    Parameters
    ----------
    x: float
        x = h*nu/(k_B * T_CMB)

    Returns
    -------
    coth_x: float
        Hyperbolic cotanget of x.
    """
    
    coth_x = (np.exp(x) + np.exp(-x)) / (np.exp(x) - np.exp(-x))
    
    return coth_x

################################################################################################################
def compton_y_to_delta_Tcmb(nu):
    """
    Get the tSZ frequency response
    
    Parameters
    ----------
    nu: float
        Frequency in Hz.

    Returns
    -------
    g_nu_with_tcmb: float
        tSZ frequency response in Tcmb units.
    """
    nu *= 1e9
    x = h * nu / k_B / Tcmb
    g_nu = x * coth(x/2.) - 4.
    g_nu_with_tcmb = Tcmb * np.mean(g_nu)

    return g_nu_with_tcmb

