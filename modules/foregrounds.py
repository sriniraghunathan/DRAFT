import numpy as np, sys, os, scipy as sc
from scipy import interpolate as intrp
from scipy import ndimage
from scipy.optimize import curve_fit
from pylab import *
import ilc

h, k_B, c=6.62607004e-34, 1.38064852e-23, 2.99792458e8
data_folder = 'data/'

################################################################################################################

def compton_y_to_delta_Tcmb(freq1, freq2 = None, Tcmb = 2.73):

    if freq1<1e4: freq1 = freq1 * 1e9
    if not freq2 is None:
        if freq2<1e4: freq2 = freq2 * 1e9
        freq = np.arange(freq1,freq2,delta_nu)
    else:
        freq = np.asarray([freq1])

    x = (h * freq) / (k_B * Tcmb)
    g_nu = x * coth(x/2.) - 4.

    return Tcmb * np.mean(g_nu)

def get_cl_dust(freq1, freq2, fg_model = 'george15', freq0 = 150, spec_index_dg_po = 1.505 - 0.077, spec_index_dg_clus = 2.51-0.2, Tcib = 20., reduce_cib_power = None):    
    
    assert fg_model in ['george15', 'reichardt21']
    el, cl_dg_po_freq0 = get_foreground_power_spt('DG-Po', freq1 = freq0, freq2 = freq0)
    el, cl_dg_clus_freq0 = get_foreground_power_spt('DG-Cl', freq1 = freq0, freq2 = freq0)
    el_norm = 3000

    #conert to Dls
    dl_fac = el * (el+1)/2/np.pi
    dl_dg_po = dl_fac * cl_dg_po_freq0
    dl_dg_clus = dl_fac * cl_dg_clus_freq0

    if reduce_cib_power: #reduce 150 GHz CIB power: useful for CMB-HD
        dl_dg_po = dl_dg_po/reduce_cib_power
        dl_dg_clus = dl_dg_clus/reduce_cib_power

    nr = ( get_dB_dT(freq0) )**2.
    dr = get_dB_dT(freq1) * get_dB_dT(freq2)

    epsilon_nu1_nu2 = nr/dr

    bnu1 = get_BnuT(freq1, temp = Tcib)
    bnu2 = get_BnuT(freq2, temp = Tcib)
    bnu0 = get_BnuT(freq0, temp = Tcib)

    etanu1_dg_po = ((1.*freq1*1e9)**spec_index_dg_po) * bnu1
    etanu2_dg_po = ((1.*freq2*1e9)**spec_index_dg_po) * bnu2
    etanu0_dg_po = ((1.*freq0*1e9)**spec_index_dg_po) * bnu0

    etanu1_dg_clus = ((1.*freq1*1e9)**spec_index_dg_clus) * bnu1
    etanu2_dg_clus = ((1.*freq2*1e9)**spec_index_dg_clus) * bnu2
    etanu0_dg_clus = ((1.*freq0*1e9)**spec_index_dg_clus) * bnu0

    dl_dg_po = dl_dg_po[el == el_norm][0] * epsilon_nu1_nu2 * (1.*etanu1_dg_po * etanu2_dg_po/etanu0_dg_po/etanu0_dg_po) * (el*1./el_norm)**2
    dl_dg_clus = dl_dg_clus[el == el_norm][0] * epsilon_nu1_nu2 * (1.*etanu1_dg_clus * etanu2_dg_clus/etanu0_dg_clus/etanu0_dg_clus) * (el*1./el_norm)**0.8

    cl_dg_po = dl_dg_po / dl_fac
    cl_dg_clus = dl_dg_clus / dl_fac

    cl_dg_po[np.isnan(cl_dg_po)] = 0.
    cl_dg_clus[np.isinf(cl_dg_clus)] = 0.

    return el, cl_dg_po, cl_dg_clus

def scale_cl_dust(el, cl_dust_freq0, freq0, freq1, freq2, beta, Tcib, el_slope, el_norm = 3000):
        
    #conert to Dls
    dl_fac = el * (el+1)/2/np.pi
    dl_dust_freq0 = dl_fac * cl_dust_freq0

    nr = ( get_dB_dT(freq0) )**2.
    dr = get_dB_dT(freq1) * get_dB_dT(freq2)

    epsilon_nu1_nu2 = nr/dr

    bnu1 = get_BnuT(freq1, temp = Tcib)
    bnu2 = get_BnuT(freq2, temp = Tcib)
    bnu0 = get_BnuT(freq0, temp = Tcib)

    etanu1 = ((1.*freq1*1e9)**beta) * bnu1
    etanu2 = ((1.*freq2*1e9)**beta) * bnu2
    etanu0 = ((1.*freq0*1e9)**beta) * bnu0

    dl_dust = dl_dust_freq0[el == el_norm][0] * epsilon_nu1_nu2 * (1.*etanu1 * etanu2 / etanu0 /etanu0) * (el*1./el_norm)**el_slope

    cl_dust = dl_dust / dl_fac
    cl_dust[np.isnan(cl_dust)] = 0.

    return el, cl_dust

def get_cl_tsz(freq1, freq2, freq0 = 150, fg_model = 'george15', reduce_tsz_power = None):

    assert fg_model in ['george15', 'reichardt21']
    el, cl_tsz_freq0 = get_foreground_power_spt('tSZ', freq1 = freq0, freq2 = freq0)

    tsz_fac_freq0 = compton_y_to_delta_Tcmb(freq0*1e9)
    tsz_fac_freq1 = compton_y_to_delta_Tcmb(freq1*1e9)
    tsz_fac_freq2 = compton_y_to_delta_Tcmb(freq2*1e9)

    scalefac = tsz_fac_freq1 * tsz_fac_freq2/ (tsz_fac_freq0**2.)

    cl_tsz = cl_tsz_freq0 * scalefac
    cl_tsz[np.isnan(cl_tsz)] = 0.
    cl_tsz[np.isinf(cl_tsz)] = 0.

    if reduce_tsz_power is not None:
        cl_tsz /= reduce_tsz_power

    return el, cl_tsz

def get_cl_tsz_cib(freq1, freq2, freq0 = 150, fg_model = 'george15', spec_index_dg_po = 1.505 - 0.077, spec_index_dg_clus = 2.51-0.2, Tcib = 20., cl_cib_dic = None, reduce_tsz_power = None, cib_flux_threshold = 1.5):

    assert fg_model in ['george15', 'reichardt21']
    if fg_model == 'george15':
        corr_coeff = 0.1
    elif fg_model == 'reichardt20':
        corr_coeff = 0.078

    el, cl_tsz_freq1_freq1 = get_cl_tsz(freq1, freq1, freq0 = freq0, fg_model = fg_model, reduce_tsz_power = reduce_tsz_power)
    if cl_cib_dic is not None:
        el, cl_dg_freq1_freq1 = cl_cib_dic[(freq1, freq1)]
    else:
        #get tSZ and CIB spectra for freq1
        el, cl_dg_po_freq1_freq1, cl_dg_clus_freq1_freq1 = get_cl_dust(freq1, freq1, freq0 = freq0, fg_model = fg_model, spec_index_dg_po = spec_index_dg_po, spec_index_dg_clus = spec_index_dg_clus, Tcib = Tcib)
        cl_dg_freq1_freq1 = cl_dg_po_freq1_freq1 + cl_dg_clus_freq1_freq1

    #get tSZ and CIB spectra for freq2
    el, cl_tsz_freq2_freq2 = get_cl_tsz(freq2, freq2, freq0 = freq0, fg_model = fg_model)
    if cl_cib_dic is not None:
        el, cl_dg_freq2_freq2 = cl_cib_dic[(freq2, freq2)]
    else:
        el, cl_dg_po_freq2_freq2, cl_dg_clus_freq2_freq2 = get_cl_dust(freq2, freq2, freq0 = freq0, fg_model = fg_model, spec_index_dg_po = spec_index_dg_po, spec_index_dg_clus = spec_index_dg_clus, Tcib = Tcib)
        cl_dg_freq2_freq2 = cl_dg_po_freq2_freq2 + cl_dg_clus_freq2_freq2
    if len(el) != len(cl_tsz_freq2_freq2):
        cl_tsz_freq1_freq1 = np.interp(el, np.arange(len(cl_tsz_freq1_freq1)),cl_tsz_freq1_freq1)
        cl_tsz_freq2_freq2 = np.interp(el, np.arange(len(cl_tsz_freq2_freq2)),cl_tsz_freq2_freq2)

    #20220325
    if freq1 >= 217 and freq2 >=217:
        corr_coeff = corr_coeff * -1.

    cl_tsz_cib = -corr_coeff * ( np.sqrt(cl_tsz_freq1_freq1 * cl_dg_freq2_freq2) + np.sqrt(cl_tsz_freq2_freq2 * cl_dg_freq1_freq1) )

    return el, cl_tsz_cib

def get_cl_radio(freq1, freq2, freq0 = 150, fg_model = 'george15', spec_index_rg = -0.9, null_highfreq_radio = 1, reduce_radio_power_150 = None):

    if fg_model == 'george15':
        el, cl_rg_freq0 = get_foreground_power_spt('RG', freq1 = freq0, freq2 = freq0)
        if reduce_radio_power_150 is not None:
            cl_rg_freq0 /= reduce_radio_power_150
        el_norm = 3000

    #conert to Dls
    dl_fac = el * (el+1)/2/np.pi
    dl_rg = dl_fac * cl_rg_freq0

    nr = ( get_dB_dT(freq0) )**2.
    dr = get_dB_dT(freq1) * get_dB_dT(freq2)

    epsilon_nu1_nu2 = nr/dr

    dl_rg = dl_rg[el == el_norm][0] * epsilon_nu1_nu2 * (1.*freq1 * freq2/freq0/freq0)**spec_index_rg * (el*1./el_norm)**2

    cl_rg = dl_rg / dl_fac

    cl_rg[np.isnan(cl_rg)] = 0.

    if null_highfreq_radio and (freq1>230 or freq2>230):
        cl_rg *= 0.

    return el, cl_rg

def smooth_cib_spectra(el, cl, min_el = 200, el_knee_guess = 2000):
    el_ori = np.copy(el)
    non_zero_ind = np.where(el>min_el)[0]
    el = el[non_zero_ind]
    cl = cl[non_zero_ind]

    el = el.astype(np.float64)
    def fitting_func_dust(p, x, DATA = None, return_fit = 0, el_slope = 1.2):
        fitfunc = lambda p, x: p[0]*( 1 + (p[1]/x)**p[2])
        #fitfunc = lambda p, x: p[0]*( 1 + (p[1]/x)**el_slope)
        if not return_fit:
            val = fitfunc(p, x) - DATA
            val[np.isinf(val)] = 0.
            val[np.isnan(val)] = 0.
            return val
        else:
            return fitfunc(p, x)

    def fitfunc(x, p1, p2, p3):        
        return p1*( 1. + (p2/x)**p3)

    import scipy.optimize as optimize    
    poi_level = np.median(cl[el>el_knee_guess])
    el_slope = 1.2 #0.8 for Dl and 1.2 for Cl
    p0 = np.asarray( [poi_level, el_knee_guess, el_slope] )
    p1, cov, infodict, success, msg = optimize.leastsq(fitting_func_dust, p0, args=(el, cl), full_output = 1)
    if p1[1]<0:
        p1 = p0
    cl_fit = fitting_func_dust(p1, el, return_fit = 1)

    cl_fit = np.interp(el_ori, el, cl_fit, left = 0., right = 0.)

    return cl_fit

def power_law(ell, A, alpha, ell_norm = 80.):
    #print(A, alpha)#, end = ' ')
    fit = A * ((ell / ell_norm) ** alpha)
    badinds = np.where((fit == np.inf) | (fit == np.nan))[0]
    fit[badinds]=0.
    return fit

def perform_power_law_fit(el, cl, ell_norm = 80):
    dl_fac = (el * (el + 1))/2/np.pi
    dl = cl * dl_fac
    badinds = np.where((dl == np.inf) | (dl == np.nan))[0]
    dl[badinds]=0.
    amp_ini = dl[el == ell_norm][0]
    alpha_ini = -.3
    delta_alpha = 0.1
    amp_low_fac, amp_high_fac = 0.95, 1.05 #0.1, 3.
    pars, cov = curve_fit(f=power_law, xdata=el, ydata=dl, p0=[amp_ini, alpha_ini], bounds = ((amp_ini*amp_low_fac, alpha_ini-delta_alpha), (amp_ini*amp_high_fac, alpha_ini+delta_alpha)))
    
    dl_fit = power_law(el, pars[0], pars[1])
    cl_fit = dl_fit / dl_fac
    cl_fit[np.isinf(cl_fit)] = 0.
    cl_fit[np.isnan(cl_fit)] = 0.

    return cl_fit

def get_cl_galactic(param_dict, component, freq1, freq2, which_spec, which_gal_mask = 0, bl_dic = None, el = None, use_power_law_fit = False, use_sed_scaling = True, freq0_for_sed_scaling = 278., ell_norm = 80., Tdust = 20., beta_dust = 1.54):

    gal_freq_dic = {20:20, 27:27, 39: 39, 93: 93, 90: 93, 143: 143, 145: 145, 150: 150, 225: 225, 220: 225, 278:278, 350: 350}

    #https://healpy.readthedocs.io/en/1.5.0/generated/healpy.sphtfunc.anafast.html#healpy.sphtfunc.anafast
    spec_inds_dic = { 'TT':0, 'EE':1, 'BB':2, 'TE':3, 'EB':4, 'TB':5} #py2

    assert component in ['dust', 'sync', 'freefree']

    try:
        which_gal_mask = param_dict['which_gal_mask']
    except:
        pass

    if component == 'dust':
        cl_gal_dic_fname = param_dict['cl_gal_dic_dust_fname']
    elif component == 'sync':
        cl_gal_dic_fname = param_dict['cl_gal_dic_sync_fname']
    elif component == 'freefree':
        cl_gal_dic_fname = param_dict['cl_gal_dic_freefree_fname']

    try:
        cl_gal_folder = param_dict['cl_gal_folder']
        cl_gal_dic_fname = '%s/%s' %(cl_gal_folder, cl_gal_dic_fname)
    except:
        pass

    if (0):##component == 'sync':
        #fix me: Forcing sync. to CUmilta's simulations
        print('\n\t\tForcing sync. to CUmilta\'s simulations\n\n')
        try:
            cl_gal_dic_fname = param_dict['cl_gal_dic_sync_fname_forced']
        except:
            pass

    #pick the requested spectra: TT, EE, BB, TE, EB, TB.
    spec_ind = spec_inds_dic[which_spec]

    if cl_gal_dic_fname.find('spt_proposal_2023_13k_sqdeg_field_mask')>-1: #20230330
        gal_freq_dic = {90: 93, 95: 93, 150: 145, 220: 225}

    freq1_to_use = gal_freq_dic[freq1]
    freq2_to_use = gal_freq_dic[freq2]

    cl_gal_dic = np.load(cl_gal_dic_fname, allow_pickle = 1, encoding = 'latin1').item()['cl_dic'][which_gal_mask]
    freq0_for_sed_scaling = 278.
    if (freq0_for_sed_scaling, freq0_for_sed_scaling) not in cl_gal_dic:
        use_sed_scaling = False
    if ( freq1_to_use >= max(list(gal_freq_dic)) or freq2_to_use >= max(list(gal_freq_dic)) ) and use_sed_scaling:
        cl_dust_freq0 = cl_gal_dic[ (freq0_for_sed_scaling, freq0_for_sed_scaling) ]
        if component == 'dust':
            cl_gal = scale_cl_dust_galactic(cl_dust_freq0, freq1, freq2 = freq2, freq0 = freq0_for_sed_scaling, Tdust = Tdust, beta_dust = beta_dust)
        else:
            cl_gal = np.zeros(cl_dust_freq0.shape)
    else:
        try:
            cl_gal = cl_gal_dic[ (freq1_to_use, freq2_to_use) ]
        except:
            cl_gal = cl_gal_dic[ (freq2_to_use, freq1_to_use) ]

        #fix me
        if np.ndim(cl_gal) == 1: #TT-only. Pol will fail.
            cl_gal = np.asarray( [cl_gal] )

    if which_spec == 'TE' and cl_gal_dic_fname.find('CUmilta')==-1:
        #force TE to be np.sqrt(TT) * np.sqrt(EE)
        cl_gal_tt, cl_gal_ee = cl_gal_dic[ (freq1, freq2) ][0], cl_gal_dic[ (freq1, freq2) ][1]

        if (1):##component == 'dust':
            rte = 0.35 #page 5 of https://arxiv.org/pdf/1801.04945.pdf: Discussion below Fig.5; also page 38 of https://readthedocs.org/projects/so-pysm-models/downloads/pdf/0.2.dev/
        else: ##elif component == 'sync':
            rte = 0.
        cl_gal = rte * np.sqrt( cl_gal_tt * cl_gal_ee )

    else:
        try:
            cl_gal = cl_gal[spec_ind]
        except:
            print('(%s,%s) not found for mask = %s in %s. Setting them to zeros.' %(freq1, freq2, which_spec, cl_gal_dic_fname))
            cl_gal = np.zeros( len(cl_gal[0]) )


    el_gal = np.arange( len(cl_gal) )

    #20210506: SED scaling/power-law fitting. By default we do not fit power law.
    beam_deconvolved = False
    if use_sed_scaling and component == 'dust':
        cl_dust_freq0 = cl_gal_dic[ (freq0_for_sed_scaling, freq0_for_sed_scaling) ][spec_ind]
        el_gal = np.arange( len(cl_dust_freq0) )
        if bl_dic is not None:        
            bl = bl_dic[freq0_for_sed_scaling]
            if len(bl) != len(cl_dust_freq0): #adjust array lengths first
                el_tmp = np.arange( len(bl) )
                cl_dust_freq0 = np.interp(el_tmp, el_gal, cl_dust_freq0, left = 0., right = 0.)
                el_gal = np.copy( el_tmp )
            cl_dust_freq0 = cl_dust_freq0 / bl**2.
            beam_deconvolved = True

        if use_power_law_fit:
            cl_dust_freq0 = perform_power_law_fit(el_gal, cl_dust_freq0, ell_norm = ell_norm)
        cl_gal = scale_cl_dust_galactic(cl_dust_freq0, freq1_to_use, freq2 = freq2_to_use, freq0 = freq0_for_sed_scaling, Tdust = Tdust, beta_dust = beta_dust)
    elif use_power_law_fit:
        if bl_dic is not None:        
            bl = bl_dic[freq0_for_sed_scaling]
            if len(bl) != len(cl_gal): #adjust array lengths first
                el_tmp = np.arange( len(bl) )
                cl_gal = np.interp(el_tmp, el_gal, cl_gal, left = 0., right = 0.)
                el_gal = np.copy( el_tmp )
            cl_gal = cl_gal / bl**2.
            beam_deconvolved = True

        cl_gal = perform_power_law_fit(el_gal, cl_gal, ell_norm = ell_norm)

    if bl_dic is not None and beam_deconvolved is False:
        bl1 = bl_dic[freq1]
        bl2 = bl_dic[freq2]

        if len(bl1) != len(cl_gal): #adjust array lengths first
            el_tmp = np.arange( len(bl1) )
            cl_gal = np.interp(el_tmp, el_gal, cl_gal, left = 0., right = 0.)
            el_gal = np.copy( el_tmp )

        cl_gal = cl_gal / (bl1 * bl2)

        if (0): #20210426 - nulling highly smoothed modes
            op_beam = bl_dic[145]
            beam_ratio = op_beam**2. / (bl1 * bl2)
            highly_deconv_inds = np.where(beam_ratio>=500)
            cl_gal[highly_deconv_inds] = 0.
            #plot(beam_ratio); ylim(1., 1000.); show(); sys.exit()

    if el is not None:
        cl_gal = np.interp(el, el_gal, cl_gal, left = 0., right = 0.)
        el_gal = np.copy( el )

    return el_gal, cl_gal

def scale_cl_dust_galactic(cl, freq1, freq2 = None, freq0 = 278., Tdust = 19.6, beta_dust = 1.6):

    if freq2 is None:
        freq2 = freq1

    nr = ( get_dB_dT(freq0) )**2.
    dr = get_dB_dT(freq1) * get_dB_dT(freq2)

    epsilon_nu1_nu2 = nr/dr

    bnu1 = get_BnuT(freq1, temp = Tdust)
    bnu2 = get_BnuT(freq2, temp = Tdust)
    bnu0 = get_BnuT(freq0, temp = Tdust)

    etanu1_dust = ((1.*freq1*1e9)**beta_dust) * bnu1
    etanu2_dust = ((1.*freq2*1e9)**beta_dust) * bnu2
    etanu0_dust = ((1.*freq0*1e9)**beta_dust) * bnu0

    cl_dust = cl * epsilon_nu1_nu2 * (1.*etanu1_dust * etanu2_dust/etanu0_dust/etanu0_dust)## * (el*1./el_norm)**el_slope

    return cl_dust

def get_cl_dust_galactic(el, freq1, freq2 = None, freq0 = 353., el_norm = 80., el_slope = -0.58, Tdust = 19.6, Adust_freq0 = 4.3, spec_index_dust = 1.6, return_dl = 0):

    if freq2 is None:
        freq2 = freq1

    nr = ( get_dB_dT(freq0) )**2.
    dr = get_dB_dT(freq1) * get_dB_dT(freq2)

    epsilon_nu1_nu2 = nr/dr

    bnu1 = get_BnuT(freq1, temp = Tdust)
    bnu2 = get_BnuT(freq2, temp = Tdust)
    bnu0 = get_BnuT(freq0, temp = Tdust)

    etanu1_dust = ((1.*freq1*1e9)**spec_index_dust) * bnu1
    etanu2_dust = ((1.*freq2*1e9)**spec_index_dust) * bnu2
    etanu0_dust = ((1.*freq0*1e9)**spec_index_dust) * bnu0

    dl_dust = Adust_freq0 * epsilon_nu1_nu2 * (1.*etanu1_dust * etanu2_dust/etanu0_dust/etanu0_dust) * (el*1./el_norm)**el_slope

    if return_dl:
        return dl_dust
    else:
        dl_fac = el * (el+1)/2/np.pi
        cl_dust = dl_dust / dl_fac
    return cl_dust

def get_foreground_power_spt(component, freq1=150, freq2=None, units='uk', lmax = None):
    """
    Foreground powers from George et al. 2015 results.

    Uses .sav file generated by Christain Reichardt.

    Parameters
    ----------
    component : str
        The foreground component to use. Must be one of
        'all', 'tSZ', 'kSZ', 'DG-Cl', 'DG-Po', 'RG', 'tSZ-CIB', 'Total', 'CMB'
    freq1 : int
        Frequency band. If `freq2` is specified, the cross-spectrum between
        the two frequencies will be returned. Otherwise autospectrum of freq1.
    freq2 : int, optional
        Frequency band for cross-spectrum with `freq1`
    units : str
        'k' or 'uk'. Note: default savfile is Dls in uK

    Returns
    -------
    fgnd_cls : array
        Power spectrum of `component` at specified frequency band.
    """
    components = [
        'all',
        'tSZ',
        'kSZ',
        'DG-Cl',
        'DG-Po',
        'RG',
        'tSZ-CIB',
        'Total',
        'CMB',
    ]
    if component not in components:
        raise ValueError(
            '{} not in list of possible foregrounds, must be one of {}'.format(
                component, components
            )
        )

    from scipy.io import readsav
    filename = '%s/george_plot_bestfit_line.sav' %(data_folder)
    data = readsav(filename)

    if freq2 is None:
        freq2 = freq1
    if freq1 == 90:
        freq1 = 95
    if freq2 == 90:
        freq2 = 95

    freqs = np.asarray(
        [(95, 95), (95, 150), (95, 220), (150, 150), (150, 220), (220, 220)]
    )
    dl_all = data['ml_dls'][(freqs[:, 0] == freq1) & (freqs[:, 1] == freq2)][0]
    labels = data['ml_dl_labels'].astype('str')
    el = np.asarray(data['ml_l'], dtype=int)

    if component == 'all':
        spec = el * 0.0
        for fg in components:
            if fg in ['all', 'tSZ-CIB', 'Total', 'CMB']:
                continue
            spec += dl_all[labels == fg][0]
    else:
        spec = dl_all[labels == component][0]

    # Changing Dls to Cls
    spec /= el * (el + 1.0) / 2.0 / np.pi
    if units.lower() == 'k':
        spec /= 1e12

    # Pad to l=0
    spec = np.concatenate((np.zeros(min(el)), spec))
    el = np.concatenate((np.arange(min(el)), el))

    if lmax is not None:
        el = el[:lmax]
        spec = spec[:lmax]

    return el, spec

def get_dB_dT(nu, nu0 = None, temp = 2.725):
    if nu<1e4: nu *= 1e9

    x=h*nu/(k_B*temp)
    dBdT = x**4. * np.exp(x) / (np.exp(x)-1)**2.

    if nu0 is not None:
        nu0 *= 1e9
        x0=h*nu0/(k_B*temp)
        dBdT0 = x0**4 * np.exp(x0) / (np.exp(x0)-1)**2.
        return  dBdT / dbdT0
    else:
        return dBdT

def get_BnuT(nu, temp = 2.725):
    if nu<1e4: nu *= 1e9
    x=h*nu/(k_B*temp)

    t1 = 2 * h * nu**3./ c**2.
    t2 = 1./ (np.exp(x)-1.)

    return t1 * t2

def coth(x):
    return (np.exp(x) + np.exp(-x)) / (np.exp(x) - np.exp(-x))