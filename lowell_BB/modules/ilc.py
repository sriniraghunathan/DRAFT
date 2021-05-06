import numpy as np, sys, os, scipy as sc, foregrounds as fg, misc, re, flatsky#, healpy as H
from pylab import *
################################################################################################################
def get_analytic_covariance(param_dict, freqarr, el = None, nl_dic = None, bl_dic = None, ignore_fg = [], which_spec = 'TT', pol_frac_per_cent_dust = 0.02, pol_frac_per_cent_radio = 0.03, pol_frac_per_cent_tsz = 0., pol_frac_per_cent_ksz = 0., include_gal = 0, null_highfreq_radio = 1, reduce_radio_power_150 = None, reduce_tsz_power = None, reduce_cib_power = None, return_fg_spectra = True, max_nl_value = 5000.):

    #ignore_fg = foreground terms that must be ignored
    possible_ignore_fg = ['cmb', 'tsz', 'ksz', 'radio', 'dust', 'noise', 'tsz_cib']
    if len(ignore_fg)>0:
        if 'cmb' in ignore_fg: ignore_fg.append('ksz')
        if not all( [ currfg in possible_ignore_fg for currfg in ignore_fg] ):
            print( '\n\t Alert: Elements of ignore_fg should be one of the following: %s\n\n' %(np.array2string(np.asarray(possible_ignore_fg))) )
            sys.exit()

    ###el_, cl_cmb = fg.get_foreground_power_spt('CMB', freq1 = param_dict['freq0'], freq2 = param_dict['freq0'])
    camb_file = '%s/%s' %(param_dict['data_folder'], param_dict['Dlfile_len'])
    el_ = np.loadtxt(camb_file, usecols = [0])
    dl_camb = np.loadtxt(camb_file, usecols = [1,2,3,4])
    Tcmb = 2.73
    cl_camb = ( Tcmb**2. * dl_camb * 2 * np.pi ) / ( el_[:,None] * (el_[:,None] + 1) )
    cl_camb *= 1e12
    cl_TT, cl_EE, cl_BB, cl_TE = cl_camb.T
    cl_cmb_dic = {'TT': cl_TT, 'EE': cl_EE, 'BB': cl_BB, 'TE': cl_TE}
    cl_cmb = cl_cmb_dic[which_spec]

    if el is None:
        el = np.copy(el_)
    cl_cmb = np.interp(el, el_, cl_cmb)

    if 'cmb' in ignore_fg: #get CMB spectra
        cl_cmb = cl_cmb * 0.

    if 'ksz' not in ignore_fg:
        el_, cl_ksz = fg.get_foreground_power_spt('kSZ', freq1 = param_dict['freq0'], freq2 = param_dict['freq0'])
        if which_spec == 'EE':
            cl_ksz = cl_ksz * pol_frac_per_cent_ksz**2.
        if which_spec == 'TE':
            cl_ksz = cl_ksz * 0.
        cl_ksz = np.interp(el, el_, cl_ksz) #interpolate to get cl for required els
    else:
        cl_ksz = np.zeros(len(el))

    cl_dic = {}
    if return_fg_spectra:
        fg_cl_dic = {}
        
    for freq1 in freqarr:
        for freq2 in freqarr:

            if (freq2, freq1) in cl_dic:
                cl_dic[(freq1, freq2)] = cl_dic[(freq2, freq1)]
                continue

            if 'tsz' not in ignore_fg: #get tsz
                el_, cl_tsz = fg.get_cl_tsz(freq1, freq2, freq0 = param_dict['freq0'], fg_model = param_dict['fg_model'], reduce_tsz_power = reduce_tsz_power)
                if which_spec == 'EE':
                    cl_tsz = cl_tsz * pol_frac_per_cent_tsz**2.
                elif which_spec == 'TE':
                    cl_tsz = cl_tsz * 0.
                cl_tsz = np.interp(el, el_, cl_tsz) #interpolate to get cl for required els
            else:
                cl_tsz = np.zeros(len(el))

            if 'radio' not in ignore_fg: #get radio
                el_, cl_radio = fg.get_cl_radio(freq1, freq2, freq0 = param_dict['freq0'], fg_model = param_dict['fg_model'], spec_index_rg = param_dict['spec_index_rg'], null_highfreq_radio = null_highfreq_radio, reduce_radio_power_150 = reduce_radio_power_150)
                if which_spec == 'EE':
                    cl_radio = cl_radio * pol_frac_per_cent_radio**2.
                elif which_spec == 'TE':
                    cl_radio = cl_radio * 0.
                cl_radio = np.interp(el, el_, cl_radio) #interpolate to get cl for required els
            else:
                cl_radio = np.zeros(len(el))

            if 'dust' not in ignore_fg: #get CIB
                el_,  cl_dg_po, cl_dg_clus = fg.get_cl_dust(freq1, freq2, freq0 = param_dict['freq0'], fg_model = param_dict['fg_model'], spec_index_dg_po = param_dict['spec_index_dg_po'], spec_index_dg_clus = param_dict['spec_index_dg_clus'], Tcib = param_dict['Tcib'], reduce_cib_power = reduce_cib_power)
                cl_dust = cl_dg_po + cl_dg_clus
                cl_dust[np.isnan(cl_dust)] = 0.
                if which_spec == 'EE':
                    cl_dust = cl_dust * pol_frac_per_cent_dust**2.
                elif which_spec == 'TE':
                    cl_dust = cl_dust * 0.
                cl_dust = np.interp(el, el_, cl_dust) #interpolate to get cl for required els
            else:
                cl_dust = np.zeros(len(el))

            #galaxy
            if include_gal:# and not pol: #get galactic dust and sync
                el_, cl_gal_dust = fg.get_cl_galactic(param_dict, 'dust', freq1, freq2, which_spec, el = el, bl_dic = bl_dic)
                el_, cl_gal_sync = fg.get_cl_galactic(param_dict, 'sync', freq1, freq2, which_spec, el = el, bl_dic = bl_dic)
                cl_gal_dust = np.interp(el, el_, cl_gal_dust) #interpolate to get cl for required els
                cl_gal_sync = np.interp(el, el_, cl_gal_sync) #interpolate to get cl for required els
            else:
                cl_gal_dust = np.zeros(len(el))
                cl_gal_sync = np.zeros(len(el))

            cl = cl_cmb + cl_ksz + cl_tsz + cl_radio + cl_dust + cl_gal_dust + cl_gal_sync

            #noise auto power spectrum
            if nl_dic is not None:
                if (freq1,freq2) in nl_dic:
                    nl = nl_dic[(freq1,freq2)]
                else:
                    nl = nl_dic[freq1]
                    if freq1 != freq2: 
                        nl = np.copy(nl) * 0.

                #interpolate to get cl for required els
                el_ = np.arange(len(nl))
                nl = np.interp(el, el_, nl) 

                if (1): #20201121: remove very large numbers because of beam deconvolution
                    ini_nl = np.median(nl[:100])
                    end_nl = np.median(nl[-100:])
                    if end_nl>ini_nl: #this implies beam deconvolution has made end nl pretty large
                        #having end_nl pretty large introduces covariance inversion issues                        
                        badinds = np.where(nl>=max_nl_value)[0]
                        nl[badinds] = max_nl_value
                        #print(ini_nl, end_nl)
                        #clf();loglog(nl); title('%s: %s,%s' %(which_spec, freq1,freq2)); show()

            else:
                nl = np.zeros(len(el))

            if 'noise' not in ignore_fg:
                if which_spec != 'TE':
                    cl = cl + np.copy(nl)

            if return_fg_spectra:
                if 'cmb' not in ignore_fg:
                    if 'cmb' not in fg_cl_dic: fg_cl_dic['cmb'] = {}
                    fg_cl_dic['cmb'][(freq1, freq2)] = fg_cl_dic['cmb'][(freq2, freq1)] = cl_cmb
                if 'ksz' not in ignore_fg:
                    if 'ksz' not in fg_cl_dic: fg_cl_dic['ksz'] = {}
                    fg_cl_dic['ksz'][(freq1, freq2)] = fg_cl_dic['ksz'][(freq2, freq1)] = cl_ksz
                if 'tsz' not in ignore_fg:
                    if 'tsz' not in fg_cl_dic: fg_cl_dic['tsz'] = {}
                    fg_cl_dic['tsz'][(freq1, freq2)] = fg_cl_dic['tsz'][(freq2, freq1)] = cl_tsz
                if 'radio' not in ignore_fg:
                    if 'radio' not in fg_cl_dic: fg_cl_dic['radio'] = {}
                    fg_cl_dic['radio'][(freq1, freq2)] = fg_cl_dic['radio'][(freq2, freq1)] = cl_radio
                if 'dust' not in ignore_fg:
                    if 'cib' not in fg_cl_dic: fg_cl_dic['cib'] = {}
                    fg_cl_dic['cib'][(freq1, freq2)] = fg_cl_dic['cib'][(freq2, freq1)] = cl_dust
                if 'noise' not in ignore_fg and which_spec != 'TE':
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

            cl_dic[(freq1, freq2)] = cl

    if return_fg_spectra:
        return el, cl_dic, fg_cl_dic
    else:
        return el, cl_dic  

################################################################################################################

def get_acap(freqarr, final_comp = 'cmb', freqcalib_fac = None, teb_len = 1):

    nc = len(freqarr)

    if freqcalib_fac is None: freqcalib_fac = np.ones(nc)

    if final_comp.lower() == 'cmb':
        freqscale_fac = np.ones(nc)

    elif final_comp.lower() == 'tsz' or final_comp.lower() == 'y':

        freqscale_fac = []
        for freq in sorted( freqarr ):
            freqscale_fac.append( fg.compton_y_to_delta_Tcmb(freq * 1e9) )

        freqscale_fac = np.asarray( freqscale_fac )

    elif final_comp.lower() == 'cib' or final_comp.lower() == 'cibpo':
        freqscale_fac = []
        for freq in sorted( freqarr ):
            freqscale_fac.append( get_cib_freq_dep(freq * 1e9) )

        freqscale_fac = np.asarray( freqscale_fac )
        freqscale_fac /= np.max(freqscale_fac)
        
    elif final_comp.lower() == 'cibclus':
        freqscale_fac = []
        for freq in sorted( freqarr ):
            freqscale_fac.append( get_cib_freq_dep(freq * 1e9, beta = 2.505) )

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
            freqscale_fac.append( get_cib_freq_dep(freq * 1e9, Tcib = misc_tcib, beta = misc_beta) )

        freqscale_fac = np.asarray( freqscale_fac )
        freqscale_fac /= np.max(freqscale_fac)

    elif final_comp.lower() == 'radio':
        freqscale_fac = []
        for freq in sorted( freqarr ):
            freqscale_fac.append( get_radio_freq_dep(freq) )

        freqscale_fac = np.asarray( freqscale_fac )

    acap = np.zeros(nc) + (freqscale_fac * freqcalib_fac) #assuming CMB is the same and calibrations factors are same for all channel

    if teb_len>1:
        acap_full = np.zeros( (teb_len, len(acap) * teb_len) )
        acap_full[0,:len(acap)] = acap
        if final_comp.lower() == 'cmb':
            acap_full[1,len(acap):] = acap
        else: #polarisation weights are zero for other foregrounds
            acap_full[1,len(acap):] = 0.

        acap_full = np.mat(acap_full).T #should be teb_len*nc x teb_len
        acap = acap_full
    else:
        acap = np.mat(acap).T #should be teb_len*nc x teb_len
    
    return acap

def get_cib_freq_dep(nu, Tcib = 20., Tcmb = 2.7255, h=6.62607004e-34, k_B=1.38064852e-23, beta = 1.505):
    if nu<1e4: nu *= 1e9

    bnu1 = fg.fn_BnuT(nu, temp = Tcib)
    dbdt = fg.fn_dB_dT(nu)
    value = (nu**beta) * bnu1 / dbdt

    return value

def get_radio_freq_dep(nu, nu0 = 150., spec_index_rg = -0.9, null_highfreq_radio = 1):

    nr = fg.fn_dB_dT(nu0)
    dr = fg.fn_dB_dT(nu)
    epsilon_nu1_nu0 = nr/dr
    scaling = (nu/nu0)**spec_index_rg
    value = epsilon_nu1_nu0 * scaling

    if null_highfreq_radio and (nu>230):
        #print('\n\tthis extrapolation does not work for high freqeuncy radio. Making cl_radio = 0 for these bands.')
        value = 0.

    return value

def get_teb_spec_combination(cl_dic):
    pspec_arr = sorted( list( cl_dic.keys() ) )

    if pspec_arr == ['TT'] or pspec_arr == ['EE'] or pspec_arr == ['BB']: #only TT is supplied
        teb_len = 1
    elif pspec_arr == sorted(['TT', 'EE']) or pspec_arr == sorted(['TT', 'EE', 'TE']): #TT/EE/TE are supplied
        teb_len = 2
    elif pspec_arr == sorted(['TT', 'EE', 'BB']) or pspec_arr == sorted(['TT', 'EE', 'BB', 'TE', 'TB', 'EB']): #TT/EE/BB are supplied
        teb_len = 3
    else:
        teb_len = 1
        pspec_arr = None

    return teb_len, pspec_arr

def create_clmat(freqarr, elcnt, cl_dic):
    """
    freqarr  = array of frequency channel
    elcnt = \el index
    cl_dic = cl_cmb + cl_FG auto and cross spectra of the frequency channels
    """
    nc = len(freqarr)
    teb_len, pspec_arr = get_teb_spec_combination(cl_dic)
    clmat = np.zeros( (teb_len * nc, teb_len * nc) )

    for pspecind, pspec in enumerate( pspec_arr ):
        curr_cl_dic = cl_dic[pspec]

        if teb_len == 1: #clmat for TT or EE or BB
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

def get_clinv(freqarr, elcnt, cl_dic, return_clmat = 0):
    clmat = np.mat( create_clmat(freqarr, elcnt, cl_dic) )
    clinv = np.linalg.pinv(clmat)

    if return_clmat:
        return clinv, clmat
    else:
        return clinv

def residual_power_weights(param_dict, freqarr, el, cl_dic, final_comp = 'cmb', freqcalib_fac = None, lmin = 10, return_weights = 0, null_comp = None):
    teb_len, pspec_arr = get_teb_spec_combination(cl_dic) #20200527 - teb
    acap = get_acap(freqarr, final_comp = final_comp, freqcalib_fac = freqcalib_fac, teb_len = teb_len)

    if null_comp is not None:
        total_comp_to_null = 0
        if np.ndim(null_comp) == 0:
            bcap = get_acap(freqarr, final_comp = null_comp, freqcalib_fac = freqcalib_fac, teb_len = teb_len)
            total_comp_to_null += 1
        else:
            bcap = None
            for curr_null_comp in null_comp:
                curr_bcap = get_acap(freqarr, final_comp = curr_null_comp, freqcalib_fac = freqcalib_fac, teb_len = teb_len)
                if bcap is None:
                    bcap = np.copy( curr_bcap )
                else:
                    bcap = np.column_stack( (bcap, curr_bcap) )
                total_comp_to_null += 1
            bcap = np.mat(bcap)
            bcap_full = np.mat( np.copy( bcap ) )

    nc = len(freqarr)
    weightsarr = np.zeros( (teb_len * nc, teb_len, len( el ) ) )
    cl_residual = np.zeros( (3, len(el)) )

    cl_residual_tmp = []

    for elcnt, currel in enumerate(el):
        if currel <= lmin: continue ## or el>=lmax: continue
        #clinv = get_clinv( freqarr, elcnt, cl_dic )
        clinv, clmat = get_clinv( freqarr, elcnt, cl_dic, return_clmat = 1 )

        if null_comp is None:
            nr = np.dot(clinv, acap)
            dr = np.dot( acap.T, np.dot(clinv, acap) )

            #drinv = sc.linalg.pinv2(dr)
            drinv = np.linalg.pinv(dr)
            weight = np.dot(nr, drinv)
        else:

            G = np.column_stack( (acap, bcap) )
            G = np.mat(G)
            total_comps = G.shape[1]
            ncap = np.zeros( total_comps )#total_comp_to_null + 1 )
            ncap[0] = 1.
            #ncap[1] = 0.5
            ncap = np.mat( ncap ).T

            nr = np.dot(clinv, G)
            dr = np.dot( G.T, np.dot(clinv, G) )

            #drinv = np.dot( sc.linalg.pinv2(dr), ncap)
            drinv = np.dot( np.linalg.pinv(dr), ncap)
            weight = np.dot(nr, drinv)

        #ILC residuals
        if teb_len>1:
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

def apply_ilc_weightsarr(maparr, weightsarr, bl_dic, nside, lmax, full_sky = True, mapparams = None, verbose = False):

    """
    maparr  = [map_90, map_150, ...]
    weightsarr = [weights_90, weights_150, ...]
    bl_dic = {90: bl_90, 150: bl_150, ...}
    mapparams = only for flatsky maps. [nx, ny, dx] --> (ny, nx) = flatskymap.shape, dx = pixel resolution in arcmins.
    """

    #rebeaming
    rebeamarr = misc.rebeam( bl_dic )

    #convert weights to 2D if flat-sky
    if not full_sky:
        assert mapparams is not None
        weightsarr_2D = []
        for currW in weightsarr:
            el = np.arange(len(currW))
            currW_2D = flatsky.cl_to_cl2d(el, currW, mapparams) 
            weightsarr_2D.append(currW_2D)
        weightsarr = np.asarray( weightsarr_2D )

    #modify weights to include rebeam
    weightsarr = weightsarr * rebeamarr

    #get the ilc combined map now
    weighted_almarr = []
    for mm in range(len(maparr)):
        if full_sky:
            map_alm = H.map2alm(maparr[mm], lmax = lmax)
            curr_weight = weightsarr[mm][:lmax]            
            map_alm_weighted = H.almxfl(map_alm, curr_weight)
        else:
            curr_map = maparr[mm]
            weightsarr[mm][np.isnan(weightsarr[mm])]=0.
            weightsarr[mm][np.isinf(weightsarr[mm])]=0.
            map_alm_weighted = np.fft.fft2(curr_map) * weightsarr[mm]
        weighted_almarr.append(map_alm_weighted)

    ilc_alm = np.sum(weighted_almarr, axis = 0)

    if full_sky:
        ilc_map = H.alm2map(ilc_alm, nside = nside, verbose = verbose, lmax = lmax)
    else:
        ilc_map = np.fft.ifft2( ilc_alm ).real

    return ilc_map

################################################################################################################
