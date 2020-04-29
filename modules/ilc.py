import numpy as np, sys, os, scipy as sc, healpy as H, foregrounds as fg, misc
from pylab import *
################################################################################################################
def get_analytic_covariance(param_dict, freqarr, nl_dic = None, bl_dic = None, ignore_fg = [], pol = 0, pol_frac_per_cent_dust = 0.02, pol_frac_per_cent_radio = 0.03, pol_frac_per_cent_tsz = 0., pol_frac_per_cent_ksz = 0., include_gal = 0):

    #ignore_fg = foreground terms that must be ignored
    possible_ignore_fg = ['cmb', 'tsz', 'ksz', 'radio', 'dust']
    if len(ignore_fg)>0:
        if not all( [ currfg in possible_ignore_fg for currfg in ignore_fg] ):
            print( '\n\t Alert: Elements of ignore_fg should be one of the following: %s\n\n' %(np.array2string(np.asarray(possible_ignore_fg))) )
            sys.exit()

    el, cl_cmb = fg.get_foreground_power_spt('CMB', freq1 = param_dict['freq0'], freq2 = param_dict['freq0'])
    el, cl_ksz = fg.get_foreground_power_spt('kSZ', freq1 = param_dict['freq0'], freq2 = param_dict['freq0'])
    if pol:
        cl_ksz = cl_ksz * pol_frac_per_cent_ksz**2.

    cl_dic = {}
    cl_ori = np.zeros(len(el))
    for freq1 in freqarr:
        for freq2 in freqarr:

            if (freq2, freq1) in cl_dic:
                cl_dic[(freq1, freq2)] = cl_dic[(freq2, freq1)]
                continue

            #get dust
            el,  cl_dg_po, cl_dg_clus = fg.get_cl_dust(freq1, freq2, freq0 = param_dict['freq0'], fg_model = param_dict['fg_model'], spec_index_dg_po = param_dict['spec_index_dg_po'], spec_index_dg_clus = param_dict['spec_index_dg_clus'], Tcib = param_dict['Tcib'])
            if pol:
                cl_dg_po = cl_dg_po * pol_frac_per_cent_dust**2.
                cl_dg_clus = cl_dg_clus * pol_frac_per_cent_dust**2.

            #get tsz
            el, cl_tsz = fg.get_cl_tsz(freq1, freq2, freq0 = param_dict['freq0'], fg_model = param_dict['fg_model'])
            if pol:
                cl_tsz = cl_tsz * pol_frac_per_cent_tsz**2.

            #get radio
            el, cl_radio = fg.get_cl_radio(freq1, freq2, freq0 = param_dict['freq0'], fg_model = param_dict['fg_model'], spec_index_rg = param_dict['spec_index_rg'])
            if pol:
                cl_radio = cl_radio * pol_frac_per_cent_radio**2.

            cl = np.copy( cl_ori )
            if 'cmb' not in ignore_fg:
                cl = cl + np.copy(cl_cmb)
            if 'ksz' not in ignore_fg:
                #print('ksz')
                cl = cl + cl_ksz             
            if 'tsz' not in ignore_fg:
                #print('tsz')
                cl = cl + cl_tsz
            if 'radio' not in ignore_fg:
                #print('radio')
                cl = cl + cl_radio
            if 'dust' not in ignore_fg:
                #print('dust')
                cl = cl + cl_dg_po + cl_dg_clus

            if include_gal:# and not pol: #get galactic dust and sync

                which_spec = 'TT'
                if pol: which_spec = 'EE'

                el, cl_gal_dust = fg.get_cl_galactic(param_dict, 'dust', freq1, freq2, which_spec, bl_dic = bl_dic, el = el)
                el, cl_gal_sync = fg.get_cl_galactic(param_dict, 'sync', freq1, freq2, which_spec, bl_dic = bl_dic, el = el)

                if (0):#not pol:

                    loglog(cl_dg_po + cl_dg_clus, label = r'EG: %s, %s' %(freq1, freq2));
                    loglog(cl_gal_dust, label = r'Dust: %s, %s' %(freq1, freq2));title(r'%s' %(param_dict['which_gal_mask']))
                    loglog(cl_gal_sync, label = r'Sync: %s, %s' %(freq1, freq2));title(r'%s' %(param_dict['which_gal_mask']))
                    #cl = cl + cl_gal_dust + cl_gal_sync
                    #loglog(cl);
                    legend(loc = 1)
                    ylim(1e-5, 1e3);#show();sys.exit()
                    xlim(20,param_dict['lmax']);ylim(1e-8,1e6);
                    show();sys.exit()
                
                cl = cl + cl_gal_dust
                cl = cl + cl_gal_sync

            #make sure cl start from el=0 rather than el=10 which is the default for SPT G15 results
            lmin = min(el)
            cl = np.concatenate( (np.zeros(lmin), cl) )

            #noise auto power spectrum
            if nl_dic is not None:

                if (freq1,freq2) in nl_dic:
                    nl = nl_dic[(freq1,freq2)]
                else:
                    nl = nl_dic[freq1]
                    if freq1 != freq2: 
                        nl = np.copy(nl) * 0.
                
                if len(cl) > len(nl):
                    cl = cl[:len(nl)]
                elif len(cl) < len(nl):
                    nl = nl[:len(cl)]

                el = np.arange(len(cl))

            else:
                nl = np.zeros(len(cl))
                #20191116 - fix this: there must be noise correlation in case of atmospheric noise
                print('\n\n\t\tfix me: there must be noise correlation in case of atmospheric noise')
                sys.exit()

            cl = cl + nl

            cl[np.isnan(cl)] = 0.
            cl[np.isinf(cl)] = 0.

            ##loglog(cl); title('%s - %s' %(freq1, freq2)); show()

            cl_dic[(freq1, freq2)] = cl


    return el, cl_dic  

################################################################################################################
def create_clmat(freqarr, elcnt, cl_dic):
    """
    freqarr  = array of frequency channel
    elcnt = \el index
    cl_dic = cl_cmb + cl_FG auto and cross spectra fr thefrequency channel
    """
    nc = len(freqarr)
    clmat = np.zeros( (nc, nc) )
    for ncnt1, freq1 in enumerate(freqarr):
        for ncnt2, freq2 in enumerate(freqarr):
            clmat[ncnt2, ncnt1] = cl_dic[(freq1, freq2)][elcnt]
    return clmat

################################################################################################################
def get_clinv(freqarr, elcnt, cl_dic):
    clmat = np.mat( create_clmat(freqarr, elcnt, cl_dic) )
    clinv = sc.linalg.pinv2(clmat)
    #clinv = sc.linalg.inv(clmat)

    return clinv

################################################################################################################
def residual_power(param_dict, freqarr, el, cl_dic, final_comp = 'cmb', freqcalib_fac = None, lmin = 0, return_weights = 0):

    try:
        lmin = param_dict['lmin']
    except:
        pass

    acap = get_acap(freqarr, final_comp = final_comp, freqcalib_fac = freqcalib_fac)
    nc = len(freqarr)

    cl_residual = np.zeros( (len(el)) )
    weightsarr = np.zeros( (nc, len( el ) ) )
    for elcnt, el in enumerate(el):
        if el <= lmin: continue ## or el>=lmax: continue
        clinv = get_clinv( freqarr, elcnt, cl_dic )
        
        nr = np.dot(clinv, acap)
        dr = np.dot( acap.T, np.dot(clinv, acap) )

        #ILC residuals
        cl_residual[elcnt] = np.asarray(1./dr).squeeze()

        #weights
        weightsarr[:, elcnt] = np.asarray(nr/dr).squeeze()

    cl_residual[np.isinf(cl_residual)] = 0.
    cl_residual[np.isnan(cl_residual)] = 0.
    
    if return_weights:
        return cl_residual, weightsarr
    else:
        return cl_residual

################################################################################################################

def coth(x):
    return (np.exp(x) + np.exp(-x)) / (np.exp(x) - np.exp(-x))

################################################################################################################

def compton_y_to_delta_Tcmb(nu, Tcmb = 2.73, h=6.62607004e-34, k_B=1.38064852e-23):

    if nu<1e3: nu *= 1e9

    x = h * nu / k_B / Tcmb
    g_nu = x * coth(x/2.) - 4.

    return Tcmb * np.mean(g_nu)

################################################################################################################

def get_acap(freqarr, final_comp = 'cmb', freqcalib_fac = None):

    nc = len(freqarr)

    if freqcalib_fac is None: freqcalib_fac = np.ones(nc)

    if final_comp.lower() == 'cmb':
        freqscale_fac = np.ones(nc)

    elif final_comp.lower() == 'tsz' or final_comp.lower() == 'y':

        freqscale_fac = []
        for freq in sorted( freqarr ):
            freqscale_fac.append( compton_y_to_delta_Tcmb(freq * 1e9) )

        freqscale_fac = np.asarray( freqscale_fac )

        if final_comp.lower() == 'tsz': #tsz at 150 GHz
            freqscale_fac = freqscale_fac / freqscale_fac[1]

    acap = np.zeros(nc) + (freqscale_fac * freqcalib_fac) #assuming CMB is the same and calibrations factors are same for all channel
    acap = np.mat(acap).T #should be nc x 1

    return acap

################################################################################################################

def get_ilc_map(final_comp, el, map_dic, bl_dic, nside, lmax, cl_dic = None, nl_dic = None, lmin = 10, freqcalib_fac = None, ignore_fg = [], full_sky = 1, estimate_covariance = 0, mapparams = None, apod_mask = None):

    """
    inputs:
    final_comp: 'CMB' or 'tsz' or 'y'

    map_dic: dicionary containing
    map_dic[95] = map_95 where map_95 is either a healpix map or a SPT3G map object or a simple flat-sky 2D numpy array
    map_dic[150] = map_150 
    map_dic[220] = map_220

    bl_dic = {} #dictionary containing Bl of different frequencies - either a 1d array for healpix maps or 2D array for SPT3G/flat-sky maps
    bl_dic[95], bl_dic[150], bl_dic[220]

    nside = healpix nside set to None for flat-sky maps

    lmax = lmax used for ILC (covariance)

    cl_dic = same as beam dic and contains Cl of different frequencies - either a 1d array for healpix maps or 2D array for SPT3G/flat-sky maps
    nl_dic = same as beam dic and contains Nl of different frequencies - either a 1d array for healpix maps or 2D array for SPT3G/flat-sky maps

    lmin = minimum \el to use

    freqcalib_fac = array containing calibration factors for different frequencies. for example: freqcalib_fac = [1., 1., 1.]

    """
    spt3g_maps = 0
    freqarr, maparr = [], []
    for keyname in sorted(map_dic):
        freqarr.append( keyname )
        curr_map = map_dic[keyname]
        '''
        if isinstance(curr_map, core.G3Frame):
            spt3g_maps = 1
            #get the apodisation mask before converting the map inot an array
            if str(apod_mask) == 'from_weight':
                apod_mask = apodmask.make_border_apodization(curr_map, radius_arcmin=90.0)
            curr_map = np.array( curr_map['T'] ) ## / core.G3Units.uK
            curr_map[np.isnan(curr_map)] = 0.
            curr_map[np.isinf(curr_map)] = 0.
        '''
        if apod_mask is not None:
            curr_map = curr_map * apod_mask

        maparr.append( curr_map )

    #get covariance
    if cl_dic is None:
        '''
        if spt3g_maps:
            el, cl_dic = get_spt3g_covariance_dic(map_dic, lmin, lmax)
        else:
            el, cl_dic = get_analytic_covariance(freqarr, nl_dic = nl_dic, ignore_fg = ignore_fg)
        '''
        if estimate_covariance:
            el, cl_dic = get_map_covariance(map_dic, lmin, lmax) 
        else:
            el, cl_dic = get_analytic_covariance(freqarr, nl_dic = nl_dic, ignore_fg = ignore_fg)


    #get weights
    weightsarr = get_multipole_weightsarr(final_comp, freqarr, el, cl_dic, lmin, freqcalib_fac)#, ignore_fg)
    weightsarr_1d = np.copy(weightsarr)

    #convert weights to 2D if flat-sky
    if not full_sky:
        assert mapparams is not None
        weightsarr_2D = []
        for currW in weightsarr:
            el = np.arange(len(currW))
            currW_2D = flatsky.cl_to_cl2d(el, currW, mapparams) 
            weightsarr_2D.append(currW_2D)
        weightsarr = np.asarray( weightsarr_2D )


    #rebeaming
    rebeamarr = misc.rebeam( bl_dic )

    #modify weights to include rebeam
    weightsarr = weightsarr * rebeamarr
    '''
    plot(weightsarr[0], 'k-'); plot(weightsarr[1], 'r-'); plot(weightsarr[2], 'g-')
    weightsarr = weightsarr * rebeamarr
    plot(weightsarr[0], 'k--'); plot(weightsarr[1], 'r-'); plot(weightsarr[2], 'g-')
    show(); sys.exit()
    sys.exit()
    '''

    #get ilc map now
    ilc_map = apply_ilc_weightsarr(maparr, weightsarr, nside, lmax, full_sky = full_sky)

    if apod_mask is not None:
        ilc_map = ilc_map * apod_mask

    return ilc_map, weightsarr_1d

################################################################################################################

def apply_ilc_weightsarr(maparr, weightsarr, nside, lmax, full_sky = 0, verbose = 0):

    '''
    clf()
    freqs = [95, 150]#, 220]
    colordic = {95: 'darkblue', 150: 'green', 220: 'darkred'}
    for frqcntr, freq in enumerate( freqs ):
        plot(weightsarr[frqcntr], color = colordic[freq], label = r'%s' %(freq))
    plot(np.sum(weightsarr, axis = 0), 'k', label = r'Sum')
    legend(loc = 1);show(); sys.exit()
    '''

    #get the ilc combined map now
    weighted_almarr = []
    for mm in range(len(maparr)):
        if full_sky:
            map_alm = H.map2alm(maparr[mm], lmax = lmax)
            curr_weight = weightsarr[mm][:lmax]
            
            map_alm_weighted = H.almxfl(map_alm, curr_weight)

            '''
            map_alm_weighted = np.zeros_like(map_alm)
            for el in range(len(curr_weight)):
                alm_inds = H.Alm.getidx(lmax, el, np.arange(el + 1))
                map_alm_weighted[alm_inds] =  curr_weight[el] * map_alm[alm_inds]
            '''
            #plot(map_alm_weighted.real, 'r');show();sys.exit()
            ###map_weighted = H.alm2map(map_alm_weighted, nside = nside, verbose = verbose, lmax = lmax)
            #clf();plot(curr_weight); show()
            #H.mollview(maparr[mm], sub = (1,2,1)); H.mollview(map_weighted, sub = (1,2,2)); show()
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

def get_multipole_weightsarr(final_comp, freqarr, el, cl_dic, lmin, freqcalib_fac):#, ignore_fg):

    acap = get_acap(freqarr, final_comp = final_comp, freqcalib_fac = freqcalib_fac)

    nc = len(freqarr)

    #get weightsarr
    if np.ndim(el) == 1:
        weightsarr = np.zeros( (nc, len( el ) ) )
        for elcnt, curr_el in enumerate( el ):
            if curr_el <= lmin: continue ## or el>=lmax: continue
            clinv = get_clinv( freqarr, elcnt, cl_dic )
            
            nr = np.dot(clinv, acap)
            dr = np.dot( acap.T, np.dot(clinv, acap) )

            weight_mat = np.asarray(nr/dr).squeeze()

            weightsarr[:, elcnt] = weight_mat
    else:
        #2d- PSDs
        from IPython import embed; embed()
        sys.exit()
        lx, ly = el
        for elcnt1, curr_lx in enumerate( lx ):
            for elcnt2, curr_ly in enumerate( ly ):
                curr_el = np.sqrt(curr_lx**2. + curr_ly**2.)
                if curr_el <= lmin: continue ## or el>=lmax: continue
                clinv = get_clinv( freqarr, elcnt, cl_dic )
            
            nr = np.dot(clinv, acap)
            dr = np.dot( acap.T, np.dot(clinv, acap) )

            weight_mat = np.asarray(nr/dr).squeeze()

            weightsarr[:, elcnt] = weight_mat

    return weightsarr

################################################################################################################
def get_map_covariance(map_dic, lmax, bl_dic, apod_mask = None):#, lmin = 2):

    print('Estimating covariance from maps now')

    freqarr = sorted(map_dic.keys())
    cl_dic = {}
    for cntr1, freq1 in enumerate( freqarr ):
        for cntr2, freq2 in enumerate( freqarr ):
            if (freq2, freq1) in cl_dic:
                cl_dic[(freq1, freq2)] = cl_dic[(freq2, freq1)]
                continue
            print((freq2, freq1))
            map1, map2 = map_dic[freq1], map_dic[freq2]

            cl = H.anafast(map1, map2, lmax = lmax - 1)# lmax + lmin - 1)
            blsquare = bl_dic[freq1] * bl_dic[freq2]
            cl /= blsquare
            #cl = cl[lmin: lmax]

            #account for the mask
            if apod_mask is not None:
                fsky = np.mean(apod_mask)**2.
                cl /= fsky

            cl[np.isnan(cl)] = 0.
            cl[np.isinf(cl)] = 0.
            cl_dic[(freq1, freq2)] = cl 
        
            #el = np.arange(lmin, lmax + lmin)
            el = np.arange(lmax)

    return el, cl_dic

################################################################################################################
'''
def get_spt3g_covariance_dic(map_dic, lmin, lmax, apod_mask = 'from_weight', return_2d = 1):

    freqarr = sorted(map_dic.keys())
    cl_dic = {}
    for cntr1, freq1 in enumerate( freqarr ):
        for cntr2, freq2 in enumerate( freqarr ):
            if (freq2, freq1) in cl_dic:
                cl_dic[(freq1, freq2)] = cl_dic[(freq2, freq1)]
                continue
            map1, map2 = map_dic[freq1], map_dic[freq2]

            if return_2d:
                mapps = map_analysis.calculate_powerspectra(map1, map2, apod_mask = apod_mask, return_2d = 1)['TT']

                mapps[np.isnan(mapps)] = 0.
                mapps[np.isinf(mapps)] = 0.
                cl_dic[(freq1, freq2)] = mapps 

                from IPython import embed; embed()

                el = mapps.get_lxly()

            else:
                cl = map_analysis.calculate_powerspectra(map1, map2, lmin = lmin, lmax = lmax, apod_mask = apod_mask)['TT']
                cl = np.concatenate( (np.zeros(lmin), cl) )
                cl[np.isnan(cl)] = 0.
                cl[np.isinf(cl)] = 0.
                cl_dic[(freq1, freq2)] = cl 
            
                el = np.arange(lmin, lmax)
                el = np.concatenate( (np.zeros(lmin), el) )

    return el, cl_dic
'''

################################################################################################################
