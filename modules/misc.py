import numpy as np, os, sys#, healpy as H
import pandas as pd
import flatsky

#################################################################################
#################################################################################
#################################################################################

def fn_get_param_dict(paramfile):
    params, paramvals = np.genfromtxt(
        paramfile, delimiter = '=', unpack = True, autostrip = True, dtype='unicode')
    param_dict = {}
    for p,pval in zip(params,paramvals):
        if pval in ['T', 'True']:
            pval = True
        elif pval in ['F', 'False']:
            pval = False
        elif pval == 'None':
            pval = None
        else:
            try:
                pval = float(pval)
                if int(pval) == float(pval):
                    pval = int(pval)
            except:
                pass
        # replace unallowed characters in paramname
        p = p.replace('(','').replace(')','')
        param_dict[p] = pval
    return param_dict

################################################################################################################

def gauss_beam(fwhm, lmax=512, pol=False):
    """Gaussian beam window function

    Computes the spherical transform of an axisimmetric gaussian beam

    For a sky of underlying power spectrum C(l) observed with beam of
    given FWHM, the measured power spectrum will be
    C(l)_meas = C(l) B(l)^2
    where B(l) is given by gaussbeam(Fwhm,Lmax).
    The polarization beam is also provided (when pol = True ) assuming
    a perfectly co-polarized beam
    (e.g., Challinor et al 2000, astro-ph/0008228)

    Parameters
    ----------
    fwhm : float
        full width half max in radians
    lmax : integer
        ell max
    pol : bool
        if False, output has size (lmax+1) and is temperature beam
        if True output has size (lmax+1, 4) with components:
        * temperature beam
        * grad/electric polarization beam
        * curl/magnetic polarization beam
        * temperature * grad beam

    Returns
    -------
    beam : array
        beam window function [0, lmax] if dim not specified
        otherwise (lmax+1, 4) contains polarized beam
    """

    sigma = fwhm / np.sqrt(8.0 * np.log(2.0))
    ell = np.arange(lmax + 1)
    sigma2 = sigma ** 2
    g = np.exp(-0.5 * ell * (ell + 1) * sigma2)

    if not pol:  # temperature-only beam
        return g
    else:  # polarization beam
        # polarization factors [1, 2 sigma^2, 2 sigma^2, sigma^2]
        pol_factor = np.exp([0.0, 2 * sigma2, 2 * sigma2, sigma2])
        return g[:, np.newaxis] * pol_factor


def get_beam_dic(freqs, beam_noise_dic, lmax, opbeam = None, make_2d = 0, mapparams = None):
    bl_dic =  {}
    for freq in freqs:
        beamval, noiseval = beam_noise_dic[freq]
        #bl_dic[freq] = H.gauss_beam(np.radians(beamval/60.), lmax=lmax-1)
        bl_dic[freq] = gauss_beam(np.radians(beamval/60.), lmax=lmax-1)

        if make_2d:
            assert mapparams is not None
            el = np.arange(len(bl_dic[freq]))
            bl_dic[freq] = flatsky.cl_to_cl2d(el, bl_dic[freq], mapparams) 

    if opbeam is not None:
        bl_dic['effective'] = H.gauss_beam(np.radians(opbeam/60.), lmax=lmax-1)

        if make_2d:
            assert mapparams is not None
            bl_dic['effective'] = flatsky.cl_to_cl2d(el, bl_dic['effective'], mapparams) 

    return bl_dic
################################################################################################################

def rebeam(bl_dic, threshold = 1000.):
    #freqarr = sorted( list(bl_dic.keys()) )
    freqarr = []
    for nu in list(bl_dic.keys()): 
        if isinstance(nu, int):
            freqarr.append(nu)
    freqarr = sorted(freqarr)

    bl_eff = bl_dic['effective']
    rebeamarr = []
    for freq in freqarr:
        if freq == 'effective': continue
        bad_inds = np.where(bl_dic[freq]<0)
        bl_dic[freq][bad_inds] = 0.
        currinvbeamval = 1./bl_dic[freq]
        currinvbeamval[currinvbeamval>threshold] = threshold
        rebeamval = bl_eff * currinvbeamval
        rebeamarr.append( rebeamval )

    return np.asarray( rebeamarr )

################################################################################################################

def healpix_rotate_coords(hmap, coord):
    """
    coord = ['C', 'G'] to convert a map in RADEC to Gal.    
    """

    #get map pixel
    pixel = np.arange(len(hmap))

    #get angles in this map first
    nside = H.get_nside(hmap)
    angles = H.pix2ang(nside, pixel)

    #roate the angles to the desired new coordinate
    rotated_angles = H.Rotator(coord=coord)(*angles)

    #get the rotated pixel values
    rotated_pixel = H.ang2pix(nside, *rotated_angles)

    #initialise new map
    rot_hmap = np.zeros(len(pixel))

    #push the original map pixel to the new map (in the rotated pixel positions)
    rot_hmap[rotated_pixel] = hmap[pixel]

    return rot_hmap

################################################################################################################

def get_bl(beamval, el):

    """
    this funciton returns bl**2
    """

    fwhm_radians = np.radians(beamval/60.)
    sigma = fwhm_radians / np.sqrt(8. * np.log(2.))
    sigma2 = sigma ** 2
    bl = np.exp(el * (el+1) * sigma2)

    #bl = np.exp(-0.5*el * (el+1) * sigma2) ##traditionally this is how it is written. But note that I do the inverse and this funciton returns bl**2
    return bl

################################################################################################################

def get_nl(noiseval, el, beamval, use_beam_window = 1, uk_to_K = 0, elknee = -1, alphaknee = 0, beamval2 = None, noiseval2 = None, elknee2 = -1, alphaknee2 = 0, rho = None, Nred1 = -1., Nred2=-1., so_like = False):

    if Nred1!= -1:
        total_years = 5.
        fsky = 0.4 #0.35
        survey_time = 1.
        obs_efficiency = 0.2
        noisy_map_eges_ign_factor = 0.15
        single_year = 365.25 * 24. * 3600. * obs_efficiency * (1.-noisy_map_eges_ign_factor)
        sky_area = 4. * np.pi * fsky
        year_scaling = single_year / total_years
        Nred1 = Nred1  * sky_area / year_scaling
        if Nred2 != -1.:
            Nred2 = Nred2  * sky_area / year_scaling

    cross_band_noise = 0
    if noiseval2 is not None and beamval2 is not None:
        assert rho is not None
        cross_band_noise = 1

    if uk_to_K: 
        noiseval = noiseval/1e6
        if cross_band_noise: noiseval2 = noiseval2/1e6

    if use_beam_window:
        bl = get_bl(beamval, el)
        if cross_band_noise: bl2 = get_bl(beamval2, el)

    delta_T_radians = noiseval * np.radians(1./60.)
    nl = np.tile(delta_T_radians**2., int(max(el)) + 1 )
    nl = np.asarray( [nl[int(l)] for l in el] )
    nl_white = np.copy(nl)

    if cross_band_noise:
        delta_T2_radians = noiseval2 * np.radians(1./60.)
        nl2 = np.tile(delta_T2_radians**2., int(max(el)) + 1 )
        nl2 = np.asarray( [nl2[int(l)] for l in el] )
        nl2_white = np.copy(nl2)

    if use_beam_window: 
        nl *= bl
        if cross_band_noise: nl2 *= bl2

    if elknee != -1.:
        if Nred1==-1:
            nl = np.copy(nl) * (1. + (elknee * 1./el)**alphaknee )
        else:
            nl = np.copy(nl) + Nred1*(elknee * 1./el)**alphaknee
            if cross_band_noise and elknee2 != -1.:
                if Nred2==-1:
                    nl2 = np.copy(nl2) * (1. + (elknee2 * 1./el)**alphaknee2 )
                else:
                    nl2 = np.copy(nl2) + Nred2*(elknee2 * 1./el)**alphaknee2

    if cross_band_noise and (elknee != -1. and elknee2 != -1.):
        ###final_nl = rho * nl**0.5 * nl2**0.5
        final_nl = rho * delta_T_radians * (elknee * 1./el)**(alphaknee/2.) * delta_T2_radians * (elknee2 * 1./el)**(alphaknee2/2.)
        #N[i,j,:] = rho * (w1*np.pi/180./60. * (ell/knee1)**(gamma1/2)) * (w2*np.pi/180./60. * (ell/knee2)**(gamma2/2))
    else:
        final_nl = np.copy(nl)

    return final_nl

def get_nl_v1(noiseval, el, beamval, use_beam_window = 1, uk_to_K = 0, elknee = -1, alphaknee = 0, beamval2 = None, noiseval2 = None, elknee2 = -1, alphaknee2 = 0, rho = None, Nred1 = -1., Nred2=-1.):

    cross_band_noise = 0
    if noiseval2 is not None and beamval2 is not None:
        assert rho is not None
        cross_band_noise = 1

    if uk_to_K: 
        noiseval = noiseval/1e6
        if cross_band_noise: noiseval2 = noiseval2/1e6

    if use_beam_window:
        bl = get_bl(beamval, el)
        if cross_band_noise: bl2 = get_bl(beamval2, el)

    delta_T_radians = noiseval * np.radians(1./60.)
    nl = np.tile(delta_T_radians**2., int(max(el)) + 1 )
    nl = np.asarray( [nl[int(l)] for l in el] )
    nl_white = np.copy(nl)

    if cross_band_noise:
        delta_T2_radians = noiseval2 * np.radians(1./60.)
        nl2 = np.tile(delta_T2_radians**2., int(max(el)) + 1 )
        nl2 = np.asarray( [nl2[int(l)] for l in el] )
        nl2_white = np.copy(nl2)

    if use_beam_window: 
        nl *= bl
        if cross_band_noise: nl2 *= bl2

    if elknee != -1.:
        if Nred1==-1:
            nl = np.copy(nl) * (1. + (elknee * 1./el)**alphaknee )
        else:
            nl = np.copy(nl) + Nred1*(elknee * 1./el)**alphaknee
            if cross_band_noise and elknee2 != -1.:
                if Nred2==-1:
                    nl2 = np.copy(nl2) * (1. + (elknee2 * 1./el)**alphaknee2 )
                else:
                    nl2 = np.copy(nl2) + Nred2*(elknee2 * 1./el)**alphaknee2

    if cross_band_noise and (elknee != -1. and elknee2 != -1.):
        ###final_nl = rho * nl**0.5 * nl2**0.5
        final_nl = rho * delta_T_radians * (elknee * 1./el)**(alphaknee/2.) * delta_T2_radians * (elknee2 * 1./el)**(alphaknee2/2.)
        #N[i,j,:] = rho * (w1*np.pi/180./60. * (ell/knee1)**(gamma1/2)) * (w2*np.pi/180./60. * (ell/knee2)**(gamma2/2))
    else:
        final_nl = np.copy(nl)

    return final_nl

################################################################################################################
def get_delta_cl(el, cl, nl, fsky = 1., delta_l = 1.):

    delta_cl = np.sqrt(2./ (2.*el + 1) / fsky / delta_l) * (cl)## + nl)

    return delta_cl
################################################################################################################

def get_apod_mask(ra_grid, dec_grid, mask_radius = 2., taper_radius = 6., in_arcmins = 1):

    import scipy as sc
    import scipy.ndimage as ndimage

    if not in_arcmins:
        ra_grid_arcmins = ra_grid * 60.
        dec_grid_arcmins = dec_grid * 60.
    else:
        ra_grid_arcmins = np.copy( ra_grid )
        dec_grid_arcmins = np.copy( dec_grid ) 

    radius = np.sqrt( (ra_grid_arcmins**2. + dec_grid_arcmins**2.) )

    mask = np.zeros( ra_grid_arcmins.shape )
    if (1): #20180118
        inds_to_mask = np.where((radius<=mask_radius)) #2arcmins - fix this for now
        mask[inds_to_mask[0], inds_to_mask[1]] = 1.

    ker=np.hanning(taper_radius)
    ker2d=np.asarray( np.sqrt(np.outer(ker,ker)) )

    mask=ndimage.convolve(mask, ker2d)
    mask/=mask.max()

    return mask

################################################################################################################
