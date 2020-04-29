import numpy as np, os, sys, healpy as H
import pandas as pd

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
def get_beam_dic(freqs, beam_noise_dic, lmax, opbeam = None, make_2d = 0, mapparams = None):
    bl_dic =  {}
    for freq in freqs:
        beamval, noiseval = beam_noise_dic[freq]
        bl_dic[freq] = H.gauss_beam(np.radians(beamval/60.), lmax=lmax-1)

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
        if freq is 'effective': continue
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

    fwhm_radians = np.radians(beamval/60.)
    sigma = fwhm_radians / np.sqrt(8. * np.log(2.))
    sigma2 = sigma ** 2
    bl = np.exp(el * (el+1) * sigma2)

    return bl

################################################################################################################

def get_nl(noiseval, el, beamval, use_beam_window = 1, uk_to_K = 0, elknee = -1, alphaknee = 0, beamval2 = None, noiseval2 = None, elknee2 = -1, alphaknee2 = 0, rho = None):

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

    if cross_band_noise:
        delta_T2_radians = noiseval2 * np.radians(1./60.)
        nl2 = np.tile(delta_T2_radians**2., int(max(el)) + 1 )
        nl2 = np.asarray( [nl2[int(l)] for l in el] )

    if use_beam_window: 
        nl *= bl
        if cross_band_noise: nl2 *= bl2

    if elknee != -1.:
        nl = np.copy(nl) * (1. + (elknee * 1./el)**alphaknee )
        if cross_band_noise and elknee2 != -1.:
            nl2 = np.copy(nl2) * (1. + (elknee2 * 1./el)**alphaknee2 )

    if cross_band_noise:
        final_nl = rho * nl**0.5 * nl2**0.5
    else:
        final_nl = np.copy(nl)

    return final_nl

################################################################################################################
def get_delta_cl(el, cl, nl, fsky = 1., delta_l = 1.):

    delta_cl = np.sqrt(2./ (2.*el + 1) / fsky / delta_l) * (cl)## + nl)

    return delta_cl
################################################################################################################
