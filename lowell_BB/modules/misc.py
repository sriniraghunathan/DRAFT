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
        if isinstance(pval, str):
            if pval.find(',')>-1 and pval.find('[')>-1 and pval.find(']')>-1: #list
                pval = pval.replace('[','').replace(']','').replace(', ',',').replace('\'','')
                pval = pval.split(',')
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
        bl_dic[freq] = gauss_beam(np.radians(beamval/60.), lmax=lmax-1) #simply does H.gauss_beam but healpy stopped working after updating my python. 

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

def get_nl(noiseval, el, beamval, use_beam_window = 1, uk_to_K = 0, elknee = -1, alphaknee = 0, beamval2 = None, noiseval2 = None, elknee2 = -1, alphaknee2 = 0, rho = None):

    cross_band_noise = 0
    if noiseval2 is not None and beamval2 is not None:
        assert rho is not None
        cross_band_noise = 1

    if uk_to_K:
        noiseval = noiseval/1e6
        if cross_band_noise: noiseval2 = noiseval2/1e6

    if use_beam_window:
        bl = gauss_beam(np.radians(beamval/60.), lmax = max(el))
        if cross_band_noise: bl2 = gauss_beam(np.radians(beamval2/60.), lmax = max(el))

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
        nl /= bl**2.
        if cross_band_noise: nl2 /= bl2**2.

    if elknee != -1.:
        nl = np.copy(nl) * (1. + (elknee * 1./el)**alphaknee )
        if cross_band_noise and elknee2 != -1.:
            nl2 = np.copy(nl2) * (1. + (elknee2 * 1./el)**alphaknee2 )

    if cross_band_noise:
        final_nl = rho * delta_T_radians * (elknee * 1./el)**(alphaknee/2.) * delta_T2_radians * (elknee2 * 1./el)**(alphaknee2/2.)
    else:
        final_nl = np.copy(nl)

    return final_nl

################################################################################################################
