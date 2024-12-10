import numpy as np, re, healpy as H, glob, sys, os, scipy as sc
import flatsky, tools
from pylab import *


def get_pixel_pixel_covariance(mapparams, el, cl_signal, cl_noise = None, bl = None, noofsims = 2500, binsize_am = 0.5, maxbin_am = 10.):
    
    nx, ny, dx = mapparams
    boxsize_am = nx * dx
    x1, x2 = -boxsize_am/2, boxsize_am/2
    ra = dec = np.linspace(x1, x2, nx)
    ra_grid, dec_grid = np.meshgrid( ra, dec )

    rad_prf_arr = []
    for n in range(noofsims):

        #signal
        curr_signal_map = flatsky.make_gaussian_realisation(mapparams, el, cl_signal, bl = bl)

        #noise
        curr_noise_map = flatsky.make_gaussian_realisation(mapparams, el, cl_noise)

        #total
        curr_sim_map = curr_signal_map + curr_noise_map



        #imshow(curr_sim_map); colorbar(); show(); sys.exit()
        rad_prf_bins, curr_rad_prf = flatsky.radial_profile(curr_sim_map, binsize=binsize_am, maxbin=maxbin_am, xy = (ra_grid, dec_grid))
        #print(curr_rad_prf); sys.exit()
        rad_prf_arr.append( curr_rad_prf )
    rad_prf_arr = np.asarray(rad_prf_arr)

    totbins = np.shape(rad_prf_arr)[1]
    rad_prf_mean = np.mean(rad_prf_arr, axis = 0)
    rad_prf_arr = rad_prf_arr - rad_prf_mean
    #sim_cov = sims.calcCov(rad_prf_arr, noofsims, npixels = totbins)
    sim_cov = np.cov(rad_prf_arr.T)
    
    return rad_prf_bins, sim_cov

def inject_source_at_centre(mapparams, source_pol_flux_in_uk, bl = None):
    nx, ny, dx = mapparams
    source_sim_map = np.zeros( (ny, nx) )
    s = int(nx/2)
    e = s + 1
    source_sim_map[s:e, s:e] = source_pol_flux_in_uk

    if bl is not None:
        bl_2d = flatsky.cl_to_cl2d(el, bl, mapparams)
        source_sim_map = np.fft.ifft2( np.fft.fft2( source_sim_map ) * bl_2d ).real

    return source_sim_map

def get_likelihood(data, model, cov):

    """
    function to calculate the likelihood given data, model, and the covariance matrix.
    """

    cov = np.mat(cov)
    cinv = sc.linalg.pinv2(cov)
    #sign, logdetval = np.linalg.slogdet(cov)
    #logdetval = logdetval * sign

    d = data.flatten()
    m = model.flatten()## - np.mean(MODEL.flatten())
    dprime = d-m

    loglval =  -0.5 * np.asarray( np.dot(dprime.T, np.dot( cinv, dprime ))).squeeze()

    return loglval


def likelihood_finer_resol(M, L, intrp_type=2):

    import scipy.optimize as optimize

    deltaM=np.diff(M)[0]
    M_ip=np.arange(min(M),max(M),deltaM/100.)

    if intrp_type == 2: #Guassian fitting
        #first guess a good parameter
        Mfit=M[np.argmax(L)]
        gau_width=abs(Mfit - M[np.argmin(abs(L))])#/2.35 * 2.
        p0=[0.,np.max(L),Mfit,gau_width]
        p1, success=optimize.leastsq(fitting_func_gaussian, p0, args=(p0, M, L))

        L_ip=fitting_func_gaussian(p1, p1, M_ip, return_fit=1)

    return M_ip, L_ip


def lnlike_to_like(M, lnlike, intrp_type=1):

    lnlike=lnlike - max(lnlike)

    '''
    if intrp_type == 1:
        deltaM=np.diff(M)[0]
        M_ip=np.arange(min(M),max(M),deltaM/100.)
        lnlike_ip=np.interp(M_ip, M, lnlike)
        M=np.copy(M_ip)
        lnlike=np.copy(lnlike_ip)
    '''

    delta_chisq=max(lnlike) - lnlike[0]
    snr=np.sqrt(2 * delta_chisq)

    L=np.exp(lnlike); L/=max(L)
    recov_mass=M[np.argmax(L)]

    if intrp_type <= 1: #no interpolation
        return M, L, recov_mass, snr

    #from IPython import embed; embed()
    M_ip, L_ip=likelihood_finer_resol(M, L, intrp_type=intrp_type)
    L_ip /= max(L_ip)
    recov_mass=M_ip[np.argmax(L_ip)]

    return M_ip, L_ip, recov_mass, snr
def get_nl(noiseval, el, elknee = -1, alphaknee = 0, beamval = None):

    delta_T_radians = noiseval * np.radians(1./60.)
    nl = np.tile(delta_T_radians**2., int(max(el)) + 1 )
    nl = np.asarray( [nl[int(l)] for l in el] )
    nl_white = np.copy(nl)

    if beamval is not None:
        bl = flatsky.gauss_beam(beamval, max(el))
        nl = nl/(bl**2.)

    if elknee != -1.:
        nl = np.copy(nl) * (1. + (elknee * 1./el)**alphaknee )

    nl[np.isnan(nl) | np.isinf(nl)] = 0.

    return nl
