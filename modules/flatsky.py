r"""
Flat-sky routines for simulating and analyzing small patches of sky in the flat-sky approximation.

The public functions are grouped into four sections:

* Fourier grid: ``get_lxly``, ``get_lxly_az_angle``, ``cl_to_cl2d``
* Filters: ``get_lpf_hpf``, ``wiener_filter``
* Simulations: ``convert_eb_qu``, ``cl2map``, ``make_gaussian_realization``
* Power spectra: ``radial_profile``, ``map2cl``

``cl_to_cl2d`` is called by ``misc.get_beam_dic`` to place 1-D beam transfer functions on a 2-D grid.
The remaining routines are standalone utilities that are currently not used elsewhere in this repository.

A flat-sky patch is described by ``mapparams = [nx, ny, dx, dy]``, where ``nx`` and ``ny`` are the number of pixels along each axis and ``dx`` and ``dy`` are the pixel sizes in arcminutes.
Multipoles follow the flat-sky convention :math:`\ell = \sqrt{\ell_x^2 + \ell_y^2}`, with the 2-D Fourier grid supplied by ``get_lxly``.
Maps and power spectra carry whatever units the input ``cl`` is given in.
"""

import numpy as np

#import sys, os, scipy as sc, healpy as H


def cl_to_cl2d(el, cl, mapparams):
    """
    converts 1d_cl to 2d_cl
    inputs:
    el = el values over which cl is defined
    cl = power spectra - cl

    mapparams = [nx, ny, dx, dy] where ny, nx = flatskymap.shape; and dy, dx are the pixel resolution in arcminutes.
    for example: [100, 100, 0.5, 0.5] is a 50' x 50' flatskymap that has dimensions 100 x 100 with dx = dy = 0.5 arcminutes.

    output:
    2d_cl
    """
    if np.any(np.diff(el) <= 0):
        raise ValueError('el must be strictly increasing, but it decreases or repeats at index %s' % (int(np.argmax(np.diff(el) <= 0)) + 1))

    lx, ly = get_lxly(mapparams)
    ell = np.sqrt(lx**2. + ly**2.)

    cl2d = np.interp(ell.flatten(), el, cl).reshape(ell.shape)  #NOTE: np.interp clamps, so modes above el.max() get cl[-1] rather than zero

    return cl2d


def get_lxly(mapparams):
    """
    returns lx, ly based on the flatskymap parameters
    input:
    mapparams = [nx, ny, dx, dy] where ny, nx = flatskymap.shape; and dy, dx are the pixel resolution in arcminutes.
    for example: [100, 100, 0.5, 0.5] is a 50' x 50' flatskymap that has dimensions 100 x 100 with dx = dy = 0.5 arcminutes.

    output:
    lx, ly
    """
    nx, ny, dx, dx = mapparams  #NOTE: dy is dropped and dx ends up holding dy
    dx = np.radians(dx/60.)

    lx, ly = np.meshgrid( np.fft.fftfreq( nx, dx ), np.fft.fftfreq( ny, dx ) )
    lx *= 2*np.pi
    ly *= 2*np.pi

    return lx, ly


def get_lxly_az_angle(lx, ly):
    """
    azimuthal angle from lx, ly

    inputs:
    lx, ly = 2d lx and ly arrays

    output:
    azimuthal angle
    """
    return 2*np.arctan2(lx, -ly)


def convert_eb_qu(map1, map2, mapparams, eb_to_qu=1):
    lx, ly = get_lxly(mapparams)
    angle = get_lxly_az_angle(lx, ly)

    map1_fft, map2_fft = np.fft.fft2(map1), np.fft.fft2(map2)
    if eb_to_qu:
        map1_mod = np.fft.ifft2( np.cos(angle) * map1_fft - np.sin(angle) * map2_fft ).real
        map2_mod = np.fft.ifft2( np.sin(angle) * map1_fft + np.cos(angle) * map2_fft ).real
    else:
        map1_mod = np.fft.ifft2( np.cos(angle) * map1_fft + np.sin(angle) * map2_fft ).real
        map2_mod = np.fft.ifft2( -np.sin(angle) * map1_fft + np.cos(angle) * map2_fft ).real

    return map1_mod, map2_mod


def get_lpf_hpf(mapparams, lmin_lmax, filter_type=0):
    """
    filter_type = 0 - low pass filter
    filter_type = 1 - high pass filter
    filter_type = 2 - band pass
    """
    lx, ly = get_lxly(mapparams)
    ell = np.sqrt(lx**2. + ly**2.)
    fft_filter = np.ones(ell.shape)
    if filter_type == 0:
        fft_filter[ell > lmin_lmax] = 0.
    elif filter_type == 1:
        fft_filter[ell < lmin_lmax] = 0.
    elif filter_type == 2:
        lmin, lmax = lmin_lmax
        fft_filter[ell < lmin] = 0.
        fft_filter[ell > lmax] = 0.
    else:
        raise ValueError('filter_type must be 0 (low pass), 1 (high pass) or 2 (band pass), got %s' % (filter_type))

    return fft_filter


def wiener_filter(mapparams, cl_signal, cl_noise, el=None):
    if el is None:
        el = np.arange(len(cl_signal))

    #nx, ny, dx, dx = mapparams  #unused

    #get 2D cl
    cl_signal2d = cl_to_cl2d(el, cl_signal, mapparams)
    cl_noise2d = cl_to_cl2d(el, cl_noise, mapparams)

    wiener_filter2d = cl_signal2d / (cl_signal2d + cl_noise2d)
    wiener_filter2d[np.isnan(wiener_filter2d)] = 0.

    return wiener_filter2d


def cl2map(mapparams, cl, el=None):
    """
    cl2map module - creates a flat sky map based on the flatskymap parameters and the input power spectra

    input:
    mapparams = [nx, ny, dx, dy] where ny, nx = flatskymap.shape; and dy, dx are the pixel resolution in arcminutes.
    for example: [100, 100, 0.5, 0.5] is a 50' x 50' flatskymap that has dimensions 100 x 100 with dx = dy = 0.5 arcminutes.

    cl: input 1d cl power spectra

    el: if None, then computed here.

    output:
    flatskymap with the given map specifications

    """
    if el is None:
        el = np.arange(len(cl))

    nx, ny, dx, dx = mapparams  #NOTE: dy is dropped and dx ends up holding dy

    #get 2D cl
    cl2d = cl_to_cl2d(el, cl, mapparams)

    #pixel area normalization
    dx_rad = np.radians(dx/60.)
    pix_area_norm = np.sqrt(1. / (dx_rad**2.))
    cl2d_sqrt_normed = np.sqrt(cl2d) * pix_area_norm

    #make a random Gaussian realization now
    gauss_reals = np.random.randn(nx, ny)  #NOTE: shape (nx, ny), but the Fourier grid from get_lxly is (ny, nx)

    #convolve with the power spectra
    flatskymap = np.fft.ifft2( np.fft.fft2(gauss_reals) * cl2d_sqrt_normed).real
    flatskymap = flatskymap - np.mean(flatskymap)

    return flatskymap


def map2cl(mapparams, flatskymap1, flatskymap2=None, binsize=None, minbin=100, maxbin=10000, mask=None, filter_2d=None):
    """
    map2cl module - get the power spectra of map/maps

    input:
    mapparams = [nx, ny, dx, dy] where ny, nx = flatskymap.shape; and dy, dx are the pixel resolution in arcminutes.
    for example: [100, 100, 0.5, 0.5] is a 50' x 50' flatskymap that has dimensions 100 x 100 with dx = dy = 0.5 arcminutes.

    flatskymap1: map1 with dimensions (ny, nx)
    flatskymap2: provide map2 with dimensions (ny, nx) cross-spectra

    binsize: el bins. computed automatically if None

    cross_power: if set, then compute the cross power between flatskymap1 and flatskymap2

    output:
    auto/cross power spectra: [el, cl, cl_err]
    """
    nx, ny, dx, dx = mapparams  #NOTE: dy is dropped and dx ends up holding dy
    dx_rad = np.radians(dx/60.)

    lx, ly = get_lxly(mapparams)

    if binsize is None:
        binsize = lx.ravel()[1] - lx.ravel()[0]

    if flatskymap2 is None:
        flatskymap_psd = abs( np.fft.fft2(flatskymap1) * dx_rad)**2 / (nx * ny)
    else: #cross spectra now
        if flatskymap1.shape != flatskymap2.shape:
            raise ValueError('flatskymap1 and flatskymap2 must have the same shape, got %s and %s' % (flatskymap1.shape, flatskymap2.shape))
        flatskymap_psd = np.fft.fft2(flatskymap1) * dx_rad * np.conj( np.fft.fft2(flatskymap2) ) * dx_rad / (nx * ny)

    rad_prf = radial_profile(flatskymap_psd, (lx, ly), bin_size=binsize, minbin=minbin, maxbin=maxbin, to_arcmins=0)
    el, cl = rad_prf[:, 0], rad_prf[:, 1]

    if mask is not None:
        fsky = np.mean(mask)  #NOTE: <W> or <W^2>?
        cl /= fsky

    if filter_2d is not None:
        rad_prf_filter_2d = radial_profile(filter_2d, (lx, ly), bin_size=binsize, minbin=minbin, maxbin=maxbin, to_arcmins=0)
        el, fl = rad_prf_filter_2d[:, 0], rad_prf_filter_2d[:, 1]
        cl /= fl  #NOTE: Unguarded, so bins where the filter is identically zero return inf

    return el, cl


def radial_profile(z, xy=None, bin_size=1., minbin=0., maxbin=10., to_arcmins=1):
    """
    get the radial profile of an image (both real and fourier space)
    """
    z = np.asarray(z)
    if xy is None:
        if z.ndim != 2:
            raise ValueError('z must be 2-D when xy is not given, got %s dimensions' % (z.ndim))
        x, y = np.indices(z.shape)
    else:
        x, y = xy

    #radius = np.hypot(X,Y) * 60.
    radius = (x**2. + y**2.) ** 0.5
    if to_arcmins: radius *= 60.

    binarr = np.arange(minbin, maxbin, bin_size)
    radprf = np.zeros((len(binarr), 3))

    hit_count = []

    for b, binval in enumerate(binarr):
        ind = np.where((radius >= binval) & (radius < binval + bin_size))
        radprf[b, 0] = (binval + bin_size/2.)
        hits = len(np.where(abs(z[ind]) > 0.)[0])  #NOTE: Counts only non-zero pixels

        if hits > 0:
            radprf[b, 1] = np.sum(z[ind]) / hits  #NOTE: radprf is float, so casting may occur
            radprf[b, 2] = np.std(z[ind])
        hit_count.append(hits)

    hit_count = np.asarray(hit_count)
    std_mean = np.sum(radprf[:, 2] * hit_count) / np.sum(hit_count)
    errval = std_mean / (hit_count)**0.5  #NOTE: inf for empty bins and std_mean is pooled across all bins, not per bin
    radprf[:, 2] = errval

    return radprf


def make_gaussian_realization(mapparams, el, cl, cl2=None, cl12=None, bl=None, qu_or_eb='qu'):
    possible_qu_or_eb = ['qu', 'eb']
    if qu_or_eb not in possible_qu_or_eb:
        raise ValueError('qu_or_eb must be one of %s, got %s' % (possible_qu_or_eb, qu_or_eb))

    nx, ny, dx, dy = mapparams
    arcmins2radians = np.radians(1/60.)

    dx *= arcmins2radians
    dy *= arcmins2radians

    #map stuff
    norm = np.sqrt(1. / (dx * dy))

    #1d to 2d now
    cltwod = cl_to_cl2d(el, cl, mapparams)

    if cl2 is not None: #for TE, etc. where two fields are correlated.
        if cl12 is None:
            raise ValueError('cl12 is required when cl2 is given')
        cltwod12 = cl_to_cl2d(el, cl12, mapparams)
        cltwod2 = cl_to_cl2d(el, cl2, mapparams)

    if cl2 is None:
        cltwod = cltwod**0.5 * norm
        cltwod[np.isnan(cltwod)] = 0.

        gauss_reals = np.random.standard_normal([nx, ny])  #NOTE: shape (nx, ny), but the Fourier grid is (ny, nx)
        SIM = np.fft.ifft2( np.copy( cltwod ) * np.fft.fft2( gauss_reals ) ).real

    else: #for TE, etc. where two fields are correlated.
        cltwod12[np.isnan(cltwod12)] = 0.
        cltwod2[np.isnan(cltwod2)] = 0.

        gauss_reals_1 = np.random.standard_normal([nx, ny])  #NOTE: shape (nx, ny), but the Fourier grid is (ny, nx)
        gauss_reals_2 = np.random.standard_normal([nx, ny])  #NOTE: as above.

        '''
        gauss_reals_1 = np.fft.fft2( gauss_reals_1 )
        gauss_reals_2 = np.fft.fft2( gauss_reals_2 )

        t1 = gauss_reals_1 * cltwod12 / cltwod2**0.5
        t2 = gauss_reals_2 * ( cltwod - (cltwod12**2. / cltwod2) )**0.5

        SIM_FFT = (t1 + t2) * norm
        SIM_FFT[np.isnan(SIM_FFT)] = 0.
        SIM = np.fft.ifft2( SIM_FFT ).real
        '''

        gauss_reals_1_fft = np.fft.fft2( gauss_reals_1 )
        gauss_reals_2_fft = np.fft.fft2( gauss_reals_2 )

        #field_1
        cltwod_tmp = np.copy( cltwod )**0.5 * norm
        SIM_FIELD_1 = np.fft.ifft2( cltwod_tmp * gauss_reals_1_fft ).real
        #SIM_FIELD_1 = np.zeros( (ny, nx) )

        #field 2 - has correlation with field_1
        t1 = np.copy( gauss_reals_1_fft ) * cltwod12 / np.copy(cltwod)**0.5
        t2 = np.copy( gauss_reals_2_fft ) * ( cltwod2 - (cltwod12**2. / np.copy(cltwod)) )**0.5
        SIM_FIELD_2_FFT = (t1 + t2) * norm
        SIM_FIELD_2_FFT[np.isnan(SIM_FIELD_2_FFT)] = 0.
        SIM_FIELD_2 = np.fft.ifft2( SIM_FIELD_2_FFT ).real

        #T and E generated. B will simply be zeroes.
        SIM_FIELD_3 = np.zeros( SIM_FIELD_2.shape )
        if qu_or_eb == 'qu': #T, Q, U: convert E/B to Q/U.
            SIM_FIELD_2, SIM_FIELD_3 = convert_eb_qu(SIM_FIELD_2, SIM_FIELD_3, mapparams, eb_to_qu=1)
        else: #T, E, B: B will simply be zeroes
            pass

        SIM = np.asarray( [SIM_FIELD_1, SIM_FIELD_2, SIM_FIELD_3] )

    if bl is not None:
        if np.ndim(bl) != 2:
            bl = cl_to_cl2d(el, bl, mapparams)
        SIM = np.fft.ifft2( np.fft.fft2(SIM) * bl).real

    SIM = SIM - np.mean(SIM)  #NOTE: For the 3-field cube this removes one global mean, not each field's own; a B map of exact zeros comes back as a non-zero constant

    return SIM
