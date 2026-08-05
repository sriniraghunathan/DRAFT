"""
Helper routines supporting the ILC noise and residual calculation in get_ilc_residuals.py.

The module is grouped into six sections:

* Argument checking: ``check_freqs_in_ghz``
* Parameter file I/O: ``get_param_dict``
* Beams: ``get_bl``, ``get_beam_dic``, ``rebeam``
* Noise: ``get_nl``
* Power spectrum uncertainties: ``get_delta_cl``
* Map-domain utilities: ``get_apod_mask``, ``healpix_rotate_coords``

``check_freqs_in_ghz`` is called by ``foregrounds.py`` and ``ilc.py``.
``get_param_dict``, ``get_beam_dic`` and ``get_nl`` are currently called by ``get_ilc_residuals.py``, and ``get_bl`` is an internal helper used by the latter two.
The remainder (``rebeam``, ``get_delta_cl``, ``get_apod_mask`` and ``healpix_rotate_coords``) are standalone utilities that are currently not used elsewhere in this repository.
"""

import warnings

import numpy as np

import flatsky


# Argument checking

def check_freqs_in_ghz(*freqs):
    r"""
    Check that frequencies look like values in GHz.

    Parameters
    ----------
    *freqs : float
        Frequencies to check. ``None`` entries are skipped.

    Raises
    ------
    ValueError
        If any frequency is not positive, or is at or above :math:`10^4` and so is presumably in Hz.
    """

    bad = [nu for nu in freqs if nu is not None and not 0. < nu < 1e4]
    if bad:
        raise ValueError('Frequencies must be positive and in GHz, got %s' % (bad))


# Parameter file I/O

def get_param_dict(paramfile):
    """
    Parse a parameter file into a dictionary.

    The file is read as ``key = value`` pairs, one per line.
    Values are converted to ``bool``, ``None``, ``int`` or ``float`` where possible and left as strings otherwise.
    Parentheses are stripped from parameter names.

    Parameters
    ----------
    paramfile : str
        Path to the parameter file, e.g. ``params.ini``.

    Returns
    -------
    param_dict : dict
        Parameter names mapped to their converted values.

    Notes
    -----
    ``T``/``True`` and ``F``/``False`` become booleans, and ``None`` becomes ``None``.
    Numeric strings become ``int`` when they have no fractional part and ``float`` otherwise.
    Since each line is split on ``=``, a value containing an ``=`` character raises a ``ValueError``.
    """

    params, paramvals = np.genfromtxt(paramfile, delimiter='=', unpack=True, autostrip=True, dtype='unicode')
    param_dict = {}
    for p, pval in zip(params, paramvals):
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
            except (ValueError, OverflowError):
                pass
        # replace unallowed characters in paramname
        p = p.replace('(', '').replace(')', '')
        param_dict[p] = pval

    return param_dict


# Beams

def get_bl(beamval, el):
    r"""
    Inverse squared Gaussian beam window function.

    Returns :math:`1 / B_\ell^2`, not :math:`B_\ell` itself:

    .. math::

        \frac{1}{B_\ell^2} = \exp\!\left[ \ell (\ell+1) \frac{\theta_b^2}{8\log2} \right] .

    This is the form needed to deconvolve the beam from a noise power spectrum.
    Callers wanting :math:`B_\ell` itself should take ``get_bl(beamval, el)**-0.5``.

    Parameters
    ----------
    beamval : float
        Beam full-width at half-maximum :math:`\theta_b` in arcmin. Converted to radians internally.
    el : array_like
        Multipole moments :math:`\ell` at which to evaluate the beam.

    Returns
    -------
    bl : ndarray
        Inverse squared beam window function, same shape as ``el``.

    Notes
    -----
    :math:`B_\ell = \exp\!\left[ -\ell(\ell+1) \sigma^2/2 \right]`, with :math:`\sigma = \frac{\theta_b}{\sqrt{8 \log2}}`, is commonly used, but this function provides the inverse squared beam, :math:`1/B_\ell^2`.
    """

    fwhm_radians = np.radians(beamval/60.)
    sigma = fwhm_radians / np.sqrt(8. * np.log(2.))
    sigma2 = sigma**2
    bl = np.exp(el * (el+1) * sigma2)

    return bl


def get_beam_dic(freqs, beam_noise_dic, lmax, opbeam=None, make_2d=0, mapparams=None):
    r"""
    Beam transfer functions :math:`B_\ell` for each frequency channel.

    Parameters
    ----------
    freqs : list of int
        Frequency channels, in GHz.
    beam_noise_dic : dict
        Maps each frequency in ``freqs`` to a ``(beam_fwhm, noise)`` pair, with the beam FWHM in arcmin.
        Only the beam value is used here.
    lmax : int
        Beams are evaluated on ``np.arange(lmax)``, i.e. :math:`\ell = 0` to ``lmax - 1``.
    opbeam : float, optional
        FWHM in arcmin of a common output beam.
        If given, an additional ``'effective'`` entry is added to the returned dictionary.
        Default is ``None``.
    make_2d : int, optional
        If non-zero, convert each 1-D :math:`B_\ell` onto a 2-D flat-sky grid using :func:`flatsky.cl_to_cl2d`.
        Requires ``mapparams``.
    mapparams : list, optional
        Flat-sky map geometry ``[nx, ny, dx, dy]``, where ``dx`` and ``dy`` are the pixel sizes in arcmin.
        Required when ``make_2d`` is set. Default is ``None``.

    Returns
    -------
    bl_dic : dict
        Beam window functions :math:`B_\ell` for each entry of ``freqs``, keyed by frequency, plus an ``'effective'`` key when ``opbeam`` is given.
        Values are 1-D arrays of length ``lmax`` or 2-D arrays when ``make_2d`` is set.
    """

    bl_dic = {}
    el = np.arange(lmax)
    for freq in freqs:
        beamval, noiseval = beam_noise_dic[freq]
        ##print(beamval, noiseval)
        #get_bl() returns the inverse squared beam 1/bl^2
        bl_dic[freq] = get_bl(beamval, el)**-0.5

        if make_2d:
            if mapparams is None:
                raise ValueError('mapparams is required when make_2d is set')
            #el = np.arange(len(bl_dic[freq]))
            bl_dic[freq] = flatsky.cl_to_cl2d(el, bl_dic[freq], mapparams)

    if opbeam is not None:
        bl_dic['effective'] = get_bl(opbeam, el)**-0.5

        if make_2d:
            if mapparams is None:
                raise ValueError('mapparams is required when make_2d is set')
            bl_dic['effective'] = flatsky.cl_to_cl2d(el, bl_dic['effective'], mapparams)

    return bl_dic


def rebeam(bl_dic, threshold=1000.):
    r"""
    Ratio of the effective beam to each channel beam.

    For every frequency channel, this returns :math:`B^\mathrm{eff}_\ell / B_\ell`.

    Parameters
    ----------
    bl_dic : dict
        Beam window functions keyed by frequency, including an ``'effective'`` beam with the common output beam.
        Typically produced by :func:`get_beam_dic` with ``opbeam`` set. Only integer keys are treated as frequency channels.
    threshold : float, optional
        Upper bound on :math:`1 / B_\ell`.
        Multipoles where the channel beam is zero or negative or where the inverse would exceed this value are clamped to it. Default is 1000.

    Returns
    -------
    rebeamarr : ndarray
        Array of shape ``(n_freq, n_el)``, ordered by ascending frequency.
    """

    if 'effective' not in bl_dic:
        raise ValueError("bl_dic must contain an 'effective' key; call get_beam_dic with opbeam set")

    #freqarr = sorted( list(bl_dic.keys()) )
    freqarr = []
    for nu in list(bl_dic.keys()):
        if isinstance(nu, (int, np.integer)):
            freqarr.append(nu)
    freqarr = sorted(freqarr)

    bl_eff = bl_dic['effective']
    rebeamarr = []
    for freq in freqarr:
        #if freq == 'effective': continue
        currbl = np.copy( np.asarray(bl_dic[freq], dtype=float) )
        currinvbeamval = np.zeros_like( currbl )
        np.divide(1., currbl, out=currinvbeamval, where=(currbl > 0.))
        currinvbeamval[currbl <= 0.] = threshold
        currinvbeamval[currinvbeamval > threshold] = threshold
        rebeamval = bl_eff * currinvbeamval
        rebeamarr.append( rebeamval )

    return np.asarray( rebeamarr )


# Noise

def get_nl(
        noiseval,
        el,
        beamval,
        use_beam_window=1,
        uk_to_K=0,
        elknee=-1,
        alphaknee=0,
        beamval2=None,
        noiseval2=None,
        elknee2=-1,
        alphaknee2=0,
        rho=None,
        Nred1=-1.,
        Nred2=-1.
        #, so_like=False  #unused
        ):
    r"""
    Noise power spectrum :math:`N_\ell` for one frequency channel or the correlated noise between two channels.

    For a single channel with a :math:`1/f` (atmospheric) component:

    .. math::

        N_\ell^X = \Delta_X^2 B_\ell^{-2} \left[ 1 + \left( \frac{\ell_\mathrm{knee}}{\ell} \right)^{\!\alpha_\mathrm{knee}} \right] ,

    where :math:`\Delta_X` is the instrumental white-noise level for temperature or polarization (:math:`X = T, P`) in μK arcmin.

    When a second channel is supplied, the correlated (atmospheric) noise between the two is returned instead:

    .. math::

        N_\ell^{X,12} = \rho \, \Delta_{X, 1} \Delta_{X, 2} \left( \frac{\ell_{\mathrm{knee}, 1}}{\ell} \right)^{\!\alpha_{\mathrm{knee}, 1} / 2} \left( \frac{\ell_{\mathrm{knee}, 2}}{\ell} \right)^{\!\alpha_{\mathrm{knee}, 2} / 2} .

    The correlated noise is however zero if either channel has its atmospheric term disabled (``elknee`` or ``elknee2`` equal to ``-1``).

    Parameters
    ----------
    noiseval : float
        White-noise level of the first channel :math:`\Delta_X` (or :math:`\Delta_{X, 1}`) in μK arcmin.
    el : array_like
        Multipole moments :math:`\ell`.
    beamval : float
        Beam full-width at half-maximum of the first channel in arcmin.
    use_beam_window : int, optional
        If non-zero, deconvolve the beam by multiplying by :func:`get_bl`, i.e. by :math:`B_\ell^{-2}`. Default is 1.
    uk_to_K : int, optional
        If non-zero, convert the supplied noise levels from μK to K. Default is 0.
    elknee : float, optional
        Atmospheric-noise knee multipole :math:`\ell_\mathrm{knee}` of the first channel.
        A value of ``-1`` disables the atmospheric term. Default is -1.
    alphaknee : float, optional
        Atmospheric-noise slope :math:`\alpha_\mathrm{knee}` of the first channel. Default is 0.
    beamval2, noiseval2 : float, optional
        Beam full-width at half-maximum in arcmin and white-noise level in μK arcmin of the second channel.
        Supplying ``noiseval2`` selects the cross-channel branch; ``rho`` is then required, and ``beamval2`` is required when ``use_beam_window`` is set.
        Defaults are ``None`` and ``None``.
    elknee2, alphaknee2 : float, optional
        Atmospheric-noise knee multipole and slope of the second channel. Defaults are -1 and 0.
    rho : float, optional
        Correlation coefficient between the atmospheric noise of the two channels.
        Required when ``noiseval2`` is given. Default is ``None``.
    Nred1, Nred2 : float, optional
        Red-noise amplitudes for the two channels, internally rescaled by survey area and duration.
        A value of ``-1`` disables this treatment. Default is -1.

    Returns
    -------
    final_nl : ndarray
        Noise power spectrum, same length as ``el``.

    Notes
    -----
    The monopole is the plain white-noise level in the single-channel case to avoid dividing by zero and zero in the cross-channel case.
    """

    if Nred1 != -1:
        total_years = 5.
        fsky = 0.4 #0.35
        #survey_time = 1.
        obs_efficiency = 0.2
        noisy_map_eges_ign_factor = 0.15
        single_year = 365.25 * 24. * 3600. * obs_efficiency * (1.-noisy_map_eges_ign_factor)
        sky_area = 4. * np.pi * fsky
        year_scaling = single_year / total_years
        Nred1 = Nred1 * sky_area / year_scaling
        if Nred2 != -1.:
            Nred2 = Nred2 * sky_area / year_scaling

    cross_band_noise = 0
    if noiseval2 is not None:
        if use_beam_window and beamval2 is None:
            raise ValueError('beamval2 is required when noiseval2 is set')
        if rho is None:
            raise ValueError('rho is required when noiseval2 is set')
        cross_band_noise = 1

    if uk_to_K:
        noiseval = noiseval/1e6
        if cross_band_noise: noiseval2 = noiseval2/1e6

    el = np.asarray(el, dtype=float)
    el_safe = np.where(el > 0., el, np.inf)

    delta_T_radians = noiseval * np.radians(1./60.)
    nl = np.full(len(el), delta_T_radians**2.)
    #nl_white = np.copy(nl)

    if cross_band_noise:
        delta_T2_radians = noiseval2 * np.radians(1./60.)
    #    nl2 = np.full(len(el), delta_T2_radians**2.)
    #    nl2_white = np.copy(nl2)

    if use_beam_window:
        nl = nl * get_bl(beamval, el)
    #    if cross_band_noise:
    #        nl2 = nl2 * get_bl(beamval2, el)

    if elknee != -1.:
        if Nred1 == -1:
            nl = np.copy(nl) * (1. + (elknee * 1./el_safe)**alphaknee)
        else:
            nl = np.copy(nl) + Nred1 * (elknee * 1./el_safe)**alphaknee
    #        if cross_band_noise and elknee2 != -1.:
    #            if Nred2 == -1:
    #                nl2 = np.copy(nl2) * (1. + (elknee2 * 1./el)**alphaknee2)
    #            else:
    #                nl2 = np.copy(nl2) + Nred2*(elknee2 * 1./el)**alphaknee2

    if cross_band_noise:
        if elknee != -1. and elknee2 != -1.:
            ###final_nl = rho * nl**0.5 * nl2**0.5
            final_nl = rho * delta_T_radians * (elknee * 1./el_safe)**(alphaknee/2.) * delta_T2_radians * (elknee2 * 1./el_safe)**(alphaknee2/2.)
            #N[i,j,:] = rho * (w1*np.pi/180./60. * (ell/knee1)**(gamma1/2)) * (w2*np.pi/180./60. * (ell/knee2)**(gamma2/2))
        else:
            #no atmosphere means no correlated noise
            final_nl = np.zeros_like(nl)
    else:
        final_nl = np.copy(nl)

    return final_nl


#def get_nl_v1(noiseval, el, beamval, use_beam_window = 1, uk_to_K = 0, elknee = -1, alphaknee = 0, beamval2 = None, noiseval2 = None, elknee2 = -1, alphaknee2 = 0, rho = None, Nred1 = -1., Nred2=-1.):
#
#    cross_band_noise = 0
#    if noiseval2 is not None and beamval2 is not None:
#        assert rho is not None
#        cross_band_noise = 1
#
#    if uk_to_K:
#        noiseval = noiseval/1e6
#        if cross_band_noise: noiseval2 = noiseval2/1e6
#
#    if use_beam_window:
#        bl = get_bl(beamval, el)
#        if cross_band_noise: bl2 = get_bl(beamval2, el)
#
#    delta_T_radians = noiseval * np.radians(1./60.)
#    nl = np.tile(delta_T_radians**2., int(max(el)) + 1 )
#    nl = np.asarray( [nl[int(l)] for l in el] )
#    nl_white = np.copy(nl)
#
#    if cross_band_noise:
#        delta_T2_radians = noiseval2 * np.radians(1./60.)
#        nl2 = np.tile(delta_T2_radians**2., int(max(el)) + 1 )
#        nl2 = np.asarray( [nl2[int(l)] for l in el] )
#        nl2_white = np.copy(nl2)
#
#    if use_beam_window:
#        nl *= bl
#        if cross_band_noise: nl2 *= bl2
#
#    if elknee != -1.:
#        if Nred1==-1:
#            nl = np.copy(nl) * (1. + (elknee * 1./el)**alphaknee )
#        else:
#            nl = np.copy(nl) + Nred1*(elknee * 1./el)**alphaknee
#            if cross_band_noise and elknee2 != -1.:
#                if Nred2==-1:
#                    nl2 = np.copy(nl2) * (1. + (elknee2 * 1./el)**alphaknee2 )
#                else:
#                    nl2 = np.copy(nl2) + Nred2*(elknee2 * 1./el)**alphaknee2
#
#    if cross_band_noise and (elknee != -1. and elknee2 != -1.):
#        ###final_nl = rho * nl**0.5 * nl2**0.5
#        final_nl = rho * delta_T_radians * (elknee * 1./el)**(alphaknee/2.) * delta_T2_radians * (elknee2 * 1./el)**(alphaknee2/2.)
#        #N[i,j,:] = rho * (w1*np.pi/180./60. * (ell/knee1)**(gamma1/2)) * (w2*np.pi/180./60. * (ell/knee2)**(gamma2/2))
#    else:
#        final_nl = np.copy(nl)
#
#    return final_nl


# Power spectrum uncertainties

#def get_delta_cl(el, cl, nl, fsky=1., delta_l=1.):
def get_delta_cl(el, cl, nl=None, fsky=1., delta_l=1.):
    r"""
    Sample-variance uncertainty on a power spectrum.

    .. math::

        \Delta C_\ell = \sqrt{ \frac{2}{(2\ell + 1)\, f_\mathrm{sky}\, \Delta\ell} }\, C_\ell

    Parameters
    ----------
    el : array_like
        Multipole moments :math:`\ell`.
    cl : array_like
        Power spectrum :math:`C_\ell`.
    nl : array_like, optional
        Noise power spectrum. Accepted but *not used*, i.e. only sample variance is returned.
        Default is ``None``.
    fsky : float, optional
        Observed sky fraction. Default is 1.
    delta_l : float, optional
        Multipole bin width :math:`\Delta \ell`. Default is 1.

    Returns
    -------
    delta_cl : ndarray
        Uncertainty at each multipole :math:`\Delta C_\ell`.

    Warns
    -----
    UserWarning
        If a non-zero ``nl`` is supplied, since it does not affect the result.
    """

    if nl is not None and np.any(np.asarray(nl) != 0.):
        warnings.warn("nl provided but not used: get_delta_cl currently returns sample variance only",
                      stacklevel=2)

    delta_cl = np.sqrt(2./(2.*el + 1) / fsky / delta_l) * (cl)## + nl)

    return delta_cl


# Map-domain utilities

def get_apod_mask(ra_grid, dec_grid, mask_radius=2., taper_radius=6., in_arcmins=1):
    """
    Circular mask with a Hann taper on a flat-sky coordinate grid.

    Parameters
    ----------
    ra_grid, dec_grid : ndarray
        Two-dimensional coordinate grids, both the same shape as the output.
    mask_radius : float, optional
        Radius in arcmin of the region set to one before tapering. Default is 2.
    taper_radius : float, optional
        Length in pixels of the Hann window convolved with the mask. Default is 6.
    in_arcmins : int, optional
        If zero, the grids are taken to be in degrees and converted to arcmin. Default is 1.

    Returns
    -------
    mask : ndarray
        Tapered mask normalized to a maximum of one, same shape as ``ra_grid``.
    """

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
        inds_to_mask = np.where((radius <= mask_radius)) #2arcmins - fix this for now
        mask[inds_to_mask[0], inds_to_mask[1]] = 1.

    ker = np.hanning(taper_radius)
    ker2d = np.asarray( np.sqrt(np.outer(ker, ker)) )

    mask = ndimage.convolve(mask, ker2d)
    mask /= mask.max()

    return mask


def healpix_rotate_coords(hmap, coord):
    """
    Rotate a HEALPix map between coordinate systems.

    Parameters
    ----------
    hmap : array_like
        HEALPix map, assumed to be in RING ordering.
    coord : list of str
        Two-element ``[from, to]`` specification passed to ``healpy.Rotator``.
        Example: ``['C', 'G']`` to go from equatorial (RA/Dec) to galactic coordinates.

    Returns
    -------
    rot_hmap : ndarray
        Map of the same length as ``hmap``, with pixel values moved to their rotated positions.

    Raises
    ------
    ImportError
        If ``healpy`` is not installed. It is imported lazily so that the rest of this module can be used without it.
    """

    try:
        import healpy as H
    except ImportError:
        raise ImportError('healpix_rotate_coords requires healpy')

    #get map pixel
    pixel = np.arange(len(hmap))

    #get angles in this map first
    nside = H.get_nside(hmap)
    angles = H.pix2ang(nside, pixel)

    #rotate the angles to the desired new coordinate
    rotated_angles = H.Rotator(coord=coord)(*angles)

    #get the rotated pixel values
    rotated_pixel = H.ang2pix(nside, *rotated_angles)

    #initialize new map
    rot_hmap = np.zeros(len(pixel))

    #push the original map pixel to the new map (in the rotated pixel positions)
    rot_hmap[rotated_pixel] = hmap[pixel]

    return rot_hmap
