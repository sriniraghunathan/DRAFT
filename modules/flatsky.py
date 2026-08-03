r"""
Flat-sky routines for simulating and analyzing small patches of sky in the flat-sky approximation.

The module is grouped into four sections:

* Fourier grid: ``get_lxly``, ``get_lxly_az_angle``, ``cl_to_cl2d``
* Filters: ``get_lpf_hpf``, ``wiener_filter``
* Simulations: ``convert_eb_qu``, ``cl2map``, ``make_gaussian_realization``
* Power spectra: ``radial_profile``, ``map2cl``

``cl_to_cl2d`` is referenced by ``misc.get_beam_dic``, which calls it only when its ``make_2d`` option is set.
The remaining routines are standalone utilities that are currently not used elsewhere in this repository.

A flat-sky patch is described by ``mapparams = [nx, ny, dx, dy]``, where ``nx`` and ``ny`` are the number of pixels along each axis and ``dx`` and ``dy`` are the pixel sizes in arcminutes.
Multipoles follow the flat-sky convention :math:`\ell = \sqrt{\ell_x^2 + \ell_y^2}`, with the 2-D Fourier grid supplied by ``get_lxly``.
Power spectra carry whatever units the input ``cl`` is given in and maps carry its square root.
"""

import numpy as np


# Fourier grid

def get_lxly(mapparams):
    r"""
    Wavenumber grids :math:`\ell_x` and :math:`\ell_y` for a flat-sky patch.

    Parameters
    ----------
    mapparams : list
        Flat-sky map geometry ``[nx, ny, dx, dy]``, where ``nx`` and ``ny`` are the pixel counts along each axis and ``dx`` and ``dy`` are the pixel sizes in arcminutes.
        ``nx`` and ``ny`` must be integers.

    Returns
    -------
    lx : ndarray
        Wavenumbers varying along the second axis, of shape ``(ny, nx)``.
    ly : ndarray
        Wavenumbers varying along the first axis, of shape ``(ny, nx)``.

    Notes
    -----
    The grids follow :func:`numpy.fft.fftfreq` ordering, so the zero mode sits at index ``[0, 0]`` and the upper half of each axis carries negative wavenumbers.

    A single pixel size is used for both axes. The unpacking binds the third and fourth entries of ``mapparams`` to the same name, so the declared ``dx`` is discarded and ``dy`` is used for both. The grid is therefore correct only for square pixels, which is what the implementation assumes throughout.
    """

    nx, ny, dx, dx = mapparams  #NOTE: dy is dropped and dx ends up holding dy
    dx = np.radians(dx/60.)

    lx, ly = np.meshgrid( np.fft.fftfreq( nx, dx ), np.fft.fftfreq( ny, dx ) )
    lx *= 2*np.pi
    ly *= 2*np.pi

    return lx, ly


def get_lxly_az_angle(lx, ly):
    r"""
    Twice the azimuthal angle of each Fourier mode.

    Parameters
    ----------
    lx, ly : array_like
        Wavenumber grids, as returned by :func:`get_lxly`.

    Returns
    -------
    angle : ndarray
        :math:`2\phi_\ell`, where :math:`\phi_\ell = \arctan(\ell_x / -\ell_y)`, with the same shape as the inputs.

    Notes
    -----
    The factor of two is included because the result is used to rotate spin-2 fields in :func:`convert_eb_qu`.
    """

    return 2*np.arctan2(lx, -ly)


def cl_to_cl2d(el, cl, mapparams):
    r"""
    Interpolate a 1-D power spectrum onto the 2-D flat-sky Fourier grid.

    Each grid point is assigned :math:`C_\ell` at :math:`\ell = \sqrt{\ell_x^2 + \ell_y^2}`, with :math:`\ell_x` and :math:`\ell_y` taken from :func:`get_lxly`.

    Parameters
    ----------
    el : array_like
        Multipoles :math:`\ell` at which ``cl`` is defined. Must be strictly increasing.
    cl : array_like
        Power spectrum :math:`C_\ell`, the same length as ``el``.
    mapparams : list
        Flat-sky map geometry ``[nx, ny, dx, dy]``, where ``nx`` and ``ny`` are the pixel counts along each axis and ``dx`` and ``dy`` are the pixel sizes in arcminutes.
        ``nx`` and ``ny`` must be integers. A float pixel count, which is what a plain :class:`numpy.ndarray` of ``mapparams`` produces, raises inside :func:`numpy.fft.fftfreq`.

    Returns
    -------
    cl2d : ndarray
        Power spectrum on the 2-D grid, of shape ``(ny, nx)``.

    Raises
    ------
    ValueError
        If ``el`` is not strictly increasing since :func:`numpy.interp` would otherwise return silently incorrect values.

    Notes
    -----
    :func:`numpy.interp` clamps outside the range of ``el``, so every grid point above ``max(el)`` is assigned ``cl[-1]`` rather than zero.
    """

    if np.any(np.diff(el) <= 0):
        raise ValueError('el must be strictly increasing, but it decreases or repeats at index %s' % (int(np.argmax(np.diff(el) <= 0)) + 1))

    lx, ly = get_lxly(mapparams)
    ell = np.sqrt(lx**2. + ly**2.)

    cl2d = np.interp(ell.flatten(), el, cl).reshape(ell.shape)  #NOTE: np.interp clamps, so modes above el.max() get cl[-1] rather than zero

    return cl2d


# Filters

def get_lpf_hpf(mapparams, lmin_lmax, filter_type=0):
    """
    Top-hat low-pass, high-pass or band-pass filter on the 2-D Fourier grid.

    Parameters
    ----------
    mapparams : list
        Flat-sky map geometry ``[nx, ny, dx, dy]``.
    lmin_lmax : float or sequence of float
        Cutoff multipole. Must be a single value for ``filter_type`` 0 or 1, and an ``(lmin, lmax)`` pair for ``filter_type`` 2.
    filter_type : int, optional
        ``0`` for low pass, ``1`` for high pass and ``2`` for band pass. Default 0.

    Returns
    -------
    fft_filter : ndarray
        Filter on the 2-D grid, of shape ``(ny, nx)``, holding one inside the pass band and zero outside.

    Raises
    ------
    ValueError
        If ``filter_type`` is not 0, 1 or 2, or if ``lmin_lmax`` does not have the shape that ``filter_type`` requires.
    """

    if filter_type in [0, 1] and np.ndim(lmin_lmax) != 0:
        raise ValueError('lmin_lmax must be a single value for filter_type %s, got shape %s' % (filter_type, np.shape(lmin_lmax)))
    if filter_type == 2 and np.shape(lmin_lmax) != (2,):
        raise ValueError('lmin_lmax must be an (lmin, lmax) pair for filter_type %s, got shape %s' % (filter_type, np.shape(lmin_lmax)))

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
    r"""
    Wiener filter on the 2-D Fourier grid.

    .. math::

        W_\ell = \frac{C_\ell^\mathrm{S}}{C_\ell^\mathrm{S} + C_\ell^\mathrm{N}}\, ,

    with both spectra placed on the grid by :func:`cl_to_cl2d`.

    Parameters
    ----------
    mapparams : list
        Flat-sky map geometry ``[nx, ny, dx, dy]``.
    cl_signal : array_like
        Signal power spectrum :math:`C_\ell^\mathrm{S}`.
    cl_noise : array_like
        Noise power spectrum :math:`C_\ell^\mathrm{N}`, the same length as ``cl_signal``.
    el : array_like, optional
        Multipoles at which both spectra are defined. Defaults to ``np.arange(len(cl_signal))``.

    Returns
    -------
    wiener_filter2d : ndarray
        Filter on the 2-D grid, of shape ``(ny, nx)``, with entries between zero and one.

    Notes
    -----
    Modes at which signal and noise both vanish are set to zero, matching the treatment of the spectrum in :func:`make_gaussian_realization`.
    """

    if el is None:
        el = np.arange(len(cl_signal))

    #nx, ny, dx, dx = mapparams  #unused

    #get 2D cl
    cl_signal2d = cl_to_cl2d(el, cl_signal, mapparams)
    cl_noise2d = cl_to_cl2d(el, cl_noise, mapparams)

    wiener_filter2d = cl_signal2d / (cl_signal2d + cl_noise2d)
    wiener_filter2d[np.isnan(wiener_filter2d)] = 0.

    return wiener_filter2d


# Simulations

def convert_eb_qu(map1, map2, mapparams, eb_to_qu=1):
    r"""
    Rotate a pair of flat-sky polarization maps between E/B and Q/U.

    The transform is a rotation by :math:`2\phi_\ell` applied in Fourier space,

    .. math::

        \begin{pmatrix} Q \\ U \end{pmatrix} = \begin{pmatrix} \cos 2\phi_\ell & -\sin 2\phi_\ell \\ \sin 2\phi_\ell & \cos 2\phi_\ell \end{pmatrix} \begin{pmatrix} E \\ B \end{pmatrix}\, ,

    with the transpose of that matrix applied in the opposite direction. The rotation is orthogonal mode by mode, but see the note below on the Nyquist row and column.

    Parameters
    ----------
    map1, map2 : array_like
        The pair to rotate, E and B when ``eb_to_qu`` is set and Q and U otherwise.
    mapparams : list
        Flat-sky map geometry ``[nx, ny, dx, dy]``.
    eb_to_qu : int, optional
        If non-zero, rotate E/B to Q/U. Otherwise rotate Q/U to E/B. Default 1.

    Returns
    -------
    map1_mod, map2_mod : ndarray
        The rotated pair, Q and U when ``eb_to_qu`` is set and E and B otherwise.

    Notes
    -----
    Only the real part of the inverse transform is retained. On an even-sized grid this discards information. A pure E-mode input, and a field with no power on the Nyquist row and column, are transformed exactly.
    """

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


def cl2map(mapparams, cl, el=None):
    r"""
    Draw a Gaussian random flat-sky map with a given power spectrum.

    White noise is generated in map space and multiplied by :math:`\sqrt{C_\ell}` on the Fourier grid, with a pixel-area normalization chosen so that the recovered spectrum matches the input.

    Parameters
    ----------
    mapparams : list
        Flat-sky map geometry ``[nx, ny, dx, dy]``.
    cl : array_like
        Power spectrum :math:`C_\ell` to realize.
    el : array_like, optional
        Multipoles at which ``cl`` is defined. Defaults to ``np.arange(len(cl))``.

    Returns
    -------
    flatskymap : ndarray
        Realization with its mean removed.

    Notes
    -----
    The random field is drawn with shape ``(nx, ny)`` while the Fourier grid from :func:`get_lxly` has shape ``(ny, nx)``, so the two agree only for a square patch and a rectangular one raises ``ValueError``.

    The normalization uses a single pixel size rather than :math:`\sqrt{\mathrm{d}x\,\mathrm{d}y}`, so for non-square pixels the amplitude differs from :func:`make_gaussian_realization` by :math:`\sqrt{\mathrm{d}x/\mathrm{d}y}`.
    Use :func:`make_gaussian_realization` for correlated fields or to apply a beam.
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


def make_gaussian_realization(mapparams, el, cl, cl2=None, cl12=None, bl=None, qu_or_eb='qu'):
    r"""
    Draw one Gaussian flat-sky field or a correlated pair followed by a field of zeros.

    A single field is drawn as in :func:`cl2map`. When ``cl2`` and ``cl12`` are supplied, a second field correlated with the first is built in Fourier space from two independent white-noise realizations,

    .. math::

        \tilde{f}_2 = \tilde{g}_1 \frac{C_\ell^{12}}{\sqrt{C_\ell^{11}}} + \tilde{g}_2 \sqrt{C_\ell^{22} - \frac{(C_\ell^{12})^2}{C_\ell^{11}}}\, ,

    so that the pair carries the requested auto and cross spectra.
    A third field of zeros is appended, standing for B when the pair is T and E.

    Parameters
    ----------
    mapparams : list
        Flat-sky map geometry ``[nx, ny, dx, dy]``.
    el : array_like
        Multipoles at which the spectra are defined.
    cl : array_like
        Auto spectrum of the first field, :math:`C_\ell^{11}`.
    cl2 : array_like, optional
        Auto spectrum of the second field, :math:`C_\ell^{22}`. Supplying it selects the correlated branch, which also requires ``cl12``.
    cl12 : array_like, optional
        Cross spectrum of the two fields, :math:`C_\ell^{12}`. Required when ``cl2`` is given.
    bl : array_like, optional
        Beam transfer function :math:`B_\ell`, either 1-D over ``el`` or already on the 2-D grid, applied to the output.
    qu_or_eb : {'qu', 'eb'}, optional
        Whether to return T, Q, U or T, E, B. Only has an effect in the correlated branch. Default ``'qu'``.

    Returns
    -------
    SIM : ndarray
        A single map of shape ``(ny, nx)`` or a cube of shape ``(3, ny, nx)`` in the correlated branch, with the mean removed.

    Raises
    ------
    ValueError
        If ``qu_or_eb`` is neither ``'qu'`` nor ``'eb'`` or if ``cl2`` is given without ``cl12``.

    Notes
    -----
    The random fields are drawn with shape ``(nx, ny)`` while the Fourier grid has shape ``(ny, nx)``, so a rectangular patch raises ``ValueError``.
    """

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


# Power spectra

def radial_profile(z, xy=None, bin_size=1., minbin=0., maxbin=10., to_arcmins=1):
    """
    Average a 2-D field in annuli of constant radius.

    Parameters
    ----------
    z : array_like
        Field to profile, real or complex. Must be 2-D when ``xy`` is not given.
    xy : tuple of ndarray, optional
        Coordinate grids ``(x, y)``, each the same shape as ``z``. Defaults to :func:`numpy.indices`, i.e. pixel indices.
    bin_size : float, optional
        Width of the radial bins, in the units of ``xy``. Default 1.
    minbin, maxbin : float, optional
        Lowest and highest radius of the binning. Defaults 0 and 10.
    to_arcmins : int, optional
        If non-zero, multiply the radius by 60 before binning, i.e. convert degrees to arcminutes. Default 1.

    Returns
    -------
    radprf : ndarray
        Array of shape ``(nbin, 3)`` holding the bin center, the mean of ``z`` over the bin and an uncertainty.

    Raises
    ------
    ValueError
        If ``z`` is not 2-D and ``xy`` is not given or if the ``xy`` grids do not have the same shape as ``z``.

    Notes
    -----
    The defaults are mutually inconsistent. ``to_arcmins`` presumes ``xy`` in degrees, whereas the fallback used when ``xy`` is omitted supplies pixel indices.
    Pass ``xy`` explicitly or set ``to_arcmins`` to zero, as :func:`map2cl` does.

    The bin mean divides by the number of non-zero entries rather than by the number of entries, so a field containing exact zeros is averaged over its non-zero subset alone.
    """

    z = np.asarray(z)
    if xy is None:
        if z.ndim != 2:
            raise ValueError('z must be 2-D when xy is not given, got %s dimensions' % (z.ndim))
        x, y = np.indices(z.shape)
    else:
        x, y = xy
        if np.shape(x) != z.shape or np.shape(y) != z.shape:
            raise ValueError('xy grids must have the same shape as z, got %s and %s for z of shape %s' % (np.shape(x), np.shape(y), z.shape))

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


def map2cl(mapparams, flatskymap1, flatskymap2=None, binsize=None, minbin=100, maxbin=10000, mask=None, filter_2d=None):
    r"""
    Auto- or cross-power spectrum of one or two flat-sky maps.

    The 2-D periodogram is formed and then averaged in annuli of :math:`\ell` by :func:`radial_profile`,

    .. math::

        \hat{C}_\ell = \frac{\mathrm{d}x^2}{n_x n_y} \left\langle \tilde{m}_1 \tilde{m}_2^* \right\rangle_\ell ,

    where :math:`\mathrm{d}x` is in radians and :math:`\tilde{m}` is the unnormalized discrete transform.

    Parameters
    ----------
    mapparams : list
        Flat-sky map geometry ``[nx, ny, dx, dy]``.
    flatskymap1 : array_like
        First map, of shape ``(ny, nx)``.
    flatskymap2 : array_like, optional
        Second map, of the same shape. If given, the cross spectrum is returned in place of the auto spectrum.
    binsize : float, optional
        Width of the :math:`\ell` bins. Defaults to the :math:`\ell_x` grid spacing, i.e. along the second axis.
    minbin, maxbin : float, optional
        Lowest and highest :math:`\ell` of the binning. Defaults 100 and 10000.
    mask : array_like, optional
        Window that the caller has already applied to the maps. Only its mean is used to rescale the spectrum. The maps are not multiplied by it here.
    filter_2d : array_like, optional
        Filter already applied to the maps, of shape ``(ny, nx)``. Its radial profile is divided out of the spectrum.

    Returns
    -------
    el : ndarray
        Bin centers.
    cl : ndarray
        Binned power spectrum.

    Raises
    ------
    ValueError
        If ``flatskymap1`` does not match the grid implied by ``mapparams`` or if ``flatskymap2`` is given and does not have the same shape as ``flatskymap1``.
    """

    nx, ny, dx, dx = mapparams  #NOTE: dy is dropped and dx ends up holding dy
    dx_rad = np.radians(dx/60.)

    lx, ly = get_lxly(mapparams)
    if np.shape(flatskymap1) != lx.shape:
        raise ValueError('flatskymap1 must have the same shape as the Fourier grid built from mapparams, got %s and %s' % (np.shape(flatskymap1), lx.shape))

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
