r"""
Foreground power spectra supporting the ILC noise and residual calculation in get_ilc_residuals.py.

Precise measurements of the CMB damping tail at multipoles :math:`\ell \gtrsim 1000` are essential for constraining :math:`N_\mathrm{eff}`, but this is also the regime where astrophysical foreground emission becomes increasingly significant relative to the primary CMB signal.
Extragalactic sources dominate at small angular scales, while galactic emission is the primary contaminant at large angular scales and near the galactic plane.
The foreground contamination in polarization is comparatively smaller than in temperature at the angular scales and frequencies relevant for primary damping-tail science.
The resulting auto- and cross-frequency spectra enter the multi-frequency covariance used in the internal linear combination of :mod:`ilc`.
See also §3.1 of |paper| for additional information and details.

The module is grouped into four sections:

* Unit conversions and spectral functions: ``coth``, ``dl_to_cl``, ``get_BnuT``, ``get_dB_dT``, ``compton_y_to_delta_Tcmb``
* Power-law fitting: ``power_law``, ``perform_power_law_fit``, ``smooth_cib_spectra``
* Galactic foregrounds: ``scale_cl_dust_galactic``, ``get_cl_dust_galactic``, ``get_cl_galactic``
* Extragalactic foregrounds: ``get_foreground_power_spt``, ``get_cl_dust_cib``, ``scale_cl_dust_cib``, ``get_cl_tsz``, ``get_cl_tsz_cib``, ``get_cl_radio``

``get_BnuT``, ``get_dB_dT``, ``get_cl_galactic``, ``get_foreground_power_spt``, ``get_cl_dust_cib``, ``get_cl_tsz``, ``get_cl_tsz_cib`` and ``get_cl_radio`` are called by ``ilc.py``.
``scale_cl_dust_cib``, ``smooth_cib_spectra`` and ``get_cl_dust_galactic`` are standalone utilities that are currently not used elsewhere in this repository.
The remaining routines are internal helpers.

**Galactic foregrounds.**
The dominant galactic signals at the relevant frequencies are thermal dust emission and synchrotron radiation.
Thermal dust arises from interstellar grains heated by the interstellar radiation field, is described by a modified blackbody spectrum, rises towards higher frequencies and is brightest towards the galactic plane where the dust column density is largest.
Synchrotron radiation is emitted by cosmic-ray electrons spiralling in the galactic magnetic field and follows a power law that falls steeply towards higher frequencies, so that it dominates the lowest bands.
The subdominant free-free and anomalous-microwave-emission components are neglected.
Both components are read from map-based PySM 3 simulations rather than computed here, with the spectra evaluated on the masked sky region and the simulation files named in the parameter dictionary.
Since those simulations carry no galactic :math:`TE`, the correlation is constructed as :math:`C_\ell^{TE} = \rho_{TE} \sqrt{C_\ell^{TT} C_\ell^{EE}}` (see §3.1 of |paper|).

**Extragalactic foregrounds.**
The dominant extragalactic components are the thermal and kinematic Sunyaev-Zel'dovich effects, the cosmic infrared background and emission from radio galaxies.
The tSZ signal is sourced by inverse Compton scattering of CMB photons off the hot electrons in the intracluster medium and the kSZ effect by the Doppler shift due to the peculiar motion of galaxy clusters. Both follow a common template.
CIB emission from dusty star-forming galaxies is decomposed into a Poisson contribution from individually unresolved sources and a spatially clustered contribution, with the frequency dependence of both described by a modified blackbody spectrum.
Radio galaxy emission, which dominates at lower frequencies, is modelled as a Poisson power spectrum whose power-law frequency scaling falls towards higher frequencies, with a clustering contribution not being included.
The tSZ-CIB correlation is available through ``get_cl_tsz_cib`` (but is not included in the forecasts of |paper|).
The brightest individual point sources and galaxy clusters are assumed to be detected and masked, so the spectra modelled here are those of the unresolved sources that remain.
That masking is captured directly at the power-spectrum level through the amplitudes of the templates rather than by injecting sources into simulated maps.
Amplitudes and frequency scalings follow the parameterization based on measurements by the South Pole Telescope (see §3.1 of |paper|).

**Conventions.**
Power spectra are in units of μK² unless stated otherwise.
Frequencies are in GHz throughout.
"""

import os
import warnings

import numpy as np
from scipy.io import readsav
from scipy.optimize import curve_fit, leastsq

import misc


# Constants and module state

#Planck constant, Boltzmann constant and speed of light (in SI units)
h, k_B, c = 6.62607004e-34, 1.38064852e-23, 2.99792458e8

#data directory, resolved relative to the repo root rather than the current working directory
data_folder = os.path.join( os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'data' )

#cache for the galactic sim spectra: get_cl_galactic is called once per band pair
_cl_gal_cache = {}


# Unit conversions and spectral functions

def coth(x):
    r"""
    Hyperbolic cotangent:

    .. math::

        \coth x = \frac{e^x + e^{-x}}{e^x - e^{-x}}, .

    Parameters
    ----------
    x : array_like
        Argument.

    Returns
    -------
    coth : ndarray
        :math:`\coth x`, same shape as ``x``.
    """

    x = np.asarray(x, dtype=float)

    return (np.exp(x) + np.exp(-x)) / (np.exp(x) - np.exp(-x))


def dl_to_cl(dl, dl_fac):
    r"""
    Convert :math:`D_\ell` to :math:`C_\ell`:

    .. math::

        C_\ell = \frac{2\pi}{\ell(\ell+1)} D_\ell\, .

    Parameters
    ----------
    dl : array_like
        Power spectrum :math:`D_\ell`.
    dl_fac : array_like
        Conversion factor :math:`\ell(\ell+1)/(2\pi)`, same shape as ``dl``.

    Returns
    -------
    cl : ndarray
        Power spectrum :math:`C_\ell`, same shape as ``dl``.
        Entries where ``dl_fac`` is zero are set to zero instead of dividing, so :math:`\ell = 0` is safe.
    """

    cl = np.zeros_like( np.asarray(dl, dtype=float) )
    np.divide(dl, dl_fac, out=cl, where=(dl_fac != 0.))

    return cl


def get_BnuT(nu, temp=2.725):  #default is 2.725 (but not used), compton_y_to_delta_Tcmb uses 2.73
    r"""
    Planck blackbody spectrum:

    .. math::

        B_\nu(T) = \frac{2 h \nu^3}{c^2} \frac{1}{e^x - 1} , \qquad x = \frac{h \nu}{k_\mathrm{B} T}\, .

    Parameters
    ----------
    nu : float
        Frequency in GHz.
    temp : float, optional
        Temperature in K. Default is 2.725.

    Returns
    -------
    bnu : float
        Blackbody intensity in :math:`\mathrm{W}\,\mathrm{m}^{-2}\,\mathrm{Hz}^{-1}\,\mathrm{sr}^{-1}`.

    Raises
    ------
    ValueError
        If ``nu`` is not positive, or is at or above :math:`10^4` and so is presumably in Hz rather than GHz.
    """

    misc.check_freqs_in_ghz(nu)
    nu = nu * 1e9
    x = h*nu/(k_B*temp)

    t1 = 2 * h * nu**3. / c**2.
    t2 = 1. / (np.exp(x)-1.)

    return t1 * t2


def get_dB_dT(nu, nu0=None, temp=2.725):  #default is 2.725, compton_y_to_delta_Tcmb uses 2.73
    r"""
    Derivative of the blackbody spectrum with respect to temperature:

    .. math::

        \frac{\mathrm{d}B_\nu}{\mathrm{d}T} \propto \frac{x^4 e^x}{(e^x - 1)^2}\, , \qquad x = \frac{h \nu}{k_\mathrm{B} T}\, ,

    up to a constant that cancels in the flux to temperature ratios used throughout this module.

    Parameters
    ----------
    nu : float
        Frequency in GHz.
    nu0 : float, optional
        Reference frequency in GHz. If given, the ratio of the value at ``nu`` to the value at ``nu0`` is returned instead. Default is ``None``.
    temp : float, optional
        Temperature in K. Default is 2.725.

    Returns
    -------
    dBdT : float
        :math:`\mathrm{d}B_\nu/\mathrm{d}T` in arbitrary units, or the dimensionless ratio to its value at ``nu0``.

    Raises
    ------
    ValueError
        If ``nu`` or ``nu0`` is not positive, or is at or above :math:`10^4` and so is presumably in Hz rather than GHz.
    """

    misc.check_freqs_in_ghz(nu, nu0)
    nu = nu * 1e9

    x = h*nu/(k_B*temp)
    dBdT = x**4. * np.exp(x) / (np.exp(x)-1)**2.

    if nu0 is not None:
        nu0 = nu0 * 1e9
        x0 = h*nu0/(k_B*temp)
        dBdT0 = x0**4 * np.exp(x0) / (np.exp(x0)-1)**2.
        return dBdT / dBdT0
    else:
        return dBdT


def compton_y_to_delta_Tcmb(freq, freq_max=None, Tcmb=2.73):  #default is 2.73, get_dB_dT and get_BnuT use 2.725
    r"""
    Convert Compton-y to a CMB temperature fluctuation.

    The thermal SZ spectral function in the non-relativistic limit is

    .. math::

        g(x) = x \coth(x/2) - 4\, , \qquad x = \frac{h \nu}{k_\mathrm{B} T_\mathrm{CMB}}\, .

    This function returns :math:`T_\mathrm{CMB} \langle g \rangle`, so that :math:`\Delta T = y` times the returned value.

    Parameters
    ----------
    freq : float
        Frequency in GHz.
    freq_max : float, optional
        Upper end of a frequency band. If given, :math:`g` is averaged across the band rather than evaluated at a single frequency.  Default is ``None``. (Currently not usable.)
    Tcmb : float, optional
        CMB temperature in K. Default is 2.73.

    Returns
    -------
    delta_tcmb : float
        Conversion factor from Compton-y to :math:`\Delta T_\mathrm{CMB}`, in the units of ``Tcmb``.

    Raises
    ------
    ValueError
        If ``freq`` or ``freq_max`` is not positive, or is at or above :math:`10^4` and so is presumably in Hz rather than GHz.
    NotImplementedError
        If ``freq_max`` is given. Band averaging is not implemented.
    """

    misc.check_freqs_in_ghz(freq, freq_max)
    freq = freq * 1e9
    if freq_max is not None:
        #freq_max = freq_max * 1e9
        #freq_arr = np.arange(freq, freq_max, delta_nu)  #fix undefined delta_nu
        raise NotImplementedError('Band averaging is not implemented yet')
    else:
        freq_arr = np.asarray([freq])

    x = (h * freq_arr) / (k_B * Tcmb)
    g_nu = x * coth(x/2.) - 4.

    return Tcmb * np.mean(g_nu)


# Power-law fitting

def power_law(ell, A, alpha, ell_norm=80.):
    r"""
    Power law in multipole :math:`\ell`:

    .. math::

        A \left( \frac{\ell}{\ell_\mathrm{norm}} \right)^{\!\alpha} .

    Parameters
    ----------
    ell : array_like
        Multipole moments :math:`\ell`.
    A : float
        Amplitude at :math:`\ell = \ell_\mathrm{norm}`.
    alpha : float
        Power-law index :math:`\alpha`.
    ell_norm : float, optional
        Multipole at which the amplitude is defined. Default is 80.

    Returns
    -------
    fit : ndarray
        Power law evaluated at ``ell``, with non-finite entries set to zero. For negative ``alpha`` this includes :math:`\ell = 0`.
    """

    #print(A, alpha)#, end = ' ')
    ell = np.asarray(ell, dtype=float)
    fit = A * ((ell / ell_norm) ** alpha)
    badinds = np.where(~np.isfinite(fit))[0]
    fit[badinds] = 0.

    return fit


def perform_power_law_fit(el, cl, ell_norm=80):
    r"""
    Fit a power law to a power spectrum.

    :math:`D_\ell` is fitted with :func:`power_law` over a deliberately narrow range: the amplitude is bounded to within 5 per cent of the measured :math:`D_\ell` at ``ell_norm`` and the index to :math:`-0.3 \pm 0.1`.

    Parameters
    ----------
    el : ndarray
        Multipole moments :math:`\ell`.
    cl : ndarray
        Power spectrum :math:`C_\ell`, same shape as ``el``.
    ell_norm : float, optional
        Multipole at which the fitted amplitude is defined. Must be present in ``el``. Default is 80.

    Returns
    -------
    cl_fit : ndarray
        Fitted spectrum :math:`C_\ell` on the input multipoles, with non-finite entries set to zero.

    Raises
    ------
    ValueError
        If ``ell_norm`` is not present in ``el``, or if :math:`D_\ell` at ``ell_norm`` is not positive.
    """

    dl_fac = (el * (el + 1))/2/np.pi
    dl = cl * dl_fac
    badinds = np.where(~np.isfinite(dl))[0]
    dl[badinds] = 0.
    if not np.any(el == ell_norm):
        raise ValueError('ell_norm=%s is not present in el' % (ell_norm))
    amp_ini = dl[el == ell_norm][0]
    if amp_ini <= 0.:
        raise ValueError('Dl at ell_norm=%s is %s, but the power-law fit requires a positive amplitude' % (ell_norm, amp_ini))
    alpha_ini = -.3
    delta_alpha = 0.1
    amp_low_fac, amp_high_fac = 0.95, 1.05  #0.1, 3.

    def power_law_at_ell_norm(ell, A, alpha):
        #bind ell_norm so the fitted amplitude is Dl(ell_norm)
        return power_law(ell, A, alpha, ell_norm=ell_norm)

    pars, cov = curve_fit(f=power_law_at_ell_norm, xdata=el, ydata=dl, p0=[amp_ini, alpha_ini], bounds=((amp_ini*amp_low_fac, alpha_ini-delta_alpha), (amp_ini*amp_high_fac, alpha_ini+delta_alpha)))

    dl_fit = power_law(el, pars[0], pars[1], ell_norm=ell_norm)
    cl_fit = dl_to_cl(dl_fit, dl_fac)
    cl_fit[np.isinf(cl_fit)] = 0.
    cl_fit[np.isnan(cl_fit)] = 0.

    return cl_fit


def smooth_cib_spectra(el, cl, min_el=200, el_knee_guess=2000):
    r"""
    Smooth a CIB spectrum with a Poisson plus clustered fit.

    Multipoles above ``min_el`` are fitted by least squares with

    .. math::

        C_\ell = A \left[ 1 + \left( \frac{\ell_\mathrm{knee}}{\ell} \right)^{\!\alpha} \right]

    and the fit is interpolated back onto the input multipoles.

    Parameters
    ----------
    el : array_like
        Multipole moments :math:`\ell`.
    cl : array_like
        Power spectrum to smooth, same shape as ``el``.
    min_el : float, optional
        Multipoles at or below this value are excluded from the fit and set to zero in the output. Default is 200.
    el_knee_guess : float, optional
        Initial guess for :math:`\ell_\mathrm{knee}`. Also the lower multipole bound used to estimate the initial Poisson level. Default is 2000.

    Returns
    -------
    cl_fit : ndarray
        Smoothed spectrum on the input multipoles, zero outside the fitted range.

    Raises
    ------
    ValueError
        If no multipole exceeds ``el_knee_guess``, in which case the Poisson level cannot be estimated.

    Notes
    -----
    If the fitted knee multipole comes out negative, the initial guess is returned instead of the fit.
    """

    el_ori = np.copy(el)
    non_zero_ind = np.where(el > min_el)[0]
    el = el[non_zero_ind]
    cl = cl[non_zero_ind]

    el = el.astype(np.float64)

    def fitting_func_dust(p, x, DATA=None, return_fit=0, el_slope=1.2):
        def fitfunc(p, x):
            return p[0]*( 1 + (p[1]/x)**p[2] )
        #fitfunc = lambda p, x: p[0]*( 1 + (p[1]/x)**el_slope)
        if not return_fit:
            val = fitfunc(p, x) - DATA
            val[np.isinf(val)] = 0.
            val[np.isnan(val)] = 0.
            return val
        else:
            return fitfunc(p, x)

    #def fitfunc(x, p1, p2, p3):
    #    return p1*( 1. + (p2/x)**p3)

    if not np.any(el > el_knee_guess):
        raise ValueError('No multipoles above el_knee_guess=%s, so cannot estimate the Poisson level' % (el_knee_guess))
    poi_level = np.median(cl[el > el_knee_guess])
    el_slope = 1.2 #0.8 for Dl and 1.2 for Cl
    p0 = np.asarray( [poi_level, el_knee_guess, el_slope] )
    p1, cov, infodict, success, msg = leastsq(fitting_func_dust, p0, args=(el, cl), full_output=1)
    if p1[1] < 0:
        p1 = p0
    cl_fit = fitting_func_dust(p1, el, return_fit=1)

    cl_fit = np.interp(el_ori, el, cl_fit, left=0., right=0.)

    return cl_fit


# Galactic foregrounds

def scale_cl_dust_galactic(cl, freq1, freq2=None, freq0=278., Tdust=19.6, beta_dust=1.6):  #defaults differ from those of get_cl_galactic (20 K and 1.54)
    r"""
    Rescale a galactic dust spectrum to a pair of frequencies.

    The modified blackbody scaling of :func:`get_cl_dust_galactic` is applied with no change of multipole shape:

    .. math::

        C_\ell^{\nu_1 \nu_2} = C_\ell^{\nu_0 \nu_0} \epsilon_{\nu_1 \nu_2} \frac{\eta_{\nu_1} \eta_{\nu_2}}{\eta_{\nu_0}^2}\, , \qquad \eta_\nu = \nu^{\beta_\mathrm{d}} B_\nu(T_\mathrm{d})\, .

    Parameters
    ----------
    cl : array_like
        Dust spectrum :math:`C_\ell` at ``freq0``.
    freq1 : float
        Frequency of the first channel in GHz.
    freq2 : float, optional
        Frequency of the second channel in GHz. Default is ``freq1``.
    freq0 : float, optional
        Reference frequency in GHz. Default is 278.
    Tdust : float, optional
        Dust temperature :math:`T_\mathrm{d}` in K. Default is 19.6.
    beta_dust : float, optional
        Dust spectral index :math:`\beta_\mathrm{d}`. Default is 1.6.

    Returns
    -------
    cl_dust : ndarray
        Rescaled spectrum, same shape as ``cl``.
    """

    if freq2 is None:
        freq2 = freq1
    misc.check_freqs_in_ghz(freq0, freq1, freq2)

    nr = ( get_dB_dT(freq0) )**2.
    dr = get_dB_dT(freq1) * get_dB_dT(freq2)

    epsilon_nu1_nu2 = nr/dr

    bnu1 = get_BnuT(freq1, temp=Tdust)
    bnu2 = get_BnuT(freq2, temp=Tdust)
    bnu0 = get_BnuT(freq0, temp=Tdust)

    etanu1_dust = ((1.*freq1*1e9)**beta_dust) * bnu1
    etanu2_dust = ((1.*freq2*1e9)**beta_dust) * bnu2
    etanu0_dust = ((1.*freq0*1e9)**beta_dust) * bnu0

    cl_dust = cl * epsilon_nu1_nu2 * (1.*etanu1_dust * etanu2_dust/etanu0_dust/etanu0_dust)## * (el*1./el_norm)**el_slope

    return cl_dust


def get_cl_dust_galactic(
        el,
        freq1,
        freq2=None,
        freq0=353.,
        el_norm=80.,
        el_slope=-0.58,
        Tdust=19.6,
        Adust_freq0=4.3,
        spec_index_dust=1.6,
        return_dl=0
        ):
    r"""
    Analytic galactic dust power spectrum.

    A power law in multipole with a modified blackbody frequency scaling,

    .. math::

        D_\ell^{\nu_1 \nu_2} = A_{\nu_0} \epsilon_{\nu_1 \nu_2} \frac{\eta_{\nu_1} \eta_{\nu_2}}{\eta_{\nu_0}^2} \left( \frac{\ell}{\ell_\mathrm{norm}} \right)^{\!\mathrm{el\_slope}} .

    Parameters
    ----------
    el : array_like
        Multipole moments :math:`\ell`.
    freq1 : float
        Frequency of the first channel in GHz.
    freq2 : float, optional
        Frequency of the second channel in GHz. Default is ``freq1``.
    freq0 : float, optional
        Reference frequency in GHz. Default is 353.
    el_norm : float, optional
        Multipole at which the amplitude is defined. Default is 80.
    el_slope : float, optional
        Power-law slope of :math:`D_\ell`. Default is -0.58.
    Tdust : float, optional
        Dust temperature in K. Default is 19.6.
    Adust_freq0 : float, optional
        Amplitude :math:`A_{\nu_0}` of :math:`D_\ell` at ``el_norm`` and ``freq0``, in μK². Default is 4.3.
    spec_index_dust : float, optional
        Dust spectral index. Default is 1.6.
    return_dl : int, optional
        If non-zero, return :math:`D_\ell` instead of :math:`C_\ell`. Default is 0.

    Returns
    -------
    cl_dust : ndarray
        Dust spectrum, :math:`C_\ell` by default or :math:`D_\ell` when ``return_dl`` is set.
    """

    if freq2 is None:
        freq2 = freq1
    misc.check_freqs_in_ghz(freq0, freq1, freq2)

    nr = ( get_dB_dT(freq0) )**2.
    dr = get_dB_dT(freq1) * get_dB_dT(freq2)

    epsilon_nu1_nu2 = nr / dr

    bnu1 = get_BnuT(freq1, temp=Tdust)
    bnu2 = get_BnuT(freq2, temp=Tdust)
    bnu0 = get_BnuT(freq0, temp=Tdust)

    etanu1_dust = ((1.*freq1*1e9)**spec_index_dust) * bnu1
    etanu2_dust = ((1.*freq2*1e9)**spec_index_dust) * bnu2
    etanu0_dust = ((1.*freq0*1e9)**spec_index_dust) * bnu0

    dl_dust = Adust_freq0 * epsilon_nu1_nu2 * (1.*etanu1_dust * etanu2_dust/etanu0_dust/etanu0_dust) * (el*1./el_norm)**el_slope

    if return_dl:
        return dl_dust

    else:
        dl_fac = el * (el+1)/2/np.pi
        cl_dust = dl_to_cl(dl_dust, dl_fac)

    return cl_dust


def get_cl_galactic(
        param_dict,
        component,
        freq1,
        freq2,
        which_spec,
        which_gal_mask=None,
        bl_dic=None,
        el=None,
        use_power_law_fit=False,
        use_sed_scaling=True,
        freq0_for_sed_scaling=278,
        ell_norm=80.,
        Tdust=20.,
        beta_dust=1.54
        ):
    r"""
    Galactic foreground power spectrum from the PySM simulations.

    The spectra are read from the simulation file named in ``param_dict``, and selected by galactic mask and spectrum type.
    For dust with ``use_sed_scaling`` set, the auto-spectrum at 278 GHz is rescaled to the requested pair of frequencies with :func:`scale_cl_dust_galactic`, rather than using the simulated spectrum for that pair.

    Parameters
    ----------
    param_dict : dict
        Parameter dictionary, as returned by ``misc.get_param_dict``.
        Must contain the file name for the requested component (``cl_gal_dic_dust_fname``, ``cl_gal_dic_sync_fname`` or ``cl_gal_dic_freefree_fname``), and may contain ``cl_gal_folder`` and ``which_gal_mask``.
    component : str
        Galactic component, one of ``'dust'``, ``'sync'`` or ``'freefree'``.
    freq1, freq2 : float
        Frequencies of the two channels in GHz. Each is mapped onto the nearest simulated band internally.
    which_spec : str
        Spectrum to return, one of ``'TT'``, ``'EE'``, ``'BB'``, ``'TE'``, ``'EB'`` or ``'TB'``.
    which_gal_mask : int, optional
        Index of the galactic mask. Default is 0, but ``param_dict['which_gal_mask']`` takes precedence when present, and a conflicting value given here is ignored with a warning. Default is ``None``.
    bl_dic : dict, optional
        Beam window functions :math:`B_\ell` keyed by frequency. If given, the spectrum is beam deconvolved. Default is ``None``.
    el : array_like, optional
        Multipoles onto which the result is interpolated. If omitted, the simulation multipoles are returned. Default is ``None``.
    use_power_law_fit : bool, optional
        If ``True``, the spectrum is replaced by the power-law fit of :func:`perform_power_law_fit`. Default is ``False``.
    use_sed_scaling : bool, optional
        If ``True``, dust spectra are obtained by rescaling the 278 GHz auto-spectrum. Default is ``True``.
    freq0_for_sed_scaling : float, optional
        Reference frequency in GHz for the SED scaling. Currently ignored: the function hardcodes 278 GHz and warns if a different value is passed. Default is 278.
    ell_norm : float, optional
        Multipole at which the power-law fit is normalized. Default is 80.
    Tdust : float, optional
        Galactic dust temperature in K. Default is 20.
    beta_dust : float, optional
        Galactic dust spectral index. Default is 1.54.

    Returns
    -------
    el_gal : ndarray
        Multipole moments :math:`\ell`.
    cl_gal : ndarray
        Galactic power spectrum.

    Raises
    ------
    ValueError
        If ``component`` is not one of the supported components.
    KeyError
        If either frequency is outside the simulated bands or if the requested band pair is absent from the file in both orderings.

    Warns
    -----
    UserWarning
        If ``freq0_for_sed_scaling`` differs from 278 or if ``which_gal_mask`` conflicts with the value in ``param_dict``.

    Notes
    -----
    Synchrotron and free-free are never SED scaled, so those spectra always come from the simulations at the nearest simulated band.
    """

    misc.check_freqs_in_ghz(freq0_for_sed_scaling, freq1, freq2)
    gal_freq_dic = {20: 20, 27: 27, 39: 39, 93: 93, 90: 93, 143: 143, 145: 145, 150: 150, 225: 225, 220: 225, 278: 278, 350: 350}

    #https://healpy.readthedocs.io/en/1.5.0/generated/healpy.sphtfunc.anafast.html#healpy.sphtfunc.anafast
    spec_inds_dic = { 'TT': 0, 'EE': 1, 'BB': 2, 'TE': 3, 'EB': 4, 'TB': 5} #py2

    if component not in ['dust', 'sync', 'freefree']:
        raise ValueError("component must be 'dust', 'sync' or 'freefree', got %s" % (component))

    if 'which_gal_mask' in param_dict:
        if which_gal_mask is not None and which_gal_mask != param_dict['which_gal_mask']:
            warnings.warn('which_gal_mask=%s ignored: param_dict specifies %s' % (which_gal_mask, param_dict['which_gal_mask']), stacklevel=2)
        which_gal_mask = param_dict['which_gal_mask']
    elif which_gal_mask is None:
        which_gal_mask = 0

    if component == 'dust':
        cl_gal_dic_fname = param_dict['cl_gal_dic_dust_fname']
    elif component == 'sync':
        cl_gal_dic_fname = param_dict['cl_gal_dic_sync_fname']
    elif component == 'freefree':
        cl_gal_dic_fname = param_dict['cl_gal_dic_freefree_fname']

    try:
        cl_gal_folder = param_dict['cl_gal_folder']
        cl_gal_dic_fname = '%s/%s' % (cl_gal_folder, cl_gal_dic_fname)
    except KeyError:
        pass

    if (0):  ##component == 'sync':  #fix
        #fix me: Forcing sync. to CUmilta's simulations
        print('\n\t\tForcing sync. to CUmilta\'s simulations\n\n')
        try:
            cl_gal_dic_fname = param_dict['cl_gal_dic_sync_fname_forced']
        except KeyError:
            pass

    #pick the requested spectra: TT, EE, BB, TE, EB, TB.
    spec_ind = spec_inds_dic[which_spec]

    if cl_gal_dic_fname.find('spt_proposal_2023_13k_sqdeg_field_mask') > -1: #20230330
        gal_freq_dic = {90: 93, 95: 93, 150: 145, 220: 225}

    freq1_to_use = gal_freq_dic[freq1]
    freq2_to_use = gal_freq_dic[freq2]

    #cl_gal_dic = np.load(cl_gal_dic_fname, allow_pickle = 1, encoding = 'latin1').item()['cl_dic'][which_gal_mask]
    if cl_gal_dic_fname not in _cl_gal_cache:
        _cl_gal_cache[cl_gal_dic_fname] = np.load(cl_gal_dic_fname, allow_pickle=1, encoding='latin1').item()['cl_dic']
    cl_gal_dic = _cl_gal_cache[cl_gal_dic_fname][which_gal_mask]
    #Next line silently overrode input, so included warning.
    if freq0_for_sed_scaling != 278:
        warnings.warn('freq0_for_sed_scaling=%s ignored: get_cl_galactic hardcodes 278 GHz.' % (freq0_for_sed_scaling), stacklevel=2)
    freq0_for_sed_scaling = 278
    if (freq0_for_sed_scaling, freq0_for_sed_scaling) not in cl_gal_dic:
        use_sed_scaling = False
    if ( freq1_to_use >= max(gal_freq_dic.values()) or freq2_to_use >= max(gal_freq_dic.values()) ) and use_sed_scaling:
        cl_dust_freq0 = cl_gal_dic[ (freq0_for_sed_scaling, freq0_for_sed_scaling) ]
        if component == 'dust':
            cl_gal = scale_cl_dust_galactic(cl_dust_freq0, freq1, freq2=freq2, freq0=freq0_for_sed_scaling, Tdust=Tdust, beta_dust=beta_dust)
        else:
            cl_gal = np.zeros(cl_dust_freq0.shape)
    else:
        try:
            cl_gal = cl_gal_dic[ (freq1_to_use, freq2_to_use) ]
        except KeyError:
            cl_gal = cl_gal_dic[ (freq2_to_use, freq1_to_use) ]

        #fix me
        if np.ndim(cl_gal) == 1: #TT-only. Pol will fail.
            cl_gal = np.asarray( [cl_gal] )

    if which_spec == 'TE' and cl_gal_dic_fname.find('CUmilta') == -1:
        #force TE to be np.sqrt(TT) * np.sqrt(EE)
        cl_gal_tt, cl_gal_ee = cl_gal_dic[ (freq1, freq2) ][0], cl_gal_dic[ (freq1, freq2) ][1]

        if (1):  ##component == 'dust':  #fix
            rte = 0.35 #page 5 of https://arxiv.org/pdf/1801.04945.pdf: Discussion below Fig.5; also page 38 of https://readthedocs.org/projects/so-pysm-models/downloads/pdf/0.2.dev/
        else: ##elif component == 'sync':
            rte = 0.
        cl_gal = rte * np.sqrt( cl_gal_tt * cl_gal_ee )

    else:
        try:
            cl_gal = cl_gal[spec_ind]
        except IndexError:
            print('(%s,%s) not found for mask = %s in %s. Setting them to zeros.' % (freq1, freq2, which_spec, cl_gal_dic_fname))
            cl_gal = np.zeros( len(cl_gal[0]) )

    el_gal = np.arange( len(cl_gal) )

    #20210506: SED scaling/power-law fitting. By default we do not fit power law.
    #Note: This currently rebuilds cl_gal
    beam_deconvolved = False
    if use_sed_scaling and component == 'dust':
        cl_dust_freq0 = cl_gal_dic[ (freq0_for_sed_scaling, freq0_for_sed_scaling) ][spec_ind]
        el_gal = np.arange( len(cl_dust_freq0) )
        if bl_dic is not None:
            bl = bl_dic[freq0_for_sed_scaling]
            if len(bl) != len(cl_dust_freq0): #adjust array lengths first
                el_tmp = np.arange( len(bl) )
                cl_dust_freq0 = np.interp(el_tmp, el_gal, cl_dust_freq0, left=0., right=0.)
                el_gal = np.copy( el_tmp )
            cl_dust_freq0 = cl_dust_freq0 / bl**2.
            beam_deconvolved = True

        if use_power_law_fit:
            cl_dust_freq0 = perform_power_law_fit(el_gal, cl_dust_freq0, ell_norm=ell_norm)
        cl_gal = scale_cl_dust_galactic(cl_dust_freq0, freq1_to_use, freq2=freq2_to_use, freq0=freq0_for_sed_scaling, Tdust=Tdust, beta_dust=beta_dust)
    elif use_power_law_fit:
        if bl_dic is not None:
            bl = bl_dic[freq0_for_sed_scaling]
            if len(bl) != len(cl_gal): #adjust array lengths first
                el_tmp = np.arange( len(bl) )
                cl_gal = np.interp(el_tmp, el_gal, cl_gal, left=0., right=0.)
                el_gal = np.copy( el_tmp )
            cl_gal = cl_gal / bl**2.
            beam_deconvolved = True

        cl_gal = perform_power_law_fit(el_gal, cl_gal, ell_norm=ell_norm)

    if bl_dic is not None and beam_deconvolved is False:
        bl1 = bl_dic[freq1]
        bl2 = bl_dic[freq2]

        if len(bl1) != len(cl_gal): #adjust array lengths first
            el_tmp = np.arange( len(bl1) )
            cl_gal = np.interp(el_tmp, el_gal, cl_gal, left=0., right=0.)
            el_gal = np.copy( el_tmp )

        cl_gal = cl_gal / (bl1 * bl2)

        if (0): #20210426 - nulling highly smoothed modes  #fix
            op_beam = bl_dic[145]
            beam_ratio = op_beam**2. / (bl1 * bl2)
            highly_deconv_inds = np.where(beam_ratio >= 500)
            cl_gal[highly_deconv_inds] = 0.
            #plot(beam_ratio); ylim(1., 1000.); show(); sys.exit()

    if el is not None:
        cl_gal = np.interp(el, el_gal, cl_gal, left=0., right=0.)
        el_gal = np.copy( el )

    return el_gal, cl_gal


# Extragalactic foregrounds

def get_foreground_power_spt(component, freq1=150, freq2=None, units='uk', lmax=None):
    r"""
    Best-fit foreground power spectra from George et al. (2015).

    The spectra are read from an IDL save file produced by Christian Reichardt, converted from :math:`D_\ell` to :math:`C_\ell` and padded with zeros down to :math:`\ell = 0`.

    Parameters
    ----------
    component : str
        Foreground component, one of ``'all'``, ``'tSZ'``, ``'kSZ'``, ``'DG-Cl'``, ``'DG-Po'``, ``'RG'``, ``'tSZ-CIB'``, ``'Total'`` or ``'CMB'``.
        ``'all'`` sums the tSZ, kSZ, clustered CIB, Poisson CIB and radio terms.
    freq1 : int, optional
        Frequency band in GHz. A value of 90 is treated as 95. Default is 150.
    freq2 : int, optional
        Second frequency band in GHz, for a cross-spectrum. Default is ``freq1``.
        The pair is sorted internally, so the order of the two frequencies does not matter.
    units : str, optional
        ``'uK'`` for μK² or ``'K'`` for K². Default is ``'uK'``.
    lmax : int, optional
        If given, the output is truncated to :math:`\ell <` ``lmax``.  Default is ``None``.

    Returns
    -------
    el : ndarray
        Multipole moments :math:`\ell`, starting at zero.
    spec : ndarray
        Power spectrum :math:`C_\ell` of ``component``.

    Raises
    ------
    ValueError
        If ``component`` is not one of the listed components, or if the frequency pair is not one of the six covered by the fit: (95, 95), (95, 150), (95, 220), (150, 150), (150, 220) and (220, 220).
    """

    components = ['all', 'tSZ', 'kSZ', 'DG-Cl', 'DG-Po', 'RG', 'tSZ-CIB', 'Total', 'CMB']
    if component not in components:
        raise ValueError('%s not in list of possible foregrounds, must be one of %s' % (component, components))
    if units.lower() not in ['uk', 'k']:
        raise ValueError("units must be 'uK' or 'K', got %s" % (units))
    if lmax is not None and lmax < 0:
        raise ValueError('lmax must be non-negative, got %s' % (lmax))

    filename = '%s/george_plot_bestfit_line.sav' % (data_folder)
    data = readsav(filename)

    if freq2 is None:
        freq2 = freq1
    if freq1 == 90:
        freq1 = 95
    if freq2 == 90:
        freq2 = 95
    if freq1 > freq2:
        freq1, freq2 = freq2, freq1
    misc.check_freqs_in_ghz(freq1, freq2)

    freqs = np.asarray([(95, 95), (95, 150), (95, 220), (150, 150), (150, 220), (220, 220)])
    #dl_all = data['ml_dls'][(freqs[:, 0] == freq1) & (freqs[:, 1] == freq2)][0]
    sel = (freqs[:, 0] == freq1) & (freqs[:, 1] == freq2)
    if not sel.any():
        raise ValueError('(%s, %s) GHz is not in the George et al. (2015) results, must be one of %s' % (freq1, freq2, [tuple(pair) for pair in freqs.tolist()]))
    dl_all = data['ml_dls'][sel][0]
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

    #changing Dls to Cls
    spec /= el * (el + 1.0) / 2.0 / np.pi
    if units.lower() == 'k':
        spec /= 1e12

    #pad to l=0
    spec = np.concatenate((np.zeros(min(el)), spec))
    el = np.concatenate((np.arange(min(el)), el))

    if lmax is not None:
        el = el[:lmax]
        spec = spec[:lmax]

    return el, spec


def get_cl_dust_cib(freq1, freq2, fg_model='george15', freq0=150, spec_index_dg_po=1.505-0.077, spec_index_dg_clus=2.51-0.2, Tcib=20., reduce_cib_power=None):
    r"""
    Cosmic infrared background power spectra (Poisson and clustered).

    The George et al. 2015 spectra at ``freq0`` are anchored at :math:`\ell = 3000` and scaled to the requested pair of frequencies with a modified blackbody,

    .. math::

        \eta_\nu = \nu^\beta B_\nu(T_\mathrm{CIB})\, ,

    and converted from flux to CMB temperature units with :math:`\epsilon_{\nu_1 \nu_2} = (\mathrm{d}B/\mathrm{d}T|_{\nu_0})^2 / (\mathrm{d}B/\mathrm{d}T|_{\nu_1} \mathrm{d}B/\mathrm{d}T|_{\nu_2})`.
    The Poisson term then follows :math:`D_\ell \propto \ell^2` and the clustered term :math:`D_\ell \propto \ell^{0.8}`.

    Parameters
    ----------
    freq1, freq2 : float
        Frequencies of the two channels in GHz.
    fg_model : str, optional
        Foreground model, ``'george15'`` or ``'reichardt21'``. Both use the same template here. Default is ``'george15'``.
    freq0 : float, optional
        Reference frequency in GHz at which the template is normalized. Default is 150.
    spec_index_dg_po : float, optional
        Modified blackbody spectral index :math:`\beta` of the Poisson term. Default is ``1.505 - 0.077``.
    spec_index_dg_clus : float, optional
        Modified blackbody spectral index :math:`\beta` of the clustered term. Default is ``2.51 - 0.2``.
    Tcib : float, optional
        CIB dust temperature in K. Default is 20.
    reduce_cib_power : float, optional
        If given, the ``freq0`` CIB power is divided by this factor before scaling. Default is ``None``.

    Returns
    -------
    el : ndarray
        Multipole moments :math:`\ell`.
    cl_dg_po : ndarray
        Poisson CIB spectrum.
    cl_dg_clus : ndarray
        Clustered CIB spectrum.

    Raises
    ------
    ValueError
        If ``fg_model`` is not one of the supported models.
    """

    misc.check_freqs_in_ghz(freq0, freq1, freq2)
    if fg_model not in ['george15', 'reichardt21']:
        raise ValueError("fg_model must be 'george15' or 'reichardt21', got %s" % (fg_model))
    el, cl_dg_po_freq0 = get_foreground_power_spt('DG-Po', freq1=freq0, freq2=freq0)
    el, cl_dg_clus_freq0 = get_foreground_power_spt('DG-Cl', freq1=freq0, freq2=freq0)
    el_norm = 3000

    #convert to Dls
    dl_fac = el * (el+1)/2/np.pi
    dl_dg_po = dl_fac * cl_dg_po_freq0
    dl_dg_clus = dl_fac * cl_dg_clus_freq0

    if reduce_cib_power is not None:  #reduce 150 GHz CIB power: useful for CMB-HD
        dl_dg_po = dl_dg_po/reduce_cib_power
        dl_dg_clus = dl_dg_clus/reduce_cib_power

    nr = ( get_dB_dT(freq0) )**2.
    dr = get_dB_dT(freq1) * get_dB_dT(freq2)

    epsilon_nu1_nu2 = nr/dr

    bnu1 = get_BnuT(freq1, temp=Tcib)
    bnu2 = get_BnuT(freq2, temp=Tcib)
    bnu0 = get_BnuT(freq0, temp=Tcib)

    etanu1_dg_po = ((1.*freq1*1e9)**spec_index_dg_po) * bnu1
    etanu2_dg_po = ((1.*freq2*1e9)**spec_index_dg_po) * bnu2
    etanu0_dg_po = ((1.*freq0*1e9)**spec_index_dg_po) * bnu0

    etanu1_dg_clus = ((1.*freq1*1e9)**spec_index_dg_clus) * bnu1
    etanu2_dg_clus = ((1.*freq2*1e9)**spec_index_dg_clus) * bnu2
    etanu0_dg_clus = ((1.*freq0*1e9)**spec_index_dg_clus) * bnu0

    dl_dg_po = dl_dg_po[el == el_norm][0] * epsilon_nu1_nu2 * (1.*etanu1_dg_po * etanu2_dg_po/etanu0_dg_po/etanu0_dg_po) * (el*1./el_norm)**2
    dl_dg_clus = dl_dg_clus[el == el_norm][0] * epsilon_nu1_nu2 * (1.*etanu1_dg_clus * etanu2_dg_clus/etanu0_dg_clus/etanu0_dg_clus) * (el*1./el_norm)**0.8

    cl_dg_po = dl_to_cl(dl_dg_po, dl_fac)
    cl_dg_clus = dl_to_cl(dl_dg_clus, dl_fac)

    cl_dg_po[np.isnan(cl_dg_po)] = 0.
    cl_dg_po[np.isinf(cl_dg_po)] = 0.
    cl_dg_clus[np.isnan(cl_dg_clus)] = 0.
    cl_dg_clus[np.isinf(cl_dg_clus)] = 0.

    return el, cl_dg_po, cl_dg_clus


def scale_cl_dust_cib(el, cl_dust_freq0, freq0, freq1, freq2, beta, Tcib, el_slope, el_norm=3000):
    r"""
    Rescale a dust power spectrum from one frequency to a pair of frequencies.

    The input spectrum at ``freq0`` is anchored at ``el_norm``, scaled with the modified blackbody of :func:`get_cl_dust_cib` and given the shape :math:`D_\ell \propto \ell^{\mathrm{el\_slope}}`.

    Parameters
    ----------
    el : array_like
        Multipole moments :math:`\ell`.
    cl_dust_freq0 : array_like
        Dust spectrum :math:`C_\ell` at ``freq0``, same shape as ``el``.
    freq0 : float
        Reference frequency in GHz.
    freq1, freq2 : float
        Frequencies of the two channels in GHz.
    beta : float
        Modified blackbody spectral index.
    Tcib : float
        Dust temperature in K.
    el_slope : float
        Power-law slope of :math:`D_\ell`.
    el_norm : float, optional
        Multipole at which the input spectrum is anchored. Must be present in ``el``. Default is 3000.

    Returns
    -------
    el : ndarray
        Multipole moments, unchanged from the input.
    cl_dust : ndarray
        Rescaled spectrum :math:`C_\ell`.

    Raises
    ------
    IndexError
        If ``el_norm`` is not present in ``el``.
    """

    misc.check_freqs_in_ghz(freq0, freq1, freq2)

    #convert to Dls
    dl_fac = el * (el+1)/2/np.pi
    dl_dust_freq0 = dl_fac * cl_dust_freq0

    nr = ( get_dB_dT(freq0) )**2.
    dr = get_dB_dT(freq1) * get_dB_dT(freq2)

    epsilon_nu1_nu2 = nr/dr

    bnu1 = get_BnuT(freq1, temp=Tcib)
    bnu2 = get_BnuT(freq2, temp=Tcib)
    bnu0 = get_BnuT(freq0, temp=Tcib)

    etanu1 = ((1.*freq1*1e9)**beta) * bnu1
    etanu2 = ((1.*freq2*1e9)**beta) * bnu2
    etanu0 = ((1.*freq0*1e9)**beta) * bnu0

    dl_dust = dl_dust_freq0[el == el_norm][0] * epsilon_nu1_nu2 * (1.*etanu1 * etanu2/etanu0/etanu0) * (el*1./el_norm)**el_slope

    cl_dust = dl_to_cl(dl_dust, dl_fac)
    cl_dust[np.isnan(cl_dust)] = 0.

    return el, cl_dust


def get_cl_tsz(freq1, freq2, freq0=150, fg_model='george15', reduce_tsz_power=None):
    r"""
    Thermal Sunyaev-Zel'dovich power spectrum.

    The template at ``freq0`` is scaled by the ratio of tSZ spectral functions,

    .. math::

        C_\ell^{\nu_1 \nu_2} = C_\ell^{\nu_0 \nu_0} \frac{g(\nu_1) g(\nu_2)}{g(\nu_0)^2}\, ,

    with :math:`g` as returned by :func:`compton_y_to_delta_Tcmb`.

    Parameters
    ----------
    freq1, freq2 : float
        Frequencies of the two channels in GHz.
    freq0 : float, optional
        Reference frequency in GHz. Default is 150.
    fg_model : str, optional
        Foreground model, ``'george15'`` or ``'reichardt21'``. Both use the same template here. Default is ``'george15'``.
    reduce_tsz_power : float, optional
        If given, the scaled spectrum is divided by this factor. Default is ``None``.

    Returns
    -------
    el : ndarray
        Multipole moments :math:`\ell`.
    cl_tsz : ndarray
        tSZ power spectrum.

    Raises
    ------
    ValueError
        If ``fg_model`` is not one of the supported models.
    """

    misc.check_freqs_in_ghz(freq0, freq1, freq2)
    if fg_model not in ['george15', 'reichardt21']:
        raise ValueError("fg_model must be 'george15' or 'reichardt21', got %s" % (fg_model))
    el, cl_tsz_freq0 = get_foreground_power_spt('tSZ', freq1=freq0, freq2=freq0)

    tsz_fac_freq0 = compton_y_to_delta_Tcmb(freq0)
    tsz_fac_freq1 = compton_y_to_delta_Tcmb(freq1)
    tsz_fac_freq2 = compton_y_to_delta_Tcmb(freq2)

    scalefac = tsz_fac_freq1 * tsz_fac_freq2 / (tsz_fac_freq0**2.)

    cl_tsz = cl_tsz_freq0 * scalefac
    cl_tsz[np.isnan(cl_tsz)] = 0.
    cl_tsz[np.isinf(cl_tsz)] = 0.

    if reduce_tsz_power is not None:
        cl_tsz /= reduce_tsz_power

    return el, cl_tsz


def get_cl_tsz_cib(freq1, freq2, freq0=150, fg_model='george15', spec_index_dg_po=1.505-0.077, spec_index_dg_clus=2.51-0.2, Tcib=20., cl_cib_dic=None, reduce_tsz_power=None):  #, cib_flux_threshold=1.5):
    r"""
    Correlated thermal SZ and CIB power spectrum:

    .. math::

        C_\ell^{\mathrm{tSZ} \times \mathrm{CIB}} = -\rho \left( \sqrt{C_\ell^{\mathrm{tSZ}, \nu_1 \nu_1} C_\ell^{\mathrm{CIB}, \nu_2 \nu_2}} + \sqrt{C_\ell^{\mathrm{tSZ}, \nu_2 \nu_2} C_\ell^{\mathrm{CIB}, \nu_1 \nu_1}} \right) .

    The correlation coefficient :math:`\rho` is 0.1 for ``'george15'`` and 0.078 for ``'reichardt21'``.

    Parameters
    ----------
    freq1, freq2 : float
        Frequencies of the two channels in GHz.
    freq0 : float, optional
        Reference frequency in GHz. Default is 150.
    fg_model : str, optional
        Foreground model, ``'george15'`` or ``'reichardt21'``. Only :math:`\rho` differs between them. Default is ``'george15'``.
    spec_index_dg_po, spec_index_dg_clus : float, optional
        Modified blackbody spectral indices of the Poisson and clustered CIB terms. Default is ``1.505 - 0.077`` and  ``2.51 - 0.2``.
    Tcib : float, optional
        CIB dust temperature in K. Default is 20.
    cl_cib_dic : dict, optional
        CIB auto-spectra keyed by ``(freq, freq)``, each an ``(el, cl)`` pair.
        If given, these replace the spectra that would otherwise come from :func:`get_cl_dust_cib`.
        Default is ``None``.
    reduce_tsz_power : float, optional
        If given, both tSZ auto-spectra are divided by this factor. Default is ``None``.

    Returns
    -------
    el : ndarray
        Multipole moments :math:`\ell`.
    cl_tsz_cib : ndarray
        Correlated tSZ x CIB spectrum.

    Raises
    ------
    ValueError
        If ``fg_model`` is not one of the supported models.
    """

    misc.check_freqs_in_ghz(freq0, freq1, freq2)
    if fg_model not in ['george15', 'reichardt21']:
        raise ValueError("fg_model must be 'george15' or 'reichardt21', got %s" % (fg_model))
    if fg_model == 'george15':
        corr_coeff = 0.1
    elif fg_model == 'reichardt21':
        corr_coeff = 0.078

    el, cl_tsz_freq1_freq1 = get_cl_tsz(freq1, freq1, freq0=freq0, fg_model=fg_model, reduce_tsz_power=reduce_tsz_power)
    if cl_cib_dic is not None:
        el, cl_dg_freq1_freq1 = cl_cib_dic[(freq1, freq1)]
    else:
        #get tSZ and CIB spectra for freq1
        el, cl_dg_po_freq1_freq1, cl_dg_clus_freq1_freq1 = get_cl_dust_cib(freq1, freq1, freq0=freq0, fg_model=fg_model, spec_index_dg_po=spec_index_dg_po, spec_index_dg_clus=spec_index_dg_clus, Tcib=Tcib)
        cl_dg_freq1_freq1 = cl_dg_po_freq1_freq1 + cl_dg_clus_freq1_freq1

    #get tSZ and CIB spectra for freq2
    el, cl_tsz_freq2_freq2 = get_cl_tsz(freq2, freq2, freq0=freq0, fg_model=fg_model, reduce_tsz_power=reduce_tsz_power)
    if cl_cib_dic is not None:
        el, cl_dg_freq2_freq2 = cl_cib_dic[(freq2, freq2)]
    else:
        el, cl_dg_po_freq2_freq2, cl_dg_clus_freq2_freq2 = get_cl_dust_cib(freq2, freq2, freq0=freq0, fg_model=fg_model, spec_index_dg_po=spec_index_dg_po, spec_index_dg_clus=spec_index_dg_clus, Tcib=Tcib)
        cl_dg_freq2_freq2 = cl_dg_po_freq2_freq2 + cl_dg_clus_freq2_freq2
    if len(cl_dg_freq1_freq1) != len(cl_dg_freq2_freq2):
        raise ValueError('The CIB spectra for %s and %s have different lengths (%s and %s)'
                         % (freq1, freq2, len(cl_dg_freq1_freq1), len(cl_dg_freq2_freq2)))
    if len(el) != len(cl_tsz_freq2_freq2):
        cl_tsz_freq1_freq1 = np.interp(el, np.arange(len(cl_tsz_freq1_freq1)), cl_tsz_freq1_freq1)
        cl_tsz_freq2_freq2 = np.interp(el, np.arange(len(cl_tsz_freq2_freq2)), cl_tsz_freq2_freq2)

    #20220325
    if freq1 >= 217 and freq2 >= 217:
        corr_coeff = corr_coeff * -1.

    cl_tsz_cib = -corr_coeff * ( np.sqrt(cl_tsz_freq1_freq1 * cl_dg_freq2_freq2) + np.sqrt(cl_tsz_freq2_freq2 * cl_dg_freq1_freq1) )

    return el, cl_tsz_cib


def get_cl_radio(freq1, freq2, freq0=150, fg_model='george15', spec_index_rg=-0.9, null_highfreq_radio=1, reduce_radio_power_150=None):
    r"""
    Radio galaxy power spectrum.

    The template at ``freq0`` is anchored at :math:`\ell = 3000` and scaled as a power law in frequency,

    .. math::

        D_\ell^{\nu_1 \nu_2} = D_{3000}^{\nu_0 \nu_0} \epsilon_{\nu_1 \nu_2} \left( \frac{\nu_1 \nu_2}{\nu_0^2} \right)^{\!\alpha_\mathrm{rg}} \left( \frac{\ell}{3000} \right)^{\!2} ,

    with :math:`\epsilon_{\nu_1 \nu_2}` the flux to CMB temperature conversion defined in :func:`get_cl_dust_cib`.

    Parameters
    ----------
    freq1, freq2 : float
        Frequencies of the two channels in GHz.
    freq0 : float, optional
        Reference frequency in GHz. Default is 150.
    fg_model : str, optional
        Foreground model. Only ``'george15'`` is implemented here. Default is ``'george15'``.
    spec_index_rg : float, optional
        Radio spectral index :math:`\alpha_\mathrm{rg}`. Default is -0.9.
    null_highfreq_radio : int, optional
        If non-zero, the spectrum is set to zero when either frequency exceeds 230 GHz. Default is 1.
    reduce_radio_power_150 : float, optional
        If given, the template is divided by this factor before scaling. Default is ``None``.

    Returns
    -------
    el : ndarray
        Multipole moments :math:`\ell`.
    cl_rg : ndarray
        Radio galaxy power spectrum.

    Raises
    ------
    ValueError
        If ``fg_model`` is not ``'george15'``.
    """

    misc.check_freqs_in_ghz(freq0, freq1, freq2)
    if fg_model != 'george15':
        raise ValueError("get_cl_radio only implements fg_model='george15', got %s" % (fg_model))

    if fg_model == 'george15':
        el, cl_rg_freq0 = get_foreground_power_spt('RG', freq1=freq0, freq2=freq0)
        if reduce_radio_power_150 is not None:
            cl_rg_freq0 /= reduce_radio_power_150
        el_norm = 3000

    #convert to Dls
    dl_fac = el * (el+1)/2/np.pi
    dl_rg = dl_fac * cl_rg_freq0

    nr = ( get_dB_dT(freq0) )**2.
    dr = get_dB_dT(freq1) * get_dB_dT(freq2)

    epsilon_nu1_nu2 = nr/dr

    dl_rg = dl_rg[el == el_norm][0] * epsilon_nu1_nu2 * (1.*freq1 * freq2/freq0/freq0)**spec_index_rg * (el*1./el_norm)**2

    cl_rg = dl_to_cl(dl_rg, dl_fac)

    cl_rg[np.isnan(cl_rg)] = 0.

    if null_highfreq_radio and (freq1 > 230 or freq2 > 230):
        cl_rg *= 0.

    return el, cl_rg
