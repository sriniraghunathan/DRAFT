r"""
Cosmological Fisher matrices from delensed CMB spectra, driven through FisherLens.

This module supplies the forecasting stage of DRAFT based on Fisher information theory.
It takes the spectra, parameter derivatives, effective noise and lensing-reconstruction noise produced by :mod:`delensing` and forms the Fisher matrix :math:`F_{\alpha\beta}` for one survey configuration, together with the separate large-scale Fisher matrix contributed by a satellite.
Combining those matrices into a total, applying priors, fixing parameters and reducing to the projected uncertainties :math:`\sigma(p_\alpha)` is the business of :func:`combine_fisher`, which is the only route to a :math:`\sigma` and is built on top of the separate matrices assembled here.
Mismodeled foreground residuals displace those parameters rather than inflating their uncertainties and that displacement is computed by :func:`parameter_bias`, whose implementation is still experimental.

The module is grouped into the following sections:

* Multipole windows: :func:`fisher_windows`
* Fisher blocks: :func:`survey_fisher_block`, :func:`satellite_fisher_block`, :func:`foreground_systematic`, :func:`bias_vector`
* Totals: :func:`prior_matrix`, :func:`invert_fisher`, :func:`combine_fisher`, :func:`parameter_bias`

The functions here are named for Fisher blocks rather than Fisher matrices to signal what they do not return.
A block is the Fisher matrix over every varied parameter for one measurement, one survey configuration or the low-ell information from a satellite, carried in a dictionary alongside the information identifying how it was formed.
It is a summand of the total, never the total itself and never a submatrix of a larger Fisher matrix.
Every block goes through :func:`combine_fisher`, including when the total consists of a single block for a single survey since that is what ensures no double-counting, e.g. of the low-ell information or the priors.

The two kinds of block differ in more than the multipoles they cover, which is why they are separate functions rather than one with a switch.
A survey block carries one Fisher matrix per spectrum type since delensing is a property of the high-ell (usually ground-based) measurement and the forecast is quoted from the delensed spectra with the lensed and unlensed ones alongside for comparison.
A satellite block carries a single matrix, formed from the lensed spectra, because delensing is essentially irrelevant for the large scales covered by a satellite; :func:`combine_fisher` therefore adds that one matrix to the total of every spectrum type.
A survey block also sums over temperature, polarization and the lensing deflection, while a satellite block inverts a one-by-one temperature covariance at each multipole, which is the likelihood of the temperature-only measurement a satellite contributes on those scales (since we add the polarization information via a tau prior).
Each survey configuration contributes its own block, while the satellite contributes one once.

FisherLens contains no sky fraction anywhere, so each block is scaled by its own :math:`f_\mathrm{sky}` here.
Power spectra are in units of μK² throughout, matching :mod:`delensing` and the ``cl_residual`` entries of the ILC product files.
Multipoles are held under ``el`` in every dictionary this module builds, as elsewhere in DRAFT, while the dictionaries handed to FisherLens keep the ``l`` that it requires.
"""

# Note:
# Several routines in this module work around behavior of FisherLens that was established empirically rather than from its documentation.
# Each such workaround carries a comment naming the pinned commit it was verified against since the reasoning is not otherwise recoverable from the code.
# The pinned commits are those recorded in the respective ``.gitmodules``; if either submodule pointer is advanced, those comments are the places to re-check.
# We do not directly employ ``paperPlots/plotTools.py`` of FisherLens, which is an analysis script that suits its own notebooks and is more specific than suitable here.

import copy
import warnings

import numpy as np

import delensing


# Constants

POL_COMBS = ['cl_TT', 'cl_TE', 'cl_EE', 'cl_dd']
"""
Spectra entering the survey Fisher sum, in the order FisherLens indexes the covariance by.

The order is significant: ``fisherTools.getGaussianCov`` locates each spectrum by its position in this list, so reordering it reorders the covariance without changing anything else.
"""

SPECTRUM_TYPES = ['unlensed', 'lensed', 'delensed']
"""
Spectrum types for which a Fisher matrix is formed by default.
"""

SATELLITE_POL_COMBS = ['cl_TT']
"""
Spectra entering the satellite Fisher block.

Temperature alone since the satellite contributes large-scale temperature information and its polarization enters the forecast as a prior on the optical depth rather than as a spectrum.
"""

SATELLITE_SPECTRUM = 'lensed'
"""
Spectrum type used for the satellite block.

Lensed rather than delensed since delensing is a property of the ground-based measurement and does not apply to the satellite's own large-scale temperature data.
"""

DEFAULT_WINDOWS = {'lmin': 30, 'lmax': 5000, 'lmax_TT': 5000, 'lmin_dd': 2, 'lmax_dd': 5000}
"""
Default multipole ranges of the survey Fisher sum.
"""

#The ILC residual is decomposed into one entry per foreground component and one for the instrumental noise. Mismodeling the noise is a real systematic, but it is not a foreground bias and letting it in here would put it under a heading where nobody would look for it.
NOT_FOREGROUNDS = ('noise',)
"""
Entries of ``fg_res_dic`` that are not foregrounds and may not be given a mismodeling fraction.
"""


# Multipole windows

def _check_windows(windows):
    """
    Validate a set of multipole windows.

    Parameters
    ----------
    windows : dict
        Windows keyed by spectrum, each a two-element sequence ``[lmin, lmax]``.

    Raises
    ------
    ValueError
        If any window is malformed, has a lower bound below 2 or does not have its lower bound below its upper bound.
    """

    for spec in sorted(windows):
        if spec == 'lmaxCov':
            continue
        window = windows[spec]
        if len(window) != 2:
            raise ValueError('The window for %s must have two entries, got %s' % (spec, window))
        low, high = int(window[0]), int(window[1])
        #A lower bound below 2 is rejected rather than clamped, because createEllsToUseDict from
        #FisherLens (652eaec) clamps it silently and mutates the caller's list while doing so.
        if low < 2:
            raise ValueError('The lower bound of the %s window must be at least 2, got %d. CLASS '
                             'returns spectra from multipole 2 and FisherLens would clamp this '
                             'silently.' % (spec, low))
        if low >= high:
            raise ValueError('The %s window must have its lower bound below its upper bound, got '
                             '[%d, %d]' % (spec, low, high))


def fisher_windows(settings=None, lmin=None, lmax=None, lmax_TT=None, lmin_dd=None, lmax_dd=None):
    r"""
    Assemble the multipole windows of the survey Fisher sum.

    Parameters
    ----------
    settings : dict, optional
        Forecasting settings, as returned by :func:`get_fisher_forecasts.load_fisher_params`, from which ``lmin``, ``lmax``, ``lmax_TT``, ``lmin_dd`` and ``lmax_dd`` are read. Default is ``None``, which uses :data:`DEFAULT_WINDOWS`.
    lmin : int, optional
        Lowest multipole of the :math:`TT`, :math:`TE` and :math:`EE` spectra. Default is ``None``, which takes the value from ``settings``.
    lmax : int, optional
        Highest multipole of the :math:`TE` and :math:`EE` spectra. Default is ``None``, which takes the value from ``settings``.
    lmax_TT : int, optional
        Highest multipole of the :math:`TT` spectrum, kept separate so the temperature can be truncated below the polarization. Default is ``None``, which takes the value from ``settings``.
    lmin_dd : int, optional
        Lowest multipole of the lensing-deflection spectrum. Default is ``None``, which takes the value from ``settings``.
    lmax_dd : int, optional
        Highest multipole of the lensing-deflection spectrum. Default is ``None``, which takes the value from ``settings``.

    Returns
    -------
    dict
        A window per spectrum of :data:`POL_COMBS`, each a two-element list ``[lmin, lmax]``, inclusive at both ends.

    Raises
    ------
    ValueError
        If any lower bound is below 2 or does not lie below its matching upper bound.
    """

    #Note:
    #A fresh dictionary with fresh lists is returned on every call, because ``fisherTools.createEllsToUseDict`` edits the windows it is given in place and returns the same object rather than a copy, so a dictionary reused across calls would be silently modified. This was verified against FisherLens ``652eaec``: a window given as ``[0, 29]`` leaves the caller's list reading ``[2, 29]``.
    #The windows are inclusive at both ends. ``getGaussianCMBFisher`` applies them as the slice ``[lmin-2 : lmax-1]`` of an array whose element zero is :math:`\ell = 2`, which retains both bounds.
    #No ``lmaxCov`` entry is set here. It is read only on the non-Gaussian path, where :func:`survey_fisher_block` supplies it, and setting it on the Gaussian path would have no effect since a covariance diagonal in multipole contributes nothing outside the windows.

    if settings is None:
        settings = DEFAULT_WINDOWS
    values = {}
    for name, given in (('lmin', lmin), ('lmax', lmax), ('lmax_TT', lmax_TT),
                        ('lmin_dd', lmin_dd), ('lmax_dd', lmax_dd)):
        if given is not None:
            values[name] = int(given)
        elif name in settings:
            values[name] = int(settings[name])
        elif name in DEFAULT_WINDOWS:
            values[name] = DEFAULT_WINDOWS[name]
        else:
            raise ValueError('No value for %s and none in the settings supplied. Its keys are %s'
                             % (name, sorted(settings.keys())))

    windows = {'cl_TT': [values['lmin'], values['lmax_TT']],
               'cl_TE': [values['lmin'], values['lmax']],
               'cl_EE': [values['lmin'], values['lmax']],
               'cl_dd': [values['lmin_dd'], values['lmax_dd']]}

    _check_windows(windows)

    return windows


# Fisher blocks

def _check_powers(powers, spectrum_types, pol_combs):
    """
    Validate the spectra handed to a block.

    Parameters
    ----------
    powers : dict
        Spectra keyed by spectrum type.
    spectrum_types : sequence of str
        Spectrum types that will be used.
    pol_combs : sequence of str
        Spectra that will be used.

    Raises
    ------
    ValueError
        If ``'unlensed'`` is absent, if a requested spectrum type is absent, if ``'lensing'`` is absent while the deflection spectrum is used or if a requested spectrum is absent from a spectrum type that will be used.
    """

    #getGaussianCov reads powersFid['unlensed']['l'] whatever spectrum types are requested, so the
    #entry must be present even when it is not among them (FisherLens 652eaec).
    if 'unlensed' not in powers:
        raise ValueError("The spectra must include an 'unlensed' entry whatever spectrum types are "
                         'requested since FisherLens reads its multipole array to set up the '
                         'covariance. The entries present are %s' % (sorted(powers.keys())))
    if 'l' not in powers['unlensed']:
        raise ValueError("The 'unlensed' spectra must hold their multipoles under 'l', as CLASS "
                         'returns them. Its entries are %s' % (sorted(powers['unlensed'].keys())))

    for spectrum in spectrum_types:
        if spectrum not in powers:
            raise ValueError('There are no %s spectra. The entries present are %s'
                             % (spectrum, sorted(powers.keys())))
    if 'cl_dd' in pol_combs and 'lensing' not in powers:
        raise ValueError("The deflection spectrum is among the spectra to sum over, so the spectra "
                         "must include a 'lensing' entry. The entries present are %s"
                         % (sorted(powers.keys())))

    for spectrum in spectrum_types:
        for spec in pol_combs:
            source = 'lensing' if spec == 'cl_dd' else spectrum
            if spec not in powers[source]:
                raise ValueError('The %s spectra have no %s entry, which is needed to build the '
                                 'covariance. Its entries are %s'
                                 % (source, spec, sorted(powers[source].keys())))


def _check_derivs(param_derivs, vary_params, spectrum_types, pol_combs):
    """
    Validate the parameter derivatives handed to a block.

    Parameters
    ----------
    param_derivs : dict
        Derivatives keyed by parameter, then spectrum type, then spectrum.
    vary_params : sequence of str
        Parameters that will be used.
    spectrum_types : sequence of str
        Spectrum types that will be used.
    pol_combs : sequence of str
        Spectra that will be used.

    Raises
    ------
    ValueError
        If a requested parameter is absent, if a requested spectrum type is absent for a parameter or if a requested spectrum is absent, including the deflection derivative under ``'lensing'``.
    """

    for name in vary_params:
        if name not in param_derivs:
            raise ValueError('There is no derivative for %s. The parameters differentiated are %s'
                             % (name, sorted(param_derivs.keys())))
        #getPowerDerivWithParams reports a missing 'lensing' entry by printing and carrying on, so
        #the gap surfaces only later as a missing key (FisherLens 652eaec).
        if 'cl_dd' in pol_combs and 'lensing' not in param_derivs[name]:
            raise ValueError("The derivatives for %s have no 'lensing' entry, so there is no "
                             'deflection derivative to sum over. The spectrum types present are '
                             '%s. FisherLens reports this omission by printing rather than '
                             'raising, so it is checked here.'
                             % (name, sorted(param_derivs[name].keys())))
        if 'cl_dd' in pol_combs and 'cl_dd' not in param_derivs[name]['lensing']:
            raise ValueError("The 'lensing' derivatives for %s have no 'cl_dd' entry. Its entries "
                             'are %s' % (name, sorted(param_derivs[name]['lensing'].keys())))
        for spectrum in spectrum_types:
            if spectrum not in param_derivs[name]:
                raise ValueError('The derivatives for %s have no %s entry. The spectrum types '
                                 'present are %s'
                                 % (name, spectrum, sorted(param_derivs[name].keys())))
            for spec in pol_combs:
                if spec == 'cl_dd':
                    continue
                if spec not in param_derivs[name][spectrum]:
                    raise ValueError('The %s derivatives for %s have no %s entry. Its entries are '
                                     '%s' % (spectrum, name, spec,
                                             sorted(param_derivs[name][spectrum].keys())))


def _check_length(what, have, needed, lmax):
    """
    Check that an array reaches the highest multipole of a window.

    Parameters
    ----------
    what : str
        Description of the array, for the message.
    have : int
        Number of entries present.
    needed : int
        Number of entries required.
    lmax : int
        Highest multipole being summed to, for the message.

    Raises
    ------
    ValueError
        If there are too few entries.
    """

    if have < needed:
        raise ValueError('Too few multipoles in %s for a sum to %d: %d entries are needed '
                         '(counting from multipole 2) and %d are present.'
                         % (what, lmax, needed, have))


def _warn_non_finite(matrix, spectrum, label):
    """
    Warn if a Fisher matrix holds a non-finite entry.

    Parameters
    ----------
    matrix : ndarray
        The matrix to check.
    spectrum : str
        Spectrum type it was formed for, for the message.
    label : str or None
        Name of the configuration, for the message.

    Warns
    -----
    UserWarning
        If any entry is not finite.
    """

    bad = int( (~np.isfinite(matrix)).sum() )
    if bad:
        warnings.warn('The %s Fisher matrix%s holds %d non-finite entries out of %d. FisherLens '
                      'fills a multipole with them when it cannot invert the covariance there and '
                      'reports that by printing rather than raising, so check its output above.'
                      % (spectrum, '' if label is None else ' for %s' % (label), bad, matrix.size),
                      stacklevel=2)


def _non_gaussian_matrices(fisher_tools,
                           powers,
                           param_derivs,
                           fisher_noise,
                           deflection,
                           vary_params,
                           spectrum_types,
                           pol_combs,
                           windows,
                           lmax_cov,
                           lensing_derivs,
                           unlensed_derivs
                           ):
    r"""
    Form the non-Gaussian Fisher matrices, one spectrum type at a time.

    Parameters
    ----------
    fisher_tools : module
        The imported ``fisherTools``.
    powers : dict
        Spectra, as for :func:`survey_fisher_block`.
    param_derivs : dict
        Parameter derivatives.
    fisher_noise : dict
        Effective noise, already re-indexed to :math:`\ell = 2`.
    deflection : ndarray
        Lensing-reconstruction noise, already re-indexed to :math:`\ell = 2`.
    vary_params : list of str
        Parameters, in the order the rows and columns take.
    spectrum_types : list of str
        Spectrum types to form a matrix for.
    pol_combs : list of str
        Spectra entering the sum.
    windows : dict
        Multipole windows.
    lmax_cov : int
        Highest multipole of the covariance.
    lensing_derivs : dict
        FisherLens ``dCldCLd``.
    unlensed_derivs : dict or None
        FisherLens ``dCldCLu``.

    Returns
    -------
    dict
        Fisher matrices keyed by spectrum type, before the sky fraction is applied.

    Raises
    ------
    ValueError
        If a derivative matrix is not two dimensional or is too small to reach ``lmax_cov``.

    Warns
    -----
    UserWarning
        If the unlensed spectra are among ``spectrum_types`` since their non-Gaussian matrix is the Gaussian one.

    Notes
    -----
    Unlike the Gaussian route, which forms every spectrum type in one call, FisherLens takes a single spectrum type per call here, so one covariance is built and factorized per type.

    The unlensed spectra have no non-Gaussian covariance. They do not depend on the lensing potential, i.e. the lensing derivative vanishes, and the second term is absent to avoid double counting the Gaussian information. Therefore, the non-Gaussian covariance is the same as the Gaussian matrix. It is returned with a warning rather than omitted.
    """

    #Note:
    #A deep copy of ``lensing_derivs`` is passed per call because ``getNonGaussianCov`` writes an identity of its own into the dictionary under ``'cl_dd'``. Reusing one dictionary across calls at different upper bounds would leave that entry at the wrong size. (FisherLens ``652eaec``)

    for name, derivs in (('lensing_derivs', lensing_derivs),
                         ('unlensed_derivs', unlensed_derivs)):
        if derivs is None:
            continue
        for key in sorted(derivs):
            array = np.asarray(derivs[key])
            if array.ndim != 2:
                raise ValueError('%s[%r] must be two dimensional, a spectrum multipole against a '
                                 'deflection multipole, got shape %s' % (name, key, array.shape))
            if min(array.shape) < lmax_cov + 1:
                raise ValueError('%s[%r] has shape %s, too small for a covariance to multipole '
                                 '%d, which needs at least %d entries on each axis'
                                 % (name, key, array.shape, lmax_cov, lmax_cov+1))

    matrices = {}
    for spectrum in spectrum_types:
        if spectrum == 'unlensed':
            warnings.warn('The unlensed spectra carry no non-Gaussian covariance since they do '
                          'not depend on the lensing potential, so their Fisher matrix is the '
                          'Gaussian one and is returned as such.', stacklevel=2)
            matrices.update( fisher_tools.getGaussianCMBFisher(powersFid=powers,
                                                               paramDerivs=param_derivs,
                                                               cmbNoiseSpectra=fisher_noise,
                                                               deflectionNoises=deflection,
                                                               cosmoParams=vary_params,
                                                               spectrumTypes=['unlensed'],
                                                               polCombsToUse=pol_combs,
                                                               ellsToUse=copy.deepcopy(windows)) )
            continue

        ells_to_use = copy.deepcopy(windows)
        ells_to_use['lmaxCov'] = int(lmax_cov)
        inv_cov_dot_derivs, deriv_stack = fisher_tools.choleskyInvCovDotParamDerivsNG(
            powersFid=powers,
            cmbNoiseSpectra=fisher_noise,
            deflectionNoiseSpectra=deflection,
            dCldCLd=copy.deepcopy(lensing_derivs),
            paramDerivs=param_derivs,
            cosmoParams=vary_params,
            dCldCLu=None if unlensed_derivs is None else copy.deepcopy(unlensed_derivs),
            ellsToUse=ells_to_use,
            polCombsToUse=pol_combs,
            spectrumType=spectrum)
        matrices[spectrum] = fisher_tools.getNonGaussianCMBFisher(
            invCovDotParamDerivs=inv_cov_dot_derivs,
            paramDerivStack=deriv_stack,
            cosmoParams=vary_params)

    return matrices


def survey_fisher_block(powers,
                        param_derivs,
                        noise,
                        deflection_noise,
                        vary_params,
                        fsky,
                        windows=None,
                        spectrum_types=None,
                        pol_combs=None,
                        deflection_from_multipole_zero=True,
                        lensing_derivs=None,
                        unlensed_derivs=None,
                        lmax_cov=None,
                        label=None,
                        fisherlens_dir=None
                        ):
    r"""
    Form the Fisher matrix of one survey configuration.

    Parameters
    ----------
    powers : dict
        Spectra from :func:`delensing.run_class`, keyed ``'unlensed'``, ``'lensed'``, ``'delensed'`` and ``'lensing'``, each holding multipoles under ``'l'`` from :math:`\ell = 2`.
    param_derivs : dict
        Derivatives from :func:`delensing.parameter_derivatives`, keyed by parameter, then by spectrum type, then by spectrum, with the deflection derivative under ``['lensing']['cl_dd']``.
    noise : dict
        Effective noise from :func:`delensing.build_effective_noise`, whose entries begin at :math:`\ell = 0`. It is re-indexed to :math:`\ell = 2` here.
    deflection_noise : array_like
        Lensing-reconstruction noise :math:`N_\ell^{dd}` of the configuration, as handed to CLASS.
    vary_params : sequence of str
        Parameters to include, in the order the rows and columns of the returned matrices take.
    fsky : float
        Sky fraction of this configuration, by which the matrices are scaled.
    windows : dict, optional
        Multipole windows from :func:`fisher_windows`. Default is ``None``, which calls it with its own defaults.
    spectrum_types : sequence of str, optional
        Spectrum types to form a matrix for. Default is ``None``, which uses :data:`SPECTRUM_TYPES`.
    pol_combs : sequence of str, optional
        Spectra entering the sum. Default is ``None``, which uses :data:`POL_COMBS`.
    deflection_from_multipole_zero : bool, optional
        Whether ``deflection_noise`` is indexed from multipole zero, as :func:`delensing.deflection_noise_for_class` returns it, in which case its first two entries are dropped. Default is ``True``.
    lensing_derivs : dict, optional
        FisherLens ``dCldCLd``, the derivatives of the spectra with respect to the lensing-deflection power, which switches on the non-Gaussian covariance. Default is ``None``, for a Gaussian covariance.
    unlensed_derivs : dict, optional
        FisherLens ``dCldCLu``, the derivatives of the spectra with respect to the unlensed spectra, adding the second non-Gaussian term. Default is ``None``.
    lmax_cov : int, optional
        Highest multipole of the non-Gaussian covariance. Default is ``None``, which uses the largest upper window bound. See the Notes.
    label : str, optional
        Name recorded with the block, such as the configuration it belongs to. Default is ``None``.
    fisherlens_dir : str, optional
        Location of the FisherLens checkout. Default is ``None``, which uses :data:`delensing.DEFAULT_FISHERLENS_DIR`.

    Returns
    -------
    dict
        A self-describing block, with the matrices under ``'fisher'`` keyed by spectrum type, the parameter order under ``'params'``, and ``'fsky'``, ``'windows'``, ``'pol_combs'``, ``'spectrum_types'``, ``'gaussian'``, ``'kind'`` and ``'label'`` recording how they were formed. ``'kind'`` reads ``'survey'``.

    Raises
    ------
    ValueError
        If ``vary_params`` is empty; if a requested parameter, spectrum type or spectrum is absent from the inputs, if ``'unlensed'`` is absent from ``powers``, if the spectra or noise do not reach the highest window bound, if a window is malformed or if ``fsky`` does not lie in :math:`(0, 1]`.

    Warns
    -----
    UserWarning
        If a returned matrix holds a non-finite entry, which indicates that FisherLens failed to invert the covariance at one or more multipoles. Also if a non-Gaussian matrix is requested for the unlensed spectra, where it is returned as the Gaussian one.

    Notes
    -----
    The sky fraction is applied here because FisherLens contains none. Each configuration therefore carries its own and the blocks of several configurations may be summed directly.

    ``'unlensed'`` must be present in ``powers`` whatever ``spectrum_types`` requests, because ``getGaussianCov`` reads its multipole array unconditionally to set up the covariance; restricting ``spectrum_types`` filters the list of matrices to form and nothing else.

    The noise is re-indexed from :math:`\ell = 0` to :math:`\ell = 2` because ``getGaussianCov`` adds it to the fiducial spectra element by element and those begin at :math:`\ell = 2`, whereas the noise written out for CLASS must begin at :math:`\ell = 0`. One array is built from zero and sliced at each use rather than assembled twice. The deflection noise is treated the same way.

    On the non-Gaussian path, ``lmax_cov`` sets a bound that ``fisherTools`` applies to both axes of the derivative matrices, so it truncates the sum over the deflection multipole as well as the range of the spectrum multipole.

    FisherLens prints to standard output from within these calls, both to report which windows it is using and to report a failed covariance inversion. That output is not suppressed since the second kind matters. A failed inversion leaves non-finite entries, which are checked for here and warned about independently.
    """

    #TODO:
    #The choice of ``lmax_cov`` needs cross-checking with a collaborator before any non-Gaussian number is quoted, and is recorded as an open item. Raising it past the largest window bound extends the sum over the deflection multipole :math:`L`, which is physical and wanted since the covariance of :math:`C_\ell` at :math:`\ell \le \ell_\mathrm{max}` does receive contributions from :math:`L > \ell_\mathrm{max}`. But the two axes share one bound, so raising it also extends the range of :math:`\ell` past the windows, where the derivative vector is zero. Zero-padding the derivative and inverting the larger covariance yields
    #d_w^\mathsf{T}\,(C^{-1})_{ww}\,d_w = d_w^\mathsf{T}\,\bigl(C_{ww} - C_{wo}C_{oo}^{-1}C_{ow}\bigr)^{-1} d_w \;\ge\; d_w^\mathsf{T}\,C_{ww}^{-1}\,d_w\, ,
    #that is more information than inverting the covariance of the band powers actually used, the excluded multipoles acting as a noise monitor. That is justified only if they are measured, which is not the case when the upper bound was chosen to avoid foregrounds and systematics. The default here is the largest window bound, which is the conservative reading; the reference example in FisherLens sets it to the buffered ``lmax_calc`` instead.

    vary_params = list(vary_params)
    if not vary_params:
        raise ValueError('There are no parameters to form a Fisher matrix over. vary_params is '
                         'empty.')
    if len(set(vary_params)) != len(vary_params):
        raise ValueError('vary_params holds a repeated name: %s' % (vary_params))
    if not 0. < float(fsky) <= 1.:
        raise ValueError('fsky must lie in (0, 1], got %s' % (fsky))

    if spectrum_types is None:
        spectrum_types = list(SPECTRUM_TYPES)
    else:
        spectrum_types = list(spectrum_types)
    if not spectrum_types:
        raise ValueError('There are no spectrum types to form a Fisher matrix for. '
                         'spectrum_types is empty.')
    pol_combs = list(POL_COMBS) if pol_combs is None else list(pol_combs)
    if not pol_combs:
        raise ValueError('There are no spectra to sum over. pol_combs is empty.')

    if windows is None:
        windows = fisher_windows()
    else:
        _check_windows(windows)
    for spec in pol_combs:
        if spec not in windows:
            raise ValueError('There is no multipole window for %s. FisherLens raises on a missing '
                             'window; the windows given cover %s'
                             % (spec, sorted(key for key in windows if key != 'lmaxCov')))

    _check_powers(powers, spectrum_types, pol_combs)
    _check_derivs(param_derivs, vary_params, spectrum_types, pol_combs)

    lmax_windows = max( int(windows[spec][1]) for spec in pol_combs )
    fisher_noise = delensing.noise_for_fisher(noise)
    deflection = np.asarray(deflection_noise, dtype=float)
    if deflection_from_multipole_zero:
        if deflection.ndim != 1:
            raise ValueError('The deflection noise must be one dimensional, got shape %s'
                             % (deflection.shape,))
        if len(deflection) < 3:
            raise ValueError('The deflection noise must cover at least multipole 2 when it is '
                             'indexed from zero, got %d entries' % (len(deflection)))
        deflection = deflection[2:]

    needed = lmax_windows - 1
    _check_length('the spectra', len(powers['unlensed']['l']), needed, lmax_windows)
    for spec in pol_combs:
        if spec == 'cl_dd':
            _check_length('the deflection noise', len(deflection), needed, lmax_windows)
        else:
            _check_length('the effective noise', len(fisher_noise[spec]), needed, lmax_windows)

    _, fisher_tools = delensing.import_fisherlens(fisherlens_dir)

    gaussian = lensing_derivs is None
    if gaussian and unlensed_derivs is not None:
        raise ValueError('Derivatives with respect to the unlensed spectra were supplied without '
                         'derivatives with respect to the lensing-deflection power, so there is no '
                         'non-Gaussian covariance to add them to. Supply lensing_derivs as well.')

    matrices = {}
    if gaussian:
        #getGaussianCMBFisher forms every requested spectrum type in one call and returns them as a
        #dictionary; a fresh copy of the windows is passed because createEllsToUseDict mutates it.
        matrices = fisher_tools.getGaussianCMBFisher(powersFid=powers,
                                                     paramDerivs=param_derivs,
                                                     cmbNoiseSpectra=fisher_noise,
                                                     deflectionNoises=deflection,
                                                     cosmoParams=vary_params,
                                                     spectrumTypes=spectrum_types,
                                                     polCombsToUse=pol_combs,
                                                     ellsToUse=copy.deepcopy(windows))
    else:
        if lmax_cov is None:
            lmax_cov = lmax_windows
        matrices = _non_gaussian_matrices(fisher_tools,
                                          powers,
                                          param_derivs,
                                          fisher_noise,
                                          deflection,
                                          vary_params,
                                          spectrum_types,
                                          pol_combs,
                                          windows,
                                          int(lmax_cov),
                                          lensing_derivs,
                                          unlensed_derivs)

    block = {'fisher': {},
             'params': vary_params,
             'fsky': float(fsky),
             'windows': copy.deepcopy(windows),
             'pol_combs': pol_combs,
             'spectrum_types': spectrum_types,
             'gaussian': gaussian,
             'kind': 'survey',
             'label': label}
    for spectrum in spectrum_types:
        matrix = np.asarray(matrices[spectrum], dtype=float) * float(fsky)
        _warn_non_finite(matrix, spectrum, label)
        block['fisher'][spectrum] = matrix
    if not gaussian:
        block['lmax_cov'] = int(lmax_cov)

    return block


def satellite_fisher_block(powers,
                           param_derivs,
                           noise,
                           vary_params,
                           fsky,
                           lmin_sat=2,
                           lmax_sat=29,
                           spectrum=None,
                           pol_combs=None,
                           label=None,
                           fisherlens_dir=None
                           ):
    r"""
    Form the large-scale temperature Fisher block contributed by the satellite.

    Parameters
    ----------
    powers : dict
        Spectra from :func:`delensing.run_class`, as for :func:`survey_fisher_block`.
    param_derivs : dict
        Derivatives from :func:`delensing.parameter_derivatives`.
    noise : dict
        Effective noise from :func:`delensing.build_effective_noise`, whose entries begin at :math:`\ell = 0`. Below the multipole at which the ILC residual is used, this is the satellite noise alone.
    vary_params : sequence of str
        Parameters to include, in the order the rows and columns take.
    fsky : float
        Sky fraction of the satellite, distinct from that of the ground-based configuration.
    lmin_sat : int, optional
        Lowest multipole of the block. Default is 2.
    lmax_sat : int, optional
        Highest multipole of the block, inclusive. Default is 29.
    spectrum : str, optional
        Spectrum type to use. Default is ``None``, which uses :data:`SATELLITE_SPECTRUM`.
    pol_combs : sequence of str, optional
        Spectra entering the block. Default is ``None``, which uses :data:`SATELLITE_POL_COMBS`.
    label : str, optional
        Name recorded with the block. Default is ``None``.
    fisherlens_dir : str, optional
        Location of the FisherLens checkout. Default is ``None``, which uses :data:`delensing.DEFAULT_FISHERLENS_DIR`.

    Returns
    -------
    dict
        A self-describing block, with a single matrix under ``'fisher'``, the parameter order under ``'params'``, and ``'fsky'``, ``'windows'``, ``'pol_combs'``, ``'spectrum'``, ``'gaussian'``, ``'kind'`` and ``'label'`` recording how it was formed. ``'kind'`` reads ``'satellite'``.

    Raises
    ------
    ValueError
        As for :func:`survey_fisher_block` and if the requested window is malformed.

    Warns
    -----
    UserWarning
        If the matrix holds a non-finite entry.

    Notes
    -----
    A single matrix is returned rather than one per spectrum type, because the block uses the lensed spectra whichever spectrum type the ground-based configuration is being evaluated for: delensing is a property of the ground-based measurement and does not apply to the large-scale temperature data of the satellite. :func:`combine_fisher` therefore adds this one matrix to the total of every spectrum type. It is added once to the total and never once per configuration.

    No non-Gaussian option is offered. The non-Gaussian covariance arises from lensing, whose contribution over :math:`\ell \in [2, 29]` in temperature is negligible.
    """

    #Note: Every window is narrowed to the requested range, not only those of ``pol_combs``. ``createEllsToUseDict`` takes the highest upper bound over every key of the dictionary it is given, whatever ``polCombsToUse`` holds, so leaving a wide window in place for a spectrum that is not used still builds the covariance out to that bound. (FisherLens ``652eaec``)
    #The covariance is inverted at the size set by pol_combs, so with its default of temperature alone a 1x1 covariance is inverted at each multipole. That is the likelihood of a temperature-only measurement, which is what the satellite contributes on these scales (and no cross-correlations to other spectra).

    vary_params = list(vary_params)
    if not vary_params:
        raise ValueError('There are no parameters to form a Fisher matrix over. vary_params is '
                         'empty.')
    if not 0. < float(fsky) <= 1.:
        raise ValueError('fsky must lie in (0, 1], got %s' % (fsky))

    spectrum = SATELLITE_SPECTRUM if spectrum is None else spectrum
    pol_combs = list(SATELLITE_POL_COMBS) if pol_combs is None else list(pol_combs)
    if not pol_combs:
        raise ValueError('There are no spectra to sum over. pol_combs is empty.')
    if 'cl_dd' in pol_combs:
        raise ValueError('The satellite block covers the largest scales in temperature and does '
                         'not carry a lensing-deflection term; remove cl_dd from pol_combs.')

    #Every window is narrowed, not only those of pol_combs, because createEllsToUseDict takes the
    #highest upper bound over every key present (FisherLens 652eaec).
    windows = {spec: [int(lmin_sat), int(lmax_sat)] for spec in POL_COMBS}
    _check_windows(windows)

    _check_powers(powers, [spectrum], pol_combs)
    _check_derivs(param_derivs, vary_params, [spectrum], pol_combs)

    fisher_noise = delensing.noise_for_fisher(noise)
    needed = int(lmax_sat) - 1
    _check_length('the spectra', len(powers['unlensed']['l']), needed, int(lmax_sat))
    for spec in pol_combs:
        _check_length('the effective noise', len(fisher_noise[spec]), needed, int(lmax_sat))

    _, fisher_tools = delensing.import_fisherlens(fisherlens_dir)

    #The deflection noise is not read when cl_dd is absent from polCombsToUse, but a correctly
    #shaped array is passed rather than a placeholder so that a future read cannot find a stand-in.
    deflection = np.zeros(len(powers['unlensed']['l']), dtype=float)

    matrices = fisher_tools.getGaussianCMBFisher(powersFid=powers,
                                                 paramDerivs=param_derivs,
                                                 cmbNoiseSpectra=fisher_noise,
                                                 deflectionNoises=deflection,
                                                 cosmoParams=vary_params,
                                                 spectrumTypes=[spectrum],
                                                 polCombsToUse=pol_combs,
                                                 ellsToUse=copy.deepcopy(windows))

    matrix = np.asarray(matrices[spectrum], dtype=float) * float(fsky)
    _warn_non_finite(matrix, spectrum, label)

    block = {'fisher': matrix,
             'params': vary_params,
             'fsky': float(fsky),
             'windows': windows,
             'pol_combs': pol_combs,
             'spectrum': spectrum,
             'gaussian': True,
             'kind': 'satellite',
             'label': label}

    return block


def foreground_systematic(ilc_dic, fractions, spectrum_types=None, pol_combs=None, lmax=None):
    r"""
    Build the spectrum-level systematic implied by mismodeling foreground components.

    Parameters
    ----------
    ilc_dic : dict
        ILC product from :func:`get_fisher_forecasts.load_ilc_product`, whose ``'fg_res_dic'`` holds the residual power of each component.
    fractions : dict
        Fraction of each component that is mismodeled, keyed by component name. A fraction of 1 means the component is not modeled at all, 0 that it is modeled perfectly.
    spectrum_types : sequence of str, optional
        Spectrum types to build an entry for. Default is ``None``, which uses :data:`SPECTRUM_TYPES`.
    pol_combs : sequence of str, optional
        Spectra to build an entry for. Default is ``None``, which uses :data:`POL_COMBS`.
    lmax : int, optional
        Highest multipole needed. Default is ``None``, which uses everything the product reaches.

    Returns
    -------
    dict
        Nested as ``powersFid`` of FisherLens is, keyed by spectrum type and then by spectrum, with the deflection entry under ``'lensing'``. Every array begins at :math:`\ell = 2`.

    Raises
    ------
    ValueError
        If the product has no ``'fg_res_dic'``, if ``fractions`` is empty or names a component the product does not hold, if it names the instrumental noise or if the product does not reach ``lmax``.

    Warns
    -----
    UserWarning
        If a fraction lies outside :math:`[0, 1]`.

    Notes
    -----
    First part of the foreground-bias implementation. The bias on the spectra is

    .. math::

        \Delta C_\ell = \mathbf{w}^T \left(\mathbf{C}_\mathrm{true} - \mathbf{C}_\mathrm{assumed}\right) \mathbf{w}\, ,

    with the ILC weights held at the values the assumed covariance produced rather than re-derived for the true one, which for a multiplicative mismodeling parameterized by the fraction :math:`f_c` per component :math:`c` reduces to a sum over the stored residuals,

    .. math::

        \Delta C_\ell^{XY} = \sum_c f_c\, C_\ell^{\mathrm{res},c,XY}\, , \qquad XY \in \{TT,\, EE\}\, ,

    with :math:`\Delta C_\ell^{TE}` and :math:`\Delta C_\ell^{dd}` zero throughout, the first because ``fg_res_dic`` holds no TE and the second because the deflection carries no ILC residual.

    Arbitrary changes of shape or of spectral energy distribution are not implemented yet: they need ``force_cl_dic``, which exists in :func:`ilc.get_analytic_covariance` but is not plumbed through ``get_covariances`` or ``run_ilc`` at the moment.

    Fractions outside :math:`[0, 1]` are permitted with a warning rather than refused. Above 1 means assuming less of a component than is present, below 0 means over-subtracting it; both are real failure modes and the formalism treats them no differently.

    The same systematic is used for every spectrum type, which may be an approximation depending on the probed scenario.
    The residual enters the observed power additively, and that channel is common to the unlensed, lensed and delensed spectra alike.
    The delensed spectra may carry a second dependence the others do not: the residual is part of the effective noise, so mismodeling it also mismodels the reconstruction noise and therefore how much lensing the delensing removes. Capturing that would mean recomputing the delensed spectra with the perturbed noise for every mismodeling scenario, which we currently do not implement at this stage.
    """

    #Note: The systematic in TE is zero throughout, because ``fg_res_dic`` holds only TT and EE. ``ilc.residual_power`` does form ``cl_residual_te``, but ``WHICH_SPEC_ARR = ['TT', 'EE']`` means it is never requested and it reaches neither ``cl_residual`` nor ``fg_res_dic``.
    #A component is looked up only in the product supplied. A product built without galactic foregrounds carries no ``'galdust'`` or ``'galsync'`` entry, so naming one is an error rather than a fraction of zero and the message says which components are present.

    if 'fg_res_dic' not in ilc_dic:
        raise ValueError("The ILC product holds no 'fg_res_dic', so the residual cannot be "
                         'decomposed into components. That entry is absent whenever the ILC was '
                         'run with save_fg_res_and_weights = 0, which build_output_dic forces '
                         'whenever noise_scalings_for_bands is set. Its keys are %s'
                         % (sorted(ilc_dic.keys())))
    fg_res_dic = ilc_dic['fg_res_dic']
    if not fractions:
        raise ValueError('No mismodeling fractions were given, so there is no systematic to form.')

    spectrum_types = list(SPECTRUM_TYPES) if spectrum_types is None else list(spectrum_types)
    pol_combs = list(POL_COMBS) if pol_combs is None else list(pol_combs)

    available = sorted({key for spec in fg_res_dic for key in fg_res_dic[spec]})
    foregrounds = [name for name in available if name not in NOT_FOREGROUNDS]
    for name in sorted(fractions):
        if name in NOT_FOREGROUNDS:
            raise ValueError('%s is the instrumental noise residual, not a foreground, and cannot '
                             'be given a mismodeling fraction here. The foregrounds this product '
                             'holds are %s.' % (name, ', '.join(foregrounds)))
        if name not in available:
            raise ValueError('The ILC product holds no %s component. It holds %s. A product built '
                             'without galactic foregrounds carries neither galdust nor galsync, so '
                             'check which product this is before assuming the name is misspelled.'
                             % (name, ', '.join(foregrounds)))
        if not 0. <= float(fractions[name]) <= 1.:
            warnings.warn('The mismodeling fraction for %s is %s, outside [0, 1]. Above one means '
                          'assuming less of it than is present and below zero means '
                          'over-subtracting it. Since both can be meaningful, this is allowed.'
                          % (name, fractions[name]), stacklevel=2)

    el = np.asarray(ilc_dic['el'])
    if el[0] != 0:
        raise ValueError('The ILC product multipoles must begin at 0, as get_ilc_residuals writes '
                         'them, got %s.' % (el[0]))
    needed = len(el)-2 if lmax is None else int(lmax)-1
    if len(el)-2 < needed:
        raise ValueError('The ILC product reaches multipole %d, too few for a sum to %s.'
                         % (el[-1], lmax))

    #The residual is stored from multipole 0, as everything on the ILC side is, and is re-indexed
    #to 2 here to match the spectra and the derivatives.
    delta = {}
    for spec, polcomb in (('TT', 'cl_TT'), ('EE', 'cl_EE')):
        if spec not in fg_res_dic:
            continue
        total = np.zeros(needed)
        for name in sorted(fractions):
            if name in fg_res_dic[spec]:
                total = total + float(fractions[name]) \
                    * np.asarray(fg_res_dic[spec][name], dtype=float)[2:needed+2]
        delta[polcomb] = total

    systematic = {}
    for spectrum in spectrum_types:
        systematic[spectrum] = {}
        for polcomb in pol_combs:
            if polcomb == 'cl_dd':
                continue
            systematic[spectrum][polcomb] = delta.get( polcomb, np.zeros(needed) ).copy()
    systematic['lensing'] = {'cl_dd': np.zeros(needed)}

    return systematic


def _masked(array, window, length):
    r"""
    Zero an array outside an inclusive multipole window.

    Parameters
    ----------
    array : array_like
        Values from :math:`\ell = 2`.
    window : sequence of int
        Inclusive bounds ``[lmin, lmax]``.
    length : int
        Number of entries to return.

    Returns
    -------
    ndarray
        Copy of ``array``, zero outside the window.
    """

    values = np.zeros(length)
    low, high = int(window[0])-2, int(window[1])-1
    values[low:high] = np.asarray(array, dtype=float)[low:high]

    return values


def bias_vector(powers,
                param_derivs,
                noise,
                deflection_noise,
                vary_params,
                fsky,
                systematic,
                windows=None,
                spectrum_types=None,
                pol_combs=None,
                deflection_from_multipole_zero=True,
                label=None,
                fisherlens_dir=None
                ):
    r"""
    Form the bias vector one survey configuration contributes.

    Parameters
    ----------
    powers : dict
        Spectra, as for :func:`survey_fisher_block`.
    param_derivs : dict
        Parameter derivatives.
    noise : dict
        Effective noise, indexed from :math:`\ell = 0`.
    deflection_noise : array_like
        Lensing-reconstruction noise.
    vary_params : sequence of str
        Parameters, in the order the entries take.
    fsky : float
        Sky fraction of this configuration.
    systematic : dict
        Spectrum-level systematic from :func:`foreground_systematic`.
    windows : dict, optional
        Multipole windows. Default is ``None``, which calls :func:`fisher_windows`.
    spectrum_types : sequence of str, optional
        Spectrum types to form a vector for. Default is ``None``, which uses :data:`SPECTRUM_TYPES`.
    pol_combs : sequence of str, optional
        Spectra entering the sum. Default is ``None``, which uses :data:`POL_COMBS`.
    deflection_from_multipole_zero : bool, optional
        Whether ``deflection_noise`` is indexed from multipole zero. Default is ``True``.
    label : str, optional
        Name recorded with the vector. Default is ``None``.
    fisherlens_dir : str, optional
        Location of the FisherLens checkout. Default is ``None``, which uses :data:`delensing.DEFAULT_FISHERLENS_DIR`.

    Returns
    -------
    dict
        A self-describing block, with the vectors under ``'bias'`` keyed by spectrum type, the parameter order under ``'params'``, and ``'fsky'``, ``'windows'``, ``'pol_combs'``, ``'spectrum_types'``, ``'kind'`` and ``'label'``. ``'kind'`` reads ``'survey_bias'``.

    Raises
    ------
    ValueError
        As for :func:`survey_fisher_block` and if the systematic lacks an entry the sum needs.

    Notes
    -----
    The sky fraction is applied as it is to a Fisher block. It very nearly cancels from the parameter bias since :math:`\Delta p = F^{-1} b` scales both, but not exactly, because a prior contributes to :math:`F` and not to :math:`b`.
    """

    #Note: ``getBiasVectorGaussian`` of FisherLens (``652eaec``) takes no multipole windows: it has no ``ellsToUse`` argument and makes no call to ``createEllsToUseDict``, so it sums from :math:`\ell = 2` to its ``lmax`` for every spectrum, with no lower bound and no way to give it one, where the Fisher sums over whatever windows :func:`fisher_windows` supplied. Both the systematic and a copy of the derivatives are therefore masked to the windows before being handed over. Masking only one is not equivalent, because the inverse covariance mixes spectra at fixed multipole, so an unmasked derivative at multipoles the windows exclude would still pull on a masked systematic through the off-diagonal terms.
    #The full set of spectra is used, with the systematic zero-filled in :math:`TE` and in the deflection, rather than restricting ``polCombsToUse`` to those that carry one. Restricting it passes straight through to ``getGaussianCov``, which would then invert a sub-block of the covariance rather than take a block of its inverse, and those differ.

    vary_params = list(vary_params)
    if not vary_params:
        raise ValueError('There are no parameters to form a bias vector over. vary_params is '
                         'empty.')
    if not 0. < float(fsky) <= 1.:
        raise ValueError('fsky must lie in (0, 1], got %s' % (fsky))

    spectrum_types = list(SPECTRUM_TYPES) if spectrum_types is None else list(spectrum_types)
    pol_combs = list(POL_COMBS) if pol_combs is None else list(pol_combs)
    if windows is None:
        windows = fisher_windows()
    else:
        _check_windows(windows)

    _check_powers(powers, spectrum_types, pol_combs)
    _check_derivs(param_derivs, vary_params, spectrum_types, pol_combs)

    lmax = max( int(windows[spec][1]) for spec in pol_combs )
    length = lmax-1
    fisher_noise = delensing.noise_for_fisher(noise)
    deflection = np.asarray(deflection_noise, dtype=float)
    if deflection_from_multipole_zero:
        deflection = deflection[2:]

    _check_length('the spectra', len(powers['unlensed']['l']), length, lmax)
    for spec in pol_combs:
        source = deflection if spec == 'cl_dd' else fisher_noise[spec]
        _check_length('the effective noise' if spec != 'cl_dd' else 'the deflection noise',
                      len(source), length, lmax)

    #Both the systematic and the derivatives are masked since getBiasVectorGaussian of
    #FisherLens applies no windows of its own and the inverse covariance mixes spectra at fixed
    #multipole.
    masked_sys, masked_derivs = {}, {}
    for spectrum in spectrum_types:
        masked_sys[spectrum] = {}
        for spec in pol_combs:
            if spec == 'cl_dd':
                continue
            if spec not in systematic.get(spectrum, {}):
                raise ValueError('The systematic has no %s entry for the %s spectra, which the sum '
                                 'needs. Its entries are %s'
                                 % (spec, spectrum, sorted( systematic.get(spectrum, {}).keys() )))
            masked_sys[spectrum][spec] = _masked(systematic[spectrum][spec], windows[spec], length)
    if 'cl_dd' in pol_combs:
        source = systematic.get('lensing', {}).get( 'cl_dd', np.zeros(length) )
        masked_sys['lensing'] = {'cl_dd': _masked(source, windows['cl_dd'], length)}

    for name in vary_params:
        masked_derivs[name] = {}
        for spectrum in spectrum_types:
            masked_derivs[name][spectrum] = {
                spec: _masked(param_derivs[name][spectrum][spec], windows[spec], length)
                for spec in pol_combs if spec != 'cl_dd'}
        if 'cl_dd' in pol_combs:
            masked_derivs[name]['lensing'] = {
                'cl_dd': _masked(param_derivs[name]['lensing']['cl_dd'], windows['cl_dd'], length)}

    _, fisher_tools = delensing.import_fisherlens(fisherlens_dir)
    vectors = fisher_tools.getBiasVectorGaussian(powersFid=powers,
                                                 paramDerivs=masked_derivs,
                                                 cmbNoiseSpectra=fisher_noise,
                                                 deflectionNoises=deflection,
                                                 cosmoParams=vary_params,
                                                 sysSpectrum=masked_sys,
                                                 spectrumTypes=spectrum_types,
                                                 polCombsToUse=pol_combs,
                                                 lmax=lmax)

    block = {'bias': {}, 'params': vary_params, 'fsky': float(fsky),
             'windows': copy.deepcopy(windows), 'pol_combs': pol_combs,
             'spectrum_types': spectrum_types, 'kind': 'survey_bias', 'label': label}
    for spectrum in spectrum_types:
        vector = np.asarray(vectors[spectrum], dtype=float) * float(fsky)
        _warn_non_finite(vector, spectrum, label)
        block['bias'][spectrum] = vector

    return block


# Totals

def prior_matrix(params, priors=None):
    r"""
    Build the Fisher matrix contributed by independent Gaussian priors.

    Parameters
    ----------
    params : sequence of str
        Parameters, in the order the rows and columns take.
    priors : dict, optional
        Gaussian prior widths keyed by parameter, as standard deviation :math:`\sigma` rather than as a variance :math:`\sigma^2`. Default is ``None``, which returns a zero matrix.

    Returns
    -------
    ndarray
        A matrix of the same shape as the Fisher matrices, holding :math:`1/\sigma^2` on the diagonal of each parameter that has a prior.

    Raises
    ------
    ValueError
        If a prior is given for a parameter that is not among ``params`` or if a width is not positive.

    Notes
    -----
    An independent Gaussian prior of width :math:`\sigma` on a parameter contributes :math:`1/\sigma^2` to its own diagonal element and nothing elsewhere, so priors are added to the summed Fisher matrix rather than to any one block.
    """

    #Note: A separate matrix is built and added rather than the diagonal being modified in place, both so that the caller's matrices are left alone and so that a correlated prior could be substituted here without changing anything else.

    params = list(params)
    matrix = np.zeros(( len(params), len(params) ))
    if not priors:
        return matrix

    for name in sorted(priors):
        if name not in params:
            raise ValueError('There is a prior on %s, which is not among the parameters %s. A '
                             'prior on a parameter that is not varied has nowhere to act.'
                             % (name, ', '.join(params)))
        width = float(priors[name])
        if not width > 0.:
            raise ValueError('The prior width for %s must be positive, got %s. It is a standard '
                             'deviation, not a variance.' % (name, priors[name]))
        index = params.index(name)
        matrix[index, index] = 1./width**2

    return matrix


def invert_fisher(matrix, params, spectrum=None):
    r"""
    Invert a Fisher matrix through a diagonal rescaling.

    Parameters
    ----------
    matrix : ndarray
        The Fisher matrix.
    params : sequence of str
        Parameter names, for the messages.
    spectrum : str, optional
        Spectrum type, for the messages. Default is ``None``.

    Returns
    -------
    covariance : ndarray
        The inverse.
    condition : float
        Condition number of the rescaled matrix, which is the meaningful one.

    Raises
    ------
    ValueError
        If a diagonal element is not positive or if the matrix cannot be inverted.

    Warns
    -----
    UserWarning
        If the rescaled condition number is large enough that the inverse should not be trusted.

    Notes
    -----
    The matrix is rescaled by its own diagonal before inversion and the result scaled back, which is exact and changes no number that comes out of it, but makes the arithmetic well behaved and the condition number meaningful.
    """

    where = '' if spectrum is None else 'the %s ' % (spectrum)
    diagonal = np.diag(matrix)
    if np.any(diagonal <= 0.):
        blind = [name for name, value in zip(params, diagonal) if value <= 0.]
        raise ValueError('In %stotal the diagonal element of %s is not positive, so the data '
                         'constrain it not at all and no uncertainty can be formed. Give it a '
                         'prior or hold it fixed.' % (where, ', '.join(blind)))

    #Without the rescaling the condition number is dominated by the units of the parameters rather than by anything about the measurement.
    scale = np.sqrt(diagonal)
    outer = np.outer(scale, scale)
    rescaled = matrix/outer
    try:
        covariance = np.linalg.inv(rescaled)/outer
    except np.linalg.LinAlgError as err:
        raise ValueError('The %stotal cannot be inverted: %s. A singular total usually means a '
                         'direction in parameter space the data do not constrain at all, which '
                         'needs either a prior or a parameter held fixed.' % (where, err)) from err

    condition = float( np.linalg.cond(rescaled) )
    if condition > 1.e12:
        warnings.warn('The %stotal has a rescaled condition number of %.2e, so some combination '
                      'of the parameters is very nearly unconstrained and the uncertainties that '
                      'follow should not be trusted.' % (where, condition), stacklevel=3)

    return covariance, condition


def combine_fisher(blocks, satellite=None, priors=None, fix_params=None, spectrum_types=None,
                   fisherlens_dir=None):
    r"""
    Sum Fisher blocks into a total and reduce it to projected uncertainties.

    Parameters
    ----------
    blocks : sequence of dict
        Survey blocks from :func:`survey_fisher_block`, one per configuration. Each is already scaled by its own sky fraction.
    satellite : dict, optional
        The satellite block from :func:`satellite_fisher_block`, added once to the total of every spectrum type. Default is ``None``.
    priors : dict, optional
        Gaussian prior widths keyed by parameter, applied once to the summed matrix. Default is ``None``.
    fix_params : sequence of str, optional
        Parameters to hold fixed, whose rows and columns are removed after the priors are added. Default is ``None``.
    spectrum_types : sequence of str, optional
        Spectrum types to form a total for. Default is ``None``, which uses those every block has.
    fisherlens_dir : str, optional
        Location of the FisherLens checkout. Default is ``None``, which uses :data:`delensing.DEFAULT_FISHERLENS_DIR`.

    Returns
    -------
    dict
        ``'fisher'``, ``'covariance'`` and ``'condition'`` keyed by spectrum type, ``'sigmas'`` keyed by spectrum type and then by parameter, the surviving parameter order under ``'params'``, and ``'params_supplied'``, ``'params_fixed'``, ``'priors'``, ``'spectrum_types'``, ``'survey_labels'`` and ``'satellite_label'`` recording how the total was formed.

    Raises
    ------
    ValueError
        If no blocks are given, if a satellite block appears among ``blocks``, if the blocks do not agree on the parameters or their order, if a requested spectrum type is missing from a block, if every parameter is fixed, or if a total cannot be inverted.

    Warns
    -----
    UserWarning
        If two blocks carry the same label, which usually means one configuration was added twice, if Gaussian and non-Gaussian blocks are mixed, or if a prior is placed on a parameter that is then fixed, where it has no effect.

    Notes
    -----
    This is the only route to a projected uncertainty and the order of operations is what makes it correct.
    Each block arrives already scaled by its own sky fraction since the configurations and sky patches are treated as independent experiments whose Fisher matrices simply add.
    The blocks are summed, the satellite is added once, the priors are added once to that sum, the fixed parameters are removed and only then is the result inverted.

    Adding a prior to an individual block instead would count it once per configuration; adding the satellite inside the loop would do the same. Both are the reason this function takes the satellite separately rather than as one more entry in ``blocks`` and rejects a satellite block that appears in the list.
    """

    blocks = list(blocks)
    if not blocks:
        raise ValueError('There are no Fisher blocks to combine. combine_fisher needs at least one '
                         'survey block; the satellite block is passed separately.')

    for position, block in enumerate(blocks):
        kind = block.get('kind')
        if kind == 'satellite':
            raise ValueError('Block %d is a satellite block. It is passed as the satellite '
                             'argument rather than among blocks so that it is added once to the '
                             'total rather than once per configuration.' % (position))
        if kind != 'survey':
            raise ValueError("Block %d has kind %r, expected 'survey'. Blocks come from "
                             'survey_fisher_block.' % (position, kind))
    if satellite is not None and satellite.get('kind') != 'satellite':
        raise ValueError("The satellite argument has kind %r, expected 'satellite'. It comes from "
                         'satellite_fisher_block.' % (satellite.get('kind')))

    labels = [block.get('label') for block in blocks]
    named = [label for label in labels if label is not None]
    if len(set(named)) != len(named):
        repeated = sorted({label for label in named if named.count(label) > 1})
        warnings.warn('More than one block is labeled %s. Blocks are summed, so a configuration '
                      'included twice is counted twice.' % (', '.join(repeated)), stacklevel=2)
    gaussian = {bool( block.get('gaussian', True) ) for block in blocks}
    if len(gaussian) > 1:
        warnings.warn('Gaussian and non-Gaussian blocks are being summed. They describe the same '
                      'measurement under different assumptions about the covariance, so their sum '
                      'describes neither.', stacklevel=2)

    params = list(blocks[0]['params'])
    for position, block in enumerate(blocks[1:], start=1):
        if list(block['params']) != params:
            raise ValueError('Block %d has parameters %s, block 0 has %s. Summing them would add '
                             'unlike rows and columns; every block must vary the same parameters '
                             'in the same order.'
                             % (position, ', '.join(block['params']), ', '.join(params)))
    if satellite is not None and list(satellite['params']) != params:
        raise ValueError('The satellite block has parameters %s, the survey blocks have %s.'
                         % (', '.join(satellite['params']), ', '.join(params)))

    common = set(blocks[0]['fisher'])
    for block in blocks[1:]:
        common &= set(block['fisher'])
    if spectrum_types is None:
        spectrum_types = [name for name in SPECTRUM_TYPES if name in common]
        if not spectrum_types:
            spectrum_types = sorted(common)
    else:
        spectrum_types = list(spectrum_types)
        for name in spectrum_types:
            if name not in common:
                raise ValueError('Not every block has a %s matrix. The types common to all of them '
                                 'are %s.' % (name, ', '.join( sorted(common) ) or 'none'))
    if not spectrum_types:
        raise ValueError('The blocks have no spectrum type in common, so there is nothing to sum.')

    fix_params = [] if fix_params is None else list(fix_params)
    for name in fix_params:
        if name not in params:
            raise ValueError('%s is to be held fixed but is not among the parameters %s.'
                             % (name, ', '.join(params)))
    if len(set(fix_params)) != len(fix_params):
        raise ValueError('fix_params holds a repeated name: %s' % (', '.join(fix_params)))
    if len(fix_params) >= len(params):
        raise ValueError('Every parameter is to be held fixed, which leaves no Fisher matrix to '
                         'invert.')
    for name in sorted(priors or {}):
        if name in fix_params:
            warnings.warn('There is a prior on %s, which is also being held fixed. Its row and '
                          'column are removed after the prior is added, so the prior has no '
                          'effect.' % (name), stacklevel=2)

    priors_matrix = prior_matrix(params, priors)
    _, fisher_tools = delensing.import_fisherlens(fisherlens_dir)

    totals, covariances, sigmas, conditions = {}, {}, {}, {}
    surviving = params
    for spectrum in spectrum_types:
        total = np.zeros(( len(params), len(params) ))
        for block in blocks:
            total = total + np.asarray(block['fisher'][spectrum], dtype=float)
        if satellite is not None:
            total = total + np.asarray(satellite['fisher'], dtype=float)
        total = total + priors_matrix

        if fix_params:
            #fixParametersInFisher indexes and removes from the list it is given, so a plain list
            #is required and an array will not do (FisherLens 652eaec).
            total, surviving = fisher_tools.fixParametersInFisher(total, list(params), fix_params,
                                                                  returnFixedParamList=True)
            surviving = list(surviving)
        _warn_non_finite(total, spectrum, 'the total')

        covariance, condition = invert_fisher(total, surviving, spectrum)
        conditions[spectrum] = condition
        variances = np.diag(covariance)
        if np.any(variances <= 0.):
            bad = [name for name, value in zip(surviving, variances) if value <= 0.]
            raise ValueError('The %s covariance has a non-positive variance for %s, so the total '
                             'is not positive definite and no uncertainty can be read from it.'
                             % (spectrum, ', '.join(bad)))

        totals[spectrum] = total
        covariances[spectrum] = covariance
        sigmas[spectrum] = {name: float(value)
                            for name, value in zip( surviving, np.sqrt(variances) )}

    return {'fisher': totals,
            'covariance': covariances,
            'sigmas': sigmas,
            'condition': conditions,
            'params': list(surviving),
            'params_supplied': params,
            'params_fixed': fix_params,
            'priors': dict(priors or {}),
            'spectrum_types': spectrum_types,
            'survey_labels': labels,
            'satellite_label': None if satellite is None else satellite.get('label', 'satellite')}


def parameter_bias(combined, vectors):
    r"""
    Reduce bias vectors to the displacement they induce in each parameter: :math:`\Delta p = F^{-1} b`.

    Parameters
    ----------
    combined : dict
        Total from :func:`combine_fisher`, whose covariance and parameter bookkeeping are used.
    vectors : sequence of dict
        Bias vectors from :func:`bias_vector`, one per configuration, each already scaled by its own sky fraction.

    Returns
    -------
    dict
        ``'bias'`` keyed by spectrum type and then by parameter, ``'bias_in_sigma'`` the same divided by the projected uncertainty, and ``'params'``, ``'fixed'``, ``'spectrum_types'`` and ``'labels'`` recording how it was formed.

    Raises
    ------
    ValueError
        If no vectors are given, if a vector is not a bias vector, if the vectors and the total do not agree on the parameters or their order, or if a spectrum type is missing from either.

    Warns
    -----
    UserWarning
        Always since the foreground-bias calculation has a limited implementation, is experimental and has not been checked against an independent calculation.

    Notes
    -----
    The displacement is :math:`\Delta p = F^{-1} b`, with the same total :math:`F` that produced the uncertainties, so the satellite block and the priors enter here through the covariance exactly as they do there.

    The satellite contributes no bias vector: Below the multipole at which the ILC residual is used the noise is the satellite's own, carrying no ILC residual to be mismodeled, so its contribution to :math:`b` is zero while its contribution to :math:`F` is not.

    The displacement is reported in units of the projected uncertainty as well as in the units of the parameter since whether a bias matters is a question about their ratio.
    """

    vectors = list(vectors)
    if not vectors:
        raise ValueError('There are no bias vectors to combine.')
    for position, vector in enumerate(vectors):
        if vector.get('kind') != 'survey_bias':
            raise ValueError("Bias vector %d has kind %r, expected 'survey_bias'. Bias vectors "
                             'come from bias_vector.' % (position, vector.get('kind')))

    #The bias is labeled experimental in the parameter file, so the warning names what that means
    #rather than leaving the reader to guess.
    warnings.warn('The foreground-bias calculation has a limited implementation, is experimental '
                  'and has not been checked against an independent calculation. Treat the numbers '
                  'as indicative.', stacklevel=2)

    supplied = list( combined['params_supplied'] )
    for position, vector in enumerate(vectors):
        if list(vector['params']) != supplied:
            raise ValueError('Bias vector %d has parameters %s, the total was formed over %s. '
                             'Summing them would add unlike entries.'
                             % (position, ', '.join(vector['params']), ', '.join(supplied)))

    spectrum_types = [name for name in combined['spectrum_types']
                      if all( name in vector['bias'] for vector in vectors )]
    if not spectrum_types:
        raise ValueError('The total and the bias vectors have no spectrum type in common. The '
                         'total holds %s and the vectors hold %s.'
                         % (', '.join(combined['spectrum_types']),
                            ', '.join( sorted(vectors[0]['bias'].keys()) )))

    fixed = list( combined.get('params_fixed', []) )
    keep = [supplied.index(name) for name in combined['params']]

    bias, in_sigma = {}, {}
    for spectrum in spectrum_types:
        total = np.zeros( len(supplied) )
        for vector in vectors:
            total = total + np.asarray(vector['bias'][spectrum], dtype=float)
        #A parameter held fixed cannot be displaced, so its entry is dropped alongside its row and
        #column of the total.
        displacement = combined['covariance'][spectrum] @ total[keep]
        bias[spectrum] = {name: float(value)
                          for name, value in zip(combined['params'], displacement)}
        in_sigma[spectrum] = {name: float(value)/combined['sigmas'][spectrum][name]
                              for name, value in zip(combined['params'], displacement)}

    return {'bias': bias,
            'bias_in_sigma': in_sigma,
            'params': list( combined['params'] ),
            'fixed': fixed,
            'spectrum_types': spectrum_types,
            'labels': [vector.get('label') for vector in vectors]}
