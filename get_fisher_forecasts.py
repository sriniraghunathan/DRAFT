r"""
Compute delensed CMB spectra and Fisher forecasts from ILC residuals.

This is the driver for the forecasting stage of DRAFT.
It reads the residual power spectra produced by the component separation in :mod:`get_ilc_residuals`, combines them with the effective satellite noise, drives the iterative delensing of CLASS_delens through FisherLens, and reduces the result to a Fisher matrix and projected parameter uncertainties.

Several ILC products are treated as independent experiments whose Fisher matrices are summed, which is how different surveys are combined into a single forecast.

Run from the repository root (example for the CMB-S4 Wide and CMB-S4 Ultra-deep surveys of the conceptual/preliminary baseline design)::

    python3 get_fisher_forecasts.py -ilc_fname products/202310xx_PBDR_config/s4wide_202310xx_pbdr_config/s4wide_202310xx_pbdr_config_ilc_galaxy1_27-39-93-145-225-278_TT-EE_galmask2_AZ_lmax6500_for7years.npy products/202310xx_PBDR_config/s4deepv3r025_202310xx_pbdr_config/s4deepv3r025_202310xx_pbdr_config_ilc_galaxy0_20-27-39-93-145-225-278_TT-EE_lmax6500_for7years.npy -fsky 0.59 0.03 -label s4_conceptual -write_cache 1

Only what does not bear on a number is given on the command line. The cosmology, the step sizes, the multipole ranges, the delensing mode, the priors and the option flags are read from ``params_fisher.ini`` so that a forecast is reproducible from a file rather than from a shell history.

The same calculation is available to other scripts through :func:`run_forecast`, which accepts the command line options as keyword arguments and returns the forecast together with the path it would be written to::

    import get_fisher_forecasts

    product, opfname = get_fisher_forecasts.run_forecast(['products/.../s4wide_..._for7years.npy', 'products/.../s4deepv3r025_..._for7years.npy'], fsky=[0.59, 0.03], label='s4_conceptual', save=False)

Each configuration costs one iterative CLASS_delens run plus two per varied parameter since the delensed spectra and their derivatives depend on the noise of that configuration.
Those derivatives are therefore the expensive part of any forecast and :func:`write_forecast_cache` records them so that ``-cache`` can read them back: a rerun at a different sky fraction, with different priors or with a parameter held fixed then costs seconds and no CLASS run at all.

Nothing need be written to disk: :func:`get_ilc_residuals.run_ilc` returns the same dictionary that :func:`load_ilc_product` reads back, so a script can run component separation and forecasting in one process::

    import get_ilc_residuals

    ilc_dic, _ = get_ilc_residuals.run_ilc('s4wide_202310xx_pbdr_config', include_gal=1, which_gal_mask=2, save=False)
    product, _ = get_fisher_forecasts.run_forecast(ilc_dic, fsky=0.59, label='s4wide', save=False)

Notes
-----
Effective noise assembly and the CLASS_delens calls live in :mod:`delensing`, and the Fisher algebra in :mod:`fisher`.
This module holds only the input and output handling, the parameter file and the command line so that neither of those modules depends on a file format.
"""

import argparse
import hashlib
import os
import sys
import warnings

import numpy as np

sys.path.append( os.path.join( os.path.dirname(os.path.abspath(__file__)), 'modules' ) )
import delensing
import fisher
import misc


# Constants

PARAMFILE = os.path.join( os.path.dirname(os.path.abspath(__file__)), 'params_fisher.ini' )
"""
Default parameter file for the forecasting stage, kept separate from the ``params.ini`` that configures the component separation.
"""

FID_PREFIX = 'fid_'
"""
Prefix marking a fiducial cosmological parameter value in the parameter file.
"""

STEP_PREFIX = 'step_'
"""
Prefix marking a derivative step size in the parameter file.
"""

#A component not named is taken as perfectly modeled. Components are looked up in the ILC product rather than against a fixed list since a product built without galactic foregrounds carries neither ``galdust`` nor ``galsync``.
FG_BIAS_PREFIX = 'fg_bias_'
"""
Prefix marking a foreground mismodeling fraction in the parameter file.
"""

PRIOR_PREFIX = 'prior_'
"""
Prefix marking a Gaussian prior width in the parameter file.
"""

#``Nl_ET`` is the temperature-polarization estimator that CLASS_delens and everything downstream call ``TE``.
#``Nl_MVpol`` has no counterpart, being the minimum-variance combination of the polarization estimators alone.
#There is no ``Nl_BB``, which the current products hold, although every entry of it is indeterminate.
LEGACY_LENSING_ESTIMATORS = {'Nl_TT': 'TT', 'Nl_EE': 'EE', 'Nl_ET': 'TE', 'Nl_TB': 'TB',
                             'Nl_EB': 'EB', 'Nl_MV': 'MV', 'Nl_MVpol': 'MVpol'}
"""
Estimator names in the lensing-noise products under ``products/``, against the names used since.
"""

#These are the cuts on the CMB fields entering the quadratic estimator, not the windows of the Fisher sum. The current settings use ``lmin``, ``lmax`` and ``lmax_TT`` for the Fisher windows and the names above for the reconstruction, so carrying the older names across unchanged would silently relabel one as the other.
LEGACY_LENSING_SETTINGS = {'lmin': 'lmin_lensing', 'lmax': 'lmax_P_lensing',
                           'lmax_tt': 'lmax_T_lensing'}
"""
Multipole cuts recorded in the older lensing-noise products, against the names used since.
"""

#Every file a run produces lives together there, because a forecast combining several configurations belongs to no single ILC product and has nowhere else to sit.
FORECAST_DIR = os.path.join( os.path.dirname(os.path.abspath(__file__)), 'results', 'forecasts' )
"""
Folder the forecast products are written under, one subfolder per forecast.
"""

FSKY_AUTO = 'auto'
"""
Value of ``-fsky`` asking for the sky fraction recorded in that ILC product rather than a supplied one.
"""

CACHE_NONE = 'none'
"""
Value of ``-cache`` asking for CLASS to be run for that configuration rather than a cache read.
"""

#The lensing convolution mixes multipoles, so spectra at ``lmax`` are only accurate if CLASS computed some way beyond it. Nothing checks how far and a buffer of a few tens is as wrong as one of zero while being far likelier as a typo, so the guard is a fraction of ``lmax`` rather than a test against zero.
#TODO
LBUFFER_FRACTION = 0.05
"""
Fraction of ``lmax`` below which ``lbuffer`` is reported as leaving too little margin (guard against a mistyped parameter file rather than a convergence criterion).
"""


# Command line

def parse_args(argv=None):
    """
    Parse the command line.

    Parameters
    ----------
    argv : list of str, optional
        Argument list. Default is ``sys.argv[1:]``.

    Returns
    -------
    argparse.Namespace
        Parsed arguments. The attribute names match the keyword arguments of :func:`run_forecast`.

    Notes
    -----
    Only the arguments that do not change a number appear here: which products to read, how much sky each covers, where the output goes and how loud the run is.
    Everything that bears on the physics, meaning the cosmology, the step sizes, the multipole ranges, the delensing mode, the priors and the option flags, is read from the parameter file so that a forecast is reproducible from a file rather than from a shell history.
    """

    parser = argparse.ArgumentParser(description='Compute delensed CMB spectra and Fisher forecasts from ILC residuals.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-ilc_fname', dest='ilc_fname', action='store', help='ILC product to forecast from, as written by get_ilc_residuals.py. Several are treated as independent experiments and combined.', nargs='+', type=str, required=True)
    parser.add_argument('-fsky', dest='fsky', action='store', help="Sky fraction of each configuration, one value per -ilc_fname, either a number or 'auto' to use the fsky_val recorded in that product. A product built with -include_gal 0 records none, so a number is required for it.", nargs='+', type=str, default=None)
    parser.add_argument('-cache', dest='cache', action='store', help="Forecast cache to read for each configuration in place of running CLASS, one value per -ilc_fname, with 'none' to run CLASS for that one. For example use to rerun with a different sky fraction or with different priors.", nargs='+', type=str, default=None)
    parser.add_argument('-write_cache', dest='write_cache', action='store', help='Write the forecast cache of each configuration so that a later run can reuse it (0 or 1).', type=int, default=0)
    parser.add_argument('-write_lensing_noise', dest='write_lensing_noise', action='store', help='Write the lensing-reconstruction noise curves of each configuration (0 or 1). Ignored for a configuration read from a cache that holds no reconstruction.', type=int, default=0)
    parser.add_argument('-label', dest='label', action='store', help='Name of this forecast, used for the output folder and the file names within it. Required when more than one ILC product is given; with one, it defaults to that product name.', type=str, default=None)
    parser.add_argument('-opfname', dest='opfname', action='store', help='Path of the forecast product, overriding the one built from -label.', type=str, default=None)
    parser.add_argument('-save', dest='save', action='store', help='Write the forecast product (0 or 1).', type=int, default=1)
    parser.add_argument('-overwrite', dest='overwrite', action='store', help='Replace existing output files rather than refusing to (0 or 1).', type=int, default=0)
    parser.add_argument('-paramfile', dest='paramfile', action='store', help='Forecasting parameter file to read.', type=str, default=PARAMFILE)
    parser.add_argument('-verbose', dest='verbose', action='store', help='Report each step as it runs (0 or 1).', type=int, default=1)
    #Note: get_ilc_residuals.py carries -interactive_mode for diagnostic plots of the beams and noise. An equivalent could be introduced here to plot the reconstruction noise curves and the effective noise against the fiducial spectra, for instance.

    return parser.parse_args(argv)


# Parameter file

def check_lbuffer(lbuffer, lmax, fraction=LBUFFER_FRACTION):
    r"""
    Warn when the buffer above ``lmax`` leaves too little margin for the lensing convolution.

    Parameters
    ----------
    lbuffer : int
        Multipoles CLASS is asked to compute beyond ``lmax``.
    lmax : int
        Highest multipole entering the Fisher sum.
    fraction : float, optional
        Fraction of ``lmax`` below which the buffer is reported. Default is :data:`LBUFFER_FRACTION`.

    Returns
    -------
    int
        The buffer below which the warning fires, for reporting.

    Raises
    ------
    ValueError
        If ``lbuffer`` is negative or ``lmax`` is not positive.

    Warns
    -----
    UserWarning
        If ``lbuffer`` falls below ``fraction`` times ``lmax``.

    Notes
    -----
    Lensing mixes neighboring multipoles, so :math:`C_\ell` at :math:`\ell_\mathrm{max}` is only accurate if CLASS computed some way past it.
    :func:`delensing.check_ilc_product` already refuses an ILC product that does not reach ``lmax + lbuffer``, but the buffer itself is not checked, so spectra that are wrong at the top of the range may be calculated incorrectly with no complaint anywhere.

    The threshold is a guard against a mistyped parameter file and not a convergence criterion.
    """

    #TODO: Test and make it a convergence criterion

    lbuffer, lmax = int(lbuffer), int(lmax)
    if lbuffer < 0:
        raise ValueError('lbuffer must not be negative, got %d' % (lbuffer))
    if lmax <= 0:
        raise ValueError('lmax must be positive, got %d' % (lmax))

    floor = int(fraction*lmax)
    if lbuffer < floor:
        warnings.warn('lbuffer = %d leaves little margin above lmax = %d for the lensing '
                      'convolution, which mixes neighboring multipoles and therefore needs '
                      'spectra computed some way beyond the highest one used. Below about %d, '
                      'which is %g of lmax, the spectra near lmax should not be trusted. This is '
                      'a guard against a mistyped parameter file rather than a measured bound.'
                      % (lbuffer, lmax, floor, fraction), stacklevel=2)

    return floor


def load_fisher_params(paramfile=PARAMFILE):
    r"""
    Read and validate the forecasting parameter file.

    The file uses the same flat ``key = value`` dialect as ``params.ini`` and is read with :func:`misc.get_param_dict`.
    Three groups of entries are collected by prefix rather than by section since that dialect has none: ``fid_<name>`` gives a fiducial cosmological parameter, ``step_<name>`` its derivative step size and ``prior_<name>`` the width of a Gaussian prior on it.
    Everything else is passed through as a setting.

    Parameters
    ----------
    paramfile : str, optional
        Path to the parameter file. Default is :data:`PARAMFILE`.

    Returns
    -------
    dict
        Settings, with the prefixed groups gathered into the entries ``'cosmo_fid'``, ``'step_sizes'`` and ``'priors'``, the varied parameters as a list under ``'vary_params'``, and every remaining entry of the file left as it was read.

    Raises
    ------
    FileNotFoundError
        If the parameter file does not exist.
    ValueError
        If ``vary_params`` is absent or propagated from :func:`delensing.validate_cosmology`.
        Also if the non-Gaussian covariance is requested without saying how the derivative matrices are to be obtained.

    Warns
    -----
    UserWarning
        Propagated from :func:`delensing.validate_cosmology` and if a prior names a parameter that is not being varied since such a prior has no effect.

    Notes
    -----
    Cosmological values and step sizes are coerced to ``float`` by :func:`delensing.validate_cosmology`.
    That matters because :func:`misc.get_param_dict` returns an ``int`` for any value written without a fractional part, so ``fid_N_eff = 3.0`` would otherwise arrive as an integer.
    Coercing here rather than asking for a particular way of writing the file keeps that detail out of the user's way.

    Multipole ranges, the lensing reconstruction cuts, the satellite block and the option flags are returned unchanged since they are consumed by :mod:`fisher` and by :func:`delensing.run_class` rather than here.
    """

    if not os.path.exists(paramfile):
        raise FileNotFoundError('No parameter file at %s' % (paramfile))

    raw = misc.get_param_dict(paramfile)

    cosmo_fid, step_sizes, priors, settings = {}, {}, {}, {}
    for key, value in raw.items():
        if key.startswith(FID_PREFIX):
            cosmo_fid[key[len(FID_PREFIX):]] = value
        elif key.startswith(STEP_PREFIX):
            step_sizes[key[len(STEP_PREFIX):]] = value
        elif key.startswith(PRIOR_PREFIX):
            priors[key[len(PRIOR_PREFIX):]] = float(value)
        else:
            settings[key] = value

    if 'vary_params' not in settings:
        raise ValueError('%s has no vary_params entry, so there is nothing to differentiate. Give it '
                         'as a space separated list of names.' % (paramfile))
    vary_params = str(settings.pop('vary_params')).split()

    cosmo_fid, step_sizes, vary_params = delensing.validate_cosmology(
        cosmo_fid, step_sizes=step_sizes, vary_params=vary_params)

    unused_priors = [name for name in priors if name not in vary_params]
    if unused_priors:
        warnings.warn('These priors name parameters that are not being varied, so they have no '
                      'effect: %s' % (', '.join(sorted(unused_priors))), stacklevel=2)

    fractions = {key[len(FG_BIAS_PREFIX):]: float(value) for key, value in settings.items()
                 if key.startswith(FG_BIAS_PREFIX)}
    for key in list(settings):
        if key.startswith(FG_BIAS_PREFIX):
            del settings[key]
    settings['fg_bias_fractions'] = fractions
    if int( settings.get('do_foreground_bias', 0) ):
        if not fractions:
            raise ValueError('%s asks for the foreground bias but names no component to mismodel. '
                             'Add one entry per component as %s<component> = <fraction>, where a '
                             'fraction of 1 means the component is not modeled at all. The '
                             'components available depend on the ILC product; one built without '
                             'galactic foregrounds carries neither galdust nor galsync.'
                             % (paramfile, FG_BIAS_PREFIX))
    elif fractions:
        warnings.warn('Foreground mismodeling fractions are set for %s but do_foreground_bias is '
                      '0, so they have no effect.' % (', '.join( sorted(fractions) )), stacklevel=2)

    if int( settings.get('do_non_gaussian', 0) ):
        if 'ng_derivatives' not in settings:
            raise ValueError(
                '%s asks for the non-Gaussian covariance but does not say how to obtain the '
                'derivative matrices. Set ng_derivatives to one of %s. There is no default '
                'because they trade a large file against a long calculation: at lmax %s plus a '
                'buffer of %s, a set is %.1f GB per derivative type.'
                % (paramfile, ', '.join(delensing.NG_DERIVATIVE_MODES), settings.get('lmax'),
                   settings.get('lbuffer'),
                   delensing.ng_derivative_bytes( int(settings.get('lmax', 5000))
                                                  + int(settings.get('lbuffer', 500)) )/1e9))
        mode = str(settings['ng_derivatives'])
        if mode not in delensing.NG_DERIVATIVE_MODES:
            raise ValueError('ng_derivatives must be one of %s, got %r'
                             % (', '.join(delensing.NG_DERIVATIVE_MODES), mode))
        settings['ng_derivatives'] = mode
    elif 'ng_derivatives' in settings:
        warnings.warn('ng_derivatives is set but do_non_gaussian is 0, so it has no effect.',
                      stacklevel=2)

    if 'lbuffer' in settings and 'lmax' in settings:
        check_lbuffer(settings['lbuffer'], settings['lmax'])

    settings['cosmo_fid'] = cosmo_fid
    settings['step_sizes'] = step_sizes
    settings['priors'] = priors
    settings['vary_params'] = vary_params

    return settings


# Input

def _load_npy_dict(fname, what='product'):
    r"""
    Read a dictionary written with :func:`numpy.save`.

    Parameters
    ----------
    fname : str
        Path to the ``.npy`` file.
    what : str, optional
        What the file is, used in the messages. Default is ``'product'``.

    Returns
    -------
    dict
        Contents of the file.

    Raises
    ------
    FileNotFoundError
        If the file does not exist.
    ValueError
        If the file does not hold a single dictionary.

    Warns
    -----
    UserWarning
        If the file can only be read with ``encoding='latin1'``, which means it was written under Python 2.

    Notes
    -----
    The default ASCII decoding is tried first and ``encoding='latin1'`` only as a fallback, rather than passing ``latin1`` unconditionally.
    """

    #Note: Every product of this stage is written the same way, with :func:`numpy.save` applied to a dictionary, so they are all read back the same way: with ``allow_pickle`` and unwrapped from the zero-dimensional object array that :func:`numpy.load` returns.

    if not os.path.exists(fname):
        raise FileNotFoundError('No %s at %s' % (what, fname))

    try:
        contents = np.load(fname, allow_pickle=True)
    except UnicodeError:
        warnings.warn('%s could not be read with the default ASCII decoding, which means it was '
                      "written under Python 2; re-reading it with encoding='latin1'." % (fname),
                      stacklevel=2)
        contents = np.load(fname, allow_pickle=True, encoding='latin1')

    if not hasattr(contents, 'item') or contents.shape != ():
        raise ValueError('The %s %s does not contain a single dictionary; got %s'
                         % (what, fname, type(contents).__name__ if not hasattr(contents, 'shape')
                            else 'an array of shape %s' % (contents.shape,)))
    product = contents.item()
    if not isinstance(product, dict):
        raise ValueError('The %s %s contains a %s rather than a dictionary'
                         % (what, fname, type(product).__name__))

    return product


def load_ilc_product(ilc_fname):
    r"""
    Read an ILC product file.

    Parameters
    ----------
    ilc_fname : str
        Path to a ``.npy`` file written by :func:`get_ilc_residuals.run_ilc`.

    Returns
    -------
    dict
        Contents of the product file. The entries used downstream are ``'el'``, ``'cl_residual'`` and, where present, ``'param_dict'``, ``'fsky_val'`` and ``'fg_res_dic'``.

    Raises
    ------
    FileNotFoundError
        If the file does not exist.
    ValueError
        If the file does not hold a single dictionary, which is how ILC products are written.

    Warns
    -----
    UserWarning
        Propagated from ``_load_npy_dict``.

    Notes
    -----
    This loader is deliberately the only place that knows the product layout; :func:`delensing.build_effective_noise` and the routines in :mod:`fisher` take the dictionary itself so that they can equally be handed the return value of :func:`get_ilc_residuals.run_ilc` without anything having been written to disk.
    """

    return _load_npy_dict(ilc_fname, 'ILC product')


def get_fsky(ilc_dic, fsky=None):
    r"""
    Determine the sky fraction to use for one configuration.

    Parameters
    ----------
    ilc_dic : dict
        ILC product, which carries ``'fsky_val'`` when galactic foregrounds were included.
    fsky : float, optional
        Sky fraction supplied by the user, which takes precedence over the stored value.
        Default is ``None``, which uses the stored value.

    Returns
    -------
    float
        Sky fraction to scale the Fisher matrix of this configuration by.

    Raises
    ------
    ValueError
        If no sky fraction is available, which happens when the ILC was run without galactic foregrounds, or if the value is outside ``(0, 1]``.

    Warns
    -----
    UserWarning
        If a supplied value overrides a stored one since the stored value is the sky fraction the foreground spectra were computed on.

    Notes
    -----
    :func:`get_ilc_residuals.build_output_dic` records ``'fsky_val'`` only when ``include_gal`` is set, so a product built without galactic foregrounds carries no sky fraction and one must be supplied.
    The stored value is the fraction remaining after the galactic mask has been intersected with the survey footprint, which is usually the relevant area for the Fisher sum.
    """

    stored = ilc_dic.get('fsky_val')

    if fsky is None:
        if stored is None:
            raise ValueError('The ILC product records no fsky_val, which is the case whenever the ILC '
                             'was run without galactic foregrounds. Supply fsky explicitly.')
        return float(stored)

    fsky = float(fsky)
    if not 0. < fsky <= 1.:
        raise ValueError('fsky must lie in (0, 1], got %s' % (fsky))
    if stored is not None and not np.isclose(fsky, float(stored)):
        warnings.warn('Using the supplied fsky = %g in place of the fsky_val = %g recorded in the ILC '
                      'product, which is the sky fraction its foreground spectra were computed on.'
                      % (fsky, float(stored)), stacklevel=2)

    return fsky


def check_cache_consistency(caches, labels=None):
    r"""
    Refuse forecast caches that do not describe the same theory.

    Parameters
    ----------
    caches : sequence of dict
        Caches from :func:`forecast_cache` or :func:`read_forecast_cache`, one per configuration.
    labels : sequence of str, optional
        Names to report the caches by. Default is ``None``, which numbers them from zero.

    Returns
    -------
    list of str
        The parameters every cache varies, in their common order.

    Raises
    ------
    ValueError
        If no caches are given, if they disagree on the varied parameters or their order, on the fiducial cosmology, on the step sizes of the varied parameters, on ``lmax_calc``, or on the unlensed or lensed spectra.

    Notes
    -----
    Fisher matrices of separate configurations are summed as independent experiments, which is only meaningful if they are derivatives of the same theory about the same fiducial point.
    """

    #Note: The unlensed and lensed spectra are compared exactly rather than approximately. They do not depend on the noise at all, so two caches built from the same cosmology at the same ``lmax_calc`` hold bit-identical copies of them however different their configurations are and any difference means the cosmology, the multipole range or the CLASS version differs. That makes them a stronger check than the recorded metadata, which only catches what somebody remembered to record.

    caches = list(caches)
    if not caches:
        raise ValueError('There are no forecast caches to check.')
    if labels is None:
        labels = ['cache %d' % (index) for index in range(len(caches))]
    labels = list(labels)

    def _settings(cache):
        return cache.get('settings') or {}

    params = list( caches[0].get('vary_params') or _settings(caches[0]).get('vary_params') or [] )
    if not params:
        raise ValueError('%s records no vary_params, so there is nothing to check the others '
                         'against.' % (labels[0]))

    for index in range(1, len(caches)):
        other = list( caches[index].get('vary_params')
                      or _settings(caches[index]).get('vary_params') or [] )
        if other != params:
            raise ValueError('%s varies %s, but %s varies %s. Fisher blocks are summed row by row, '
                             'so every configuration must vary the same parameters in the same '
                             'order.' % (labels[index], ', '.join(other) or 'nothing', labels[0],
                                         ', '.join(params)))

    reference = _settings(caches[0])
    for index in range(1, len(caches)):
        current = _settings(caches[index])
        for name in sorted( reference.get('cosmo_fid', {}) ):
            here = current.get('cosmo_fid', {}).get(name)
            there = reference['cosmo_fid'][name]
            if here is None or float(here) != float(there):
                raise ValueError('%s has %s = %s in its fiducial cosmology where %s has %s. '
                                 'Summing Fisher matrices expanded about different fiducial '
                                 'points is meaningless.'
                                 % (labels[index], name, here, labels[0], there))
        for name in params:
            here = current.get('step_sizes', {}).get(name)
            there = reference.get('step_sizes', {}).get(name)
            if here is not None and there is not None and float(here) != float(there):
                raise ValueError('%s used a step size of %s for %s where %s used %s. The two '
                                 'derivatives are of the same quantity but carry different '
                                 'truncation errors, so they are not the same derivative.'
                                 % (labels[index], here, name, labels[0], there))
        if caches[index].get('lmax_calc') != caches[0].get('lmax_calc'):
            raise ValueError('%s was computed to lmax_calc = %s and %s to %s.'
                             % (labels[index], caches[index].get('lmax_calc'), labels[0],
                                caches[0].get('lmax_calc')))

    for spectrum in ('unlensed', 'lensed'):
        first = ( caches[0].get('powers') or {} ).get(spectrum, {})
        for index in range(1, len(caches)):
            current = ( caches[index].get('powers') or {} ).get(spectrum, {})
            for key in sorted(first):
                if key not in current or not np.array_equal( np.asarray(first[key]),
                                                             np.asarray(current[key]) ):
                    raise ValueError("%s and %s disagree on the %s '%s' spectrum, which does not "
                                     'depend on the noise and so must be identical across '
                                     'configurations of the same cosmology. The two caches were '
                                     'not produced by the same theory calculation.'
                                     % (labels[0], labels[index], spectrum, key))

    return params


# Output

def write_lensing_noise_curves(curves, fname, overwrite=False):
    r"""
    Write a lensing-reconstruction noise product to disk.

    Parameters
    ----------
    curves : dict
        Product assembled by :func:`delensing.lensing_noise_curves`.
    fname : str
        Destination path. Parent directories are created as needed.
    overwrite : bool, optional
        Replace an existing file. Default is ``False``, which refuses rather than discarding a previous result.

    Returns
    -------
    str
        The path written.

    Raises
    ------
    FileExistsError
        If the destination exists and ``overwrite`` is not set.
    ValueError
        If ``curves`` lacks the entries a reader would expect.

    Notes
    -----
    The layout differs from the lensing-noise products under ``products/``, which an earlier pipeline wrote, so a reader of those needs updating rather than reusing.
    The differences are deliberate:

    * multipoles are under ``'el'``, as in the ILC products, rather than ``'els'``;
    * the reconstruction noise is grouped under ``'nl_dd'`` and keyed by estimator, rather than spread over top-level ``Nl_TT``, ``Nl_EB`` and so on;
    * the estimator combining temperature and polarization is ``'TE'``, as CLASS_delens names it, not ``Nl_ET``;
    * arrays are real and finite, rather than ``complex128`` carrying ``nan`` below the first reconstructed multipole;
    * everything is in the deflection convention that CLASS_delens works in and the Fisher sum uses, with :math:`C_\ell^{dd} = \ell(\ell+1) C_\ell^{\phi\phi}`, rather than the convergence convention.

    Convert to convergence with :math:`N_\ell^{\kappa\kappa} = \ell(\ell+1) N_\ell^{dd}/4`; ``'cl_kk'`` is stored alongside so a signal-to-noise ratio needs no conversion at all.
    """

    for key in ('el', 'nl_dd', 'cl_dd'):
        if key not in curves:
            raise ValueError("The product has no '%s' entry; keys are %s"
                             % (key, sorted(curves.keys())))
    if os.path.exists(fname) and not overwrite:
        raise FileExistsError('%s exists already; pass overwrite=True to replace it' % (fname))

    parent = os.path.dirname( os.path.abspath(fname) )
    if parent:
        os.makedirs(parent, exist_ok=True)
    np.save(fname, curves)

    return fname if fname.endswith('.npy') else fname + '.npy'


def convert_legacy_lensing_noise(product, fname=None):
    r"""
    Convert a lensing-noise product from the layout used under ``products/``.

    Parameters
    ----------
    product : dict
        Contents of an older product, holding ``'els'``, ``'cl_kk'`` and the ``Nl_`` entries.
    fname : str, optional
        Path it came from, for the messages. Default is ``None``.

    Returns
    -------
    dict
        The same information in the layout :func:`delensing.lensing_noise_curves` produces, with ``'legacy'`` set to record the conversion.

    Raises
    ------
    ValueError
        If the product holds none of the expected estimators or has no multipole array.

    Warns
    -----
    UserWarning
        If an estimator has an imaginary part that is not zero, which would mean the conversion is discarding information rather than a stored type.

    Notes
    -----
    The noise is stored in the convergence convention and is converted to the deflection convention used since, by

    .. math::

        N_\ell^{dd} = \frac{4}{\ell(\ell+1)}\, N_\ell^{\kappa\kappa}\, ,

    and likewise for the signal, from which ``'cl_phiphi'`` and ``'cl_dd'`` are formed. The stored ``'cl_kk'`` is carried through unchanged. Multipoles below 2 are dropped, both because the conversion is not defined there and because the current products start at 2.

    The arrays are stored as ``complex128`` with the imaginary part zero and are returned as real.

    The older products name the recorded multipole cuts ``lmin``, ``lmax`` and ``lmax_tt``, but they are the cuts on the fields entering the quadratic estimator, whereas ``lmin``, ``lmax`` and ``lmax_TT`` in the current settings are the windows of the Fisher sum. They are renamed here through :data:`LEGACY_LENSING_SETTINGS` rather than carried across since carrying them across would put reconstruction cuts where a later stage reads Fisher windows.

    Reading one of these files emits several ``VisibleDeprecationWarning`` from NumPy about ``align=0``. That comes from unpickling arrays written under an older NumPy and is not suppressed here since suppressing it would also hide the same warning raised for other reasons.
    """

    where = '' if fname is None else ' in %s' % (fname)
    if 'els' not in product:
        raise ValueError('There is no els entry%s, so this is not a lensing-noise product in the '
                         'older layout. Its keys are %s' % (where, sorted(product.keys())))
    present = [key for key in LEGACY_LENSING_ESTIMATORS if key in product]
    if not present:
        raise ValueError('There are no Nl_ entries%s, so there is no reconstruction noise to '
                         'convert. Its keys are %s' % (where, sorted(product.keys())))

    el = np.asarray(product['els'], dtype=float)
    keep = el >= 2.
    el = el[keep]
    to_deflection = 4. / ( el*(el+1.) )

    converted = {'el': el.astype(int) if np.allclose(el, np.round(el)) else el,
                 'legacy': True}

    for key in present:
        values = np.asarray(product[key])
        if np.iscomplexobj(values):
            finite = np.isfinite(values.imag)
            if not np.allclose(values.imag[finite], 0.):
                warnings.warn('%s%s has an imaginary part that is not zero, so taking the real '
                              'part discards information rather than a stored type.'
                              % (key, where), stacklevel=2)
            values = values.real
        converted.setdefault('nl_dd', {})
        converted['nl_dd'][ LEGACY_LENSING_ESTIMATORS[key] ] = \
            np.asarray(values, dtype=float)[keep] * to_deflection

    if 'cl_kk' in product:
        cl_kk = np.asarray(product['cl_kk'], dtype=float)[keep]
        converted['cl_kk'] = cl_kk
        converted['cl_dd'] = cl_kk*to_deflection
        converted['cl_phiphi'] = converted['cl_dd'] / ( el*(el+1.) )

    converted['estimators'] = sorted( converted['nl_dd'].keys() )

    #Recorded so that a converted product has the same shape as a current one. It is not checked
    #against the individual estimators the way lensing_noise_curves checks its own: that guards
    #against a CLASS_delens branch which sets the entry equal to a single estimator, and these
    #files predate CLASS_delens. Theirs was verified to be the exact inverse-variance combination
    #of TT, TE, EE, TB and EB, with MVpol that of EE and EB.
    if 'MV' in converted['nl_dd']:
        converted['mv_estimator'] = 'MV'

    #Renamed rather than carried across: these are reconstruction cuts, not Fisher windows.
    for old, new in LEGACY_LENSING_SETTINGS.items():
        if old in product:
            converted[new] = product[old]
    for key in ('noise_T_uKarmins', 'beam_arcmins'):
        if key in product:
            converted[key] = product[key]

    return converted


def read_lensing_noise_curves(fname):
    r"""
    Read a lensing-reconstruction noise product written by :func:`write_lensing_noise_curves`.

    Parameters
    ----------
    fname : str
        Path to the product.

    Returns
    -------
    dict
        Contents of the product.

    Raises
    ------
    FileNotFoundError
        If the file does not exist.
    ValueError
        If the file does not hold a single dictionary.

    Warns
    -----
    UserWarning
        If the file is in the older layout under ``products/`` and has been converted on read.

    Notes
    -----
    The older products under ``products/`` are recognized by their ``'els'`` entry and converted by :func:`convert_legacy_lensing_noise`. The conversion changes both the names and the convention of what is returned, which is why it is warned about.
    """

    product = _load_npy_dict(fname, 'lensing-noise product')
    if 'els' in product and 'nl_dd' not in product:
        warnings.warn('%s is in the older layout used under products/, and has been converted: the '
                      'reconstruction noise is grouped under nl_dd, renamed as CLASS_delens names '
                      'it, made real and converted from the convergence convention to the '
                      'deflection one. See convert_legacy_lensing_noise.' % (fname), stacklevel=2)
        return convert_legacy_lensing_noise(product, fname)

    return product


def forecast_cache(powers, param_derivs, noise, deflection_noise, reconstruction=None,
                   settings=None, extra=None):
    r"""
    Gather everything needed to rebuild a Fisher matrix without rerunning CLASS.

    Parameters
    ----------
    powers : dict
        Fiducial spectra from :func:`delensing.run_class`.
    param_derivs : dict
        Derivatives from :func:`delensing.parameter_derivatives`.
    noise : dict
        Effective noise from :func:`delensing.build_effective_noise`, in its canonical form with multipoles under ``'el'`` starting at zero.
    deflection_noise : array_like
        Reconstruction noise as it was handed to CLASS, from :func:`delensing.deflection_noise_for_class`.
    reconstruction : dict, optional
        Per-estimator reconstruction noise from :func:`delensing.run_class`. Default is ``None``.
    settings : dict, optional
        Forecasting settings, recorded so the cache is self describing. Default is ``None``.
    extra : dict, optional
        Further entries, such as the ILC product the noise came from. Default is ``None``.

    Returns
    -------
    dict
        Entries ``'powers'``, ``'param_derivs'``, ``'noise'``, ``'deflection_noise'`` and, where given, ``'reconstruction'``, ``'settings'`` and whatever ``extra`` holds.

    Notes
    -----
    The expensive part of a forecast is the CLASS runs and nothing in the Fisher stage needs to repeat them: changing the multipole windows, the sky fraction, the priors or which parameters are held fixed all act on quantities this cache already holds.
    Keeping it therefore turns a scan over any of those into seconds rather than potentially hours.

    The noise is stored in the canonical form rather than either of the re-indexed views since :func:`delensing.noise_for_fisher` and :func:`delensing.noise_for_class` derive those cheaply and storing one of them would invite using it in the wrong place.

    Written and read with :func:`write_forecast_cache` and :func:`read_forecast_cache`. These files are large since they hold every derivative at every multipole, so they belong under ``results/`` rather than in the repository.
    """

    cache = {'powers': powers,
             'param_derivs': param_derivs,
             'noise': noise,
             'deflection_noise': np.asarray(deflection_noise)}
    if reconstruction is not None:
        cache['reconstruction'] = reconstruction
    if settings is not None:
        cache['settings'] = settings
    if extra is not None:
        cache.update(extra)

    return cache


def write_forecast_cache(cache, fname, overwrite=False):
    r"""
    Write a forecast cache to disk.

    Parameters
    ----------
    cache : dict
        Assembled by :func:`forecast_cache`.
    fname : str
        Destination path. Parent directories are created as needed.
    overwrite : bool, optional
        Replace an existing file. Default is ``False``.

    Returns
    -------
    str
        The path written.

    Raises
    ------
    FileExistsError
        If the destination exists and ``overwrite`` is not set.
    ValueError
        If the cache lacks an entry the Fisher stage needs.
    """

    for key in ('powers', 'param_derivs', 'noise', 'deflection_noise'):
        if key not in cache:
            raise ValueError("The cache has no '%s' entry; keys are %s" % (key, sorted(cache.keys())))
    if os.path.exists(fname) and not overwrite:
        raise FileExistsError('%s exists already; pass overwrite=True to replace it' % (fname))

    parent = os.path.dirname( os.path.abspath(fname) )
    if parent:
        os.makedirs(parent, exist_ok=True)
    np.save(fname, cache)

    return fname if fname.endswith('.npy') else fname + '.npy'


def read_forecast_cache(fname):
    r"""
    Read a forecast cache written by :func:`write_forecast_cache`.

    Parameters
    ----------
    fname : str
        Path to the cache.

    Returns
    -------
    dict
        Contents of the cache.

    Raises
    ------
    FileNotFoundError
        If the file does not exist.
    ValueError
        If the file does not hold a single dictionary.
    """

    return _load_npy_dict(fname, 'forecast cache')


def forecast_product(total, configurations, settings=None, windows=None, biases=None, extra=None):
    r"""
    Assemble the forecast for output.

    Parameters
    ----------
    total : dict
        Combined Fisher matrix and uncertainties from :func:`fisher.combine_fisher`.
    configurations : sequence of dict
        One entry per configuration, recording at least its ``'ilc_fname'``, ``'fsky'`` and ``'label'``, and where written its ``'cache_fname'`` and ``'lensing_noise_fname'``.
    settings : dict, optional
        Settings the forecast was run under, from :func:`load_fisher_params`. Default is ``None``.
    windows : dict, optional
        Multipole windows from :func:`fisher.fisher_windows`. Default is ``None``.
    biases : dict, optional
        Foreground-induced parameter biases from :func:`fisher.parameter_bias`. Default is ``None``.
    extra : dict, optional
        Further entries to record. Default is ``None``.

    Returns
    -------
    dict
        Every entry of ``total`` plus ``'configurations'`` and, where given, ``'settings'``, ``'windows'``, ``'bias'``, ``'bias_in_sigma'`` and whatever ``extra`` holds.

    Raises
    ------
    ValueError
        If ``total`` is not a combined Fisher matrix or no configurations are given.
    """

    #Note: The uncertainties are kept for every spectrum type rather than for the delensed spectra alone. They cost nothing once the derivatives exist, and the lensed and unlensed values are what the delensed number is quoted against.

    for key in ('fisher', 'covariance', 'sigmas', 'params'):
        if key not in total:
            raise ValueError("The total has no '%s' entry; keys are %s. It comes from "
                             'fisher.combine_fisher.' % (key, sorted(total.keys())))
    configurations = list(configurations)
    if not configurations:
        raise ValueError('A forecast product records the configurations it was built from and '
                         'none were given.')

    product = dict(total)
    product['configurations'] = configurations
    if settings is not None:
        product['settings'] = settings
    if windows is not None:
        product['windows'] = windows
    if biases is not None:
        product['bias'] = biases.get('bias')
        product['bias_in_sigma'] = biases.get('bias_in_sigma')
        product['bias_labels'] = biases.get('labels')
    if extra is not None:
        product.update(extra)

    return product


def write_forecast_product(product, fname, overwrite=False):
    r"""
    Write a forecast product to disk.

    Parameters
    ----------
    product : dict
        Assembled by :func:`forecast_product`.
    fname : str
        Destination path. Parent directories are created as needed.
    overwrite : bool, optional
        Replace an existing file. Default is ``False``.

    Returns
    -------
    str
        The path written.

    Raises
    ------
    FileExistsError
        If the destination exists and ``overwrite`` is not set.
    ValueError
        If the product lacks an entry that makes it a forecast.
    """

    for key in ('fisher', 'covariance', 'sigmas', 'params', 'configurations'):
        if key not in product:
            raise ValueError("The forecast has no '%s' entry; keys are %s."
                             % (key, sorted(product.keys())))
    if os.path.exists(fname) and not overwrite:
        raise FileExistsError('%s exists already; pass overwrite=True to replace it.' % (fname))

    parent = os.path.dirname( os.path.abspath(fname) )
    if parent:
        os.makedirs(parent, exist_ok=True)
    np.save(fname, product)

    return fname if fname.endswith('.npy') else fname + '.npy'


def read_forecast_product(fname):
    r"""
    Read a forecast product written by :func:`write_forecast_product`.

    Parameters
    ----------
    fname : str
        Path to the product.

    Returns
    -------
    dict
        Contents of the product.

    Raises
    ------
    FileNotFoundError
        If the file does not exist.
    ValueError
        If the file does not hold a single dictionary.
    """

    return _load_npy_dict(fname, 'forecast')


def array_fingerprint(*arrays):
    r"""
    Reduce one or more arrays to a short string that changes when their contents do.

    Parameters
    ----------
    *arrays : array_like or dict or None
        Arrays to fingerprint. A dictionary is folded in by sorted key, so the result does not depend on insertion order. ``None`` contributes nothing.

    Returns
    -------
    str or None
        Sixteen hexadecimal characters, or ``None`` if nothing was supplied.

    Notes
    -----
    Used to record which noise a set of non-Gaussian derivative matrices was computed with so that a stored set can be refused when it does not match the one being asked for.
    Everything is coerced to ``float64`` and to C order before hashing so that a view, a byte order or a stride cannot change the answer for the same numbers.
    This identifies contents, not provenance; it is not a checksum against corruption and makes no cryptographic claim.
    """

    digest = hashlib.sha256()
    supplied = False
    for array in arrays:
        if array is None:
            continue
        supplied = True
        if isinstance(array, dict):
            for key in sorted(array):
                digest.update( str(key).encode() )
                digest.update( np.ascontiguousarray(array[key], dtype=float).tobytes() )
        else:
            digest.update( np.ascontiguousarray(array, dtype=float).tobytes() )

    return digest.hexdigest()[:16] if supplied else None


def ng_derivative_product(lensing_derivs, unlensed_derivs, cosmo_fid, lmax_calc, derivative_type,
                          cmb_noise=None, deflection_noise=None, class_parameters=None,
                          extra=None):
    r"""
    Gather a set of non-Gaussian derivative matrices with the record of how it was made.

    Parameters
    ----------
    lensing_derivs : dict
        ``dCldCLd`` from :func:`delensing.lensing_derivatives`.
    unlensed_derivs : dict or None
        ``dCldCLu``, or ``None`` when it was not requested.
    cosmo_fid : dict
        Cosmology the matrices were computed at.
    lmax_calc : int
        Highest multipole they cover.
    derivative_type : str
        One of :data:`delensing.DERIVATIVE_TYPES`.
    cmb_noise : dict, optional
        Effective noise the run was given, recorded as a fingerprint rather than in full. Default is ``None``.
    deflection_noise : array_like, optional
        Reconstruction noise the run was given, likewise. Default is ``None``.
    class_parameters : dict, optional
        The deck CLASS recorded, from :func:`delensing.class_run_parameters`. Default is ``None``.
    extra : dict, optional
        Further entries for the manifest. Default is ``None``.

    Returns
    -------
    dict
        Entries ``'lensing_derivs'``, ``'unlensed_derivs'`` where present, and ``'manifest'``.

    Notes
    -----
    The manifest records the DRAFT-side inputs rather than the CLASS deck because the deck names the cosmology as CLASS names it and recovering a DRAFT cosmology from it would mean reproducing the translation inside ``classWrapTools``. The deck is stored alongside as provenance when it is supplied, but is not what :func:`check_ng_derivatives` compares.

    The noise enters only as a fingerprint. Storing it in full would add little and invite reading it back as though it were the canonical copy, which belongs in the forecast cache.

    The lensed matrices do not depend on the noise, so one product serves every configuration; the delensed ones do, so a product is specific to the configuration named in its manifest.
    """

    manifest = {'cosmo_fid': dict(cosmo_fid),
                'lmax_calc': int(lmax_calc),
                'derivative_type': str(derivative_type),
                'include_unlensed': unlensed_derivs is not None,
                'noise_fingerprint': array_fingerprint(cmb_noise),
                'deflection_fingerprint': array_fingerprint(deflection_noise)}
    if class_parameters is not None:
        manifest['class_parameters'] = dict(class_parameters)
    if extra is not None:
        manifest.update(extra)

    product = {'lensing_derivs': lensing_derivs, 'manifest': manifest}
    if unlensed_derivs is not None:
        product['unlensed_derivs'] = unlensed_derivs

    return product


def write_ng_derivatives(product, fname, overwrite=False):
    r"""
    Write a set of non-Gaussian derivative matrices to disk.

    Parameters
    ----------
    product : dict
        Assembled by :func:`ng_derivative_product`.
    fname : str
        Destination path. Parent directories are created as needed.
    overwrite : bool, optional
        Replace an existing file. Default is ``False``.

    Returns
    -------
    str
        The path written.

    Raises
    ------
    FileExistsError
        If the destination exists and ``overwrite`` is not set.
    ValueError
        If the product lacks an entry or holds a matrix that is not square in the multipole.

    Notes
    -----
    These files are the largest DRAFT writes, some 2.4 GB per derivative type at lmax = 5000 and quadratic in it.
    The matrices are kept at ``float64``, which is what CLASS writes and what they are held at in memory.
    """

    for key in ('lensing_derivs', 'manifest'):
        if key not in product:
            raise ValueError("The product has no '%s' entry; keys are %s"
                             % (key, sorted(product.keys())))
    lmax_calc = int(product['manifest']['lmax_calc'])
    for group in ('lensing_derivs', 'unlensed_derivs'):
        for key, value in product.get(group, {}).items():
            shape = np.asarray(value).shape
            if shape != (lmax_calc+1, lmax_calc+1):
                raise ValueError('%s[%r] has shape %s, expected (%d, %d) for a set reaching '
                                 'multipole %d'
                                 % (group, key, shape, lmax_calc+1, lmax_calc+1, lmax_calc))
    if os.path.exists(fname) and not overwrite:
        raise FileExistsError('%s exists already; pass overwrite=True to replace it' % (fname))

    parent = os.path.dirname( os.path.abspath(fname) )
    if parent:
        os.makedirs(parent, exist_ok=True)
    np.save(fname, product)

    return fname if fname.endswith('.npy') else fname + '.npy'


def read_ng_derivatives(fname):
    r"""
    Read a set of non-Gaussian derivative matrices written by :func:`write_ng_derivatives`.

    Parameters
    ----------
    fname : str
        Path to the product.

    Returns
    -------
    dict
        Contents of the product.

    Raises
    ------
    FileNotFoundError
        If the file does not exist.
    ValueError
        If the file does not hold a product with a manifest.
    """

    product = _load_npy_dict(fname, 'non-Gaussian derivative product')
    if 'manifest' not in product or 'lensing_derivs' not in product:
        raise ValueError('%s does not hold a set of derivative matrices written by '
                         'write_ng_derivatives; its keys are %s'
                         % (fname, sorted(product.keys())))

    return product


def check_ng_derivatives(product, cosmo_fid, lmax_calc, derivative_type, include_unlensed=True,
                         cmb_noise=None, deflection_noise=None):
    r"""
    Report why a stored set of derivative matrices does not match the one being asked for.

    Parameters
    ----------
    product : dict
        Read by :func:`read_ng_derivatives`.
    cosmo_fid : dict
        Cosmology being asked for.
    lmax_calc : int
        Highest multipole being asked for.
    derivative_type : str
        One of :data:`delensing.DERIVATIVE_TYPES`.
    include_unlensed : bool, optional
        Whether the derivatives with respect to the unlensed spectra are needed. Default is ``True``.
    cmb_noise : dict, optional
        Effective noise being asked for. Default is ``None``, which does not compare it.
    deflection_noise : array_like, optional
        Reconstruction noise being asked for. Default is ``None``.

    Returns
    -------
    list of str
        Complaints, one per mismatch, empty when the stored set can be used as it is.

    Notes
    -----
    A stored set is refused rather than silently reused because nothing about the matrices themselves reveals the cosmology or the noise they were computed at and CLASS leaves its own output in place indefinitely.
    The noise is compared only when it is supplied since the lensed matrices do not depend on it and a lensed product is deliberately shared across configurations.
    """

    manifest = product.get('manifest', {})
    complaints = []

    if int( manifest.get('lmax_calc', -1) ) != int(lmax_calc):
        complaints.append('it reaches multipole %s, not %d'
                          % (manifest.get('lmax_calc'), lmax_calc))
    if str( manifest.get('derivative_type') ) != str(derivative_type):
        complaints.append('it holds %s derivatives, not %s'
                          % (manifest.get('derivative_type'), derivative_type))
    if include_unlensed and 'unlensed_derivs' not in product:
        complaints.append('it has no derivatives with respect to the unlensed spectra')

    stored_cosmo = manifest.get('cosmo_fid', {})
    for name in sorted( set(stored_cosmo) | set(cosmo_fid) ):
        if name not in stored_cosmo:
            complaints.append('it was computed without %s' % (name))
        elif name not in cosmo_fid:
            complaints.append('it was computed with %s, which is not being asked for' % (name))
        elif not np.isclose(float(stored_cosmo[name]), float(cosmo_fid[name]), rtol=1e-12, atol=0.):
            complaints.append('it was computed at %s = %s, not %s'
                              % (name, stored_cosmo[name], cosmo_fid[name]))

    for label, stored, asked in (('noise', manifest.get('noise_fingerprint'),
                                  array_fingerprint(cmb_noise)),
                                 ('deflection noise', manifest.get('deflection_fingerprint'),
                                  array_fingerprint(deflection_noise))):
        if asked is not None and stored != asked:
            complaints.append('it was computed with a different %s (%s, not %s)'
                              % (label, stored, asked))

    return complaints


# Driver

def obtain_ng_derivatives(mode,
                          cosmo_fid,
                          lmax_calc,
                          derivative_type,
                          fname=None,
                          cmb_noise=None,
                          deflection_noise=None,
                          recon_mask=None,
                          include_unlensed=True,
                          root_name=None,
                          paths=None,
                          overwrite=False,
                          verbose=True
                          ):
    r"""
    Obtain the non-Gaussian derivative matrices by the route the parameter file asks for.

    Parameters
    ----------
    mode : str
        One of :data:`delensing.NG_DERIVATIVE_MODES`, as ``ng_derivatives`` in the parameter file.
    cosmo_fid : dict
        Cosmology to compute at or to check a stored set against.
    lmax_calc : int
        Highest multipole.
    derivative_type : str
        One of :data:`delensing.DERIVATIVE_TYPES`.
    fname : str, optional
        Product to write or read. Required for ``'save'`` and ``'load'``, unused for ``'recompute'``. Default is ``None``.
    cmb_noise : dict, optional
        Effective noise, needed for the delensed derivatives. Default is ``None``.
    deflection_noise : array_like, optional
        Reconstruction noise, likewise. Default is ``None``.
    recon_mask : dict, optional
        Reconstruction cuts from :func:`delensing.reconstruction_mask`. Default is ``None``.
    include_unlensed : bool, optional
        Whether the derivatives with respect to the unlensed spectra are wanted. Default is ``True``.
    root_name : str, optional
        Base name for the CLASS run. Default is ``None``, which generates a unique one.
    paths : dict, optional
        Resolved paths from :func:`delensing.resolve_paths`. Default is ``None``.
    overwrite : bool, optional
        Replace an existing product when saving. Default is ``False``.
    verbose : bool, optional
        Whether to report what is being done. Default is ``True``.

    Returns
    -------
    lensing_derivs : dict
        ``dCldCLd``.
    unlensed_derivs : dict or None
        ``dCldCLu``, or ``None`` when it was not requested.

    Raises
    ------
    ValueError
        If ``mode`` is not recognized, if a path is needed and was not given, or if a stored set does not match what is being asked for.

    Notes
    -----
    The three modes trade a large file against a long calculation, which is why the parameter file has no default for them: ``'recompute'`` runs CLASS and keeps nothing, ``'save'`` runs CLASS and writes the matrices for reuse, ``'load'`` reads them back without running CLASS at all.

    A stored set is checked against the cosmology, the multipole range, the derivative type and the noise before it is used, and refused with the reason rather than reused silently. That check is the whole point of the manifest: CLASS leaves its output in place indefinitely and the matrices themselves carry no record of what produced them.
    """

    mode = str(mode)
    if mode not in delensing.NG_DERIVATIVE_MODES:
        raise ValueError('ng_derivatives must be one of %s, got %r'
                         % (', '.join(delensing.NG_DERIVATIVE_MODES), mode))
    if mode in ('save', 'load') and fname is None:
        raise ValueError("ng_derivatives is %r, which needs a path to %s, but none was given."
                         % (mode, 'write to' if mode == 'save' else 'read from'))

    if mode == 'load':
        if verbose:
            print('Reading %s derivatives from %s' % (derivative_type, fname))
        product = read_ng_derivatives(fname)
        complaints = check_ng_derivatives(product, cosmo_fid, lmax_calc, derivative_type,
                                          include_unlensed=include_unlensed, cmb_noise=cmb_noise,
                                          deflection_noise=deflection_noise)
        if complaints:
            raise ValueError('%s does not hold the derivatives being asked for: %s. Recompute '
                             'them or point ng_derivatives at the right product.'
                             % (fname, '; '.join(complaints)))

        return product['lensing_derivs'], product.get('unlensed_derivs')

    if root_name is None:
        root_name = 'draft_derivs_%s' % (derivative_type)
    lensing_derivs, unlensed_derivs = delensing.lensing_derivatives(
        cosmo_fid, lmax_calc, derivative_type=derivative_type, cmb_noise=cmb_noise,
        deflection_noise=deflection_noise, recon_mask=recon_mask,
        include_unlensed=include_unlensed, root_name=root_name, paths=paths, verbose=verbose)

    if mode == 'save':
        if paths is None:
            paths = delensing.resolve_paths()
        try:
            class_parameters = delensing.class_run_parameters(paths['class_data_dir'], root_name)
        except FileNotFoundError:
            class_parameters = None
        product = ng_derivative_product(lensing_derivs, unlensed_derivs, cosmo_fid, lmax_calc,
                                        derivative_type, cmb_noise=cmb_noise,
                                        deflection_noise=deflection_noise,
                                        class_parameters=class_parameters,
                                        extra={'root_name': root_name})
        written = write_ng_derivatives(product, fname, overwrite=overwrite)
        if verbose:
            print('  wrote %s (%.1f GB)'
                  % (written, os.path.getsize(written)/1e9))

    return lensing_derivs, unlensed_derivs


def _per_configuration(value, count, what, sentinel=None):
    """
    Spread a per-configuration argument over the configurations.

    Parameters
    ----------
    value : object
        ``None``, a single value when there is one configuration or a sequence of ``count`` values.
    count : int
        Number of configurations.
    what : str
        Name of the argument, for the error message.
    sentinel : str, optional
        String standing for "not given here", replaced by ``None``. Default is ``None``.

    Returns
    -------
    list
        One entry per configuration.

    Raises
    ------
    ValueError
        If a sequence of the wrong length is given or a single value for several configurations.
    """

    if value is None:
        return [None] * count
    if isinstance(value, (str, float, int)) or not hasattr(value, '__len__'):
        if count != 1:
            raise ValueError('%s was given once but there are %d configurations. It is a property '
                             'of each one, so give %d values in the same order as the ILC '
                             'products.' % (what, count, count))
        entries = [value]
    else:
        entries = list(value)
    if len(entries) != count:
        raise ValueError('%s has %d values but there are %d configurations. Give one per ILC '
                         'product, in the same order.' % (what, len(entries), count))

    return [None if (sentinel is not None and isinstance(entry, str) and entry == sentinel)
            else entry for entry in entries]


def report_forecast(product, spectrum='delensed'):
    r"""
    Print the projected uncertainties of a forecast.

    Parameters
    ----------
    product : dict
        Forecast from :func:`forecast_product` or :func:`read_forecast_product`.
    spectrum : str, optional
        Spectrum type to report. Default is ``'delensed'``; ``'lensed'`` and ``'unlensed'`` also available.

    Raises
    ------
    ValueError
        If the product holds no uncertainties for ``spectrum``.
    """

    if spectrum not in product.get('sigmas', {}):
        raise ValueError("The forecast holds no %r uncertainties; it has %s"
                         % (spectrum, ', '.join( sorted( product.get('sigmas', {}) ) ) or 'none'))
    params = list(product['params'])

    print('\nProjected 1 sigma uncertainties, %s spectra' % (spectrum))
    for position, record in enumerate(product['configurations']):
        print('  %-14s %s at fsky %.4f%s'
              % ('configurations' if position == 0 else '', record['label'], record['fsky'],
                 '   (from a cache)' if record.get('from_cache') else ''))
    print('  %-14s %s   priors %s   fixed %s'
          % ('total', product.get('satellite_label') or 'no low-ell block',
             product.get('priors') or 'none', ', '.join( product.get('params_fixed') or [] ) or 'none'))
    print()
    print('  %-11s %s' % ('', ' '.join( '%12s' % (name) for name in params )))
    print('  %-11s %s' % ('sigma', ' '.join( '%12.5g' % (product['sigmas'][spectrum][name])
                                             for name in params )))
    if product.get('bias') is not None and spectrum in product['bias']:
        print('  %-11s %s' % ('bias', ' '.join( '%12.5g' % (product['bias'][spectrum][name])
                                                for name in params )))
        print('  %-11s %s' % ('bias/sigma',
                              ' '.join( '%12.4f' % (product['bias_in_sigma'][spectrum][name])
                                        for name in params )))


def run_forecast(ilc,
                 fsky=None,
                 cache=None,
                 label=None,
                 opfname=None,
                 paramfile=PARAMFILE,
                 settings=None,
                 fix_params=None,
                 write_cache=False,
                 write_lensing_noise=False,
                 save=True,
                 overwrite=False,
                 paths=None,
                 verbose=True
                 ):
    r"""
    Forecast one or more survey configurations and combine them into projected uncertainties.

    Parameters
    ----------
    ilc : str or dict or sequence
        ILC product to forecast from, either a path to a ``.npy`` product or the dictionary :func:`get_ilc_residuals.run_ilc` returns, or a sequence of either.
        Several are treated as independent experiments whose Fisher matrices are summed.
    fsky : float or str or sequence, optional
        Sky fraction of each configuration, one per entry of ``ilc``, or ``None`` or :data:`FSKY_AUTO` to use the ``fsky_val`` recorded in that product. Default is ``None``.
    cache : str or sequence, optional
        Forecast cache to read for each configuration in place of running CLASS, one per entry of ``ilc``, with ``None`` or :data:`CACHE_NONE` to run CLASS for that one. Default is ``None``.
    label : str, optional
        Name of this forecast. Default is ``None``, which uses the product name when there is one configuration and is required when there are several.
    opfname : str, optional
        Path of the forecast product. Default is ``None``, which builds one under :data:`FORECAST_DIR`.
    paramfile : str, optional
        Parameter file to read. Default is :data:`PARAMFILE`. Ignored when ``settings`` is given.
    settings : dict, optional
        Settings from :func:`load_fisher_params`, in place of reading a parameter file. Default is ``None``.
    fix_params : sequence of str, optional
        Parameters to hold fixed after the total is formed, passed to :func:`fisher.combine_fisher`. Default is ``None``.
    write_cache : bool, optional
        Write the forecast cache of each configuration. Default is ``False``.
    write_lensing_noise : bool, optional
        Write the lensing-noise curves of each configuration. Default is ``False``.
    save : bool, optional
        Write the forecast product. Default is ``True``.
    overwrite : bool, optional
        Replace existing output files. Default is ``False``.
    paths : dict, optional
        Resolved paths from :func:`delensing.resolve_paths`. Default is ``None``, which resolves them and, when any configuration is to be run rather than read from a cache, checks the setup with :func:`delensing.check_setup` first.
    verbose : bool, optional
        Report each step as it runs. Default is ``True``.

    Returns
    -------
    product : dict
        Contents of the forecast product, from :func:`forecast_product`.
    opfname : str
        Path the result was or would have been written to.

    Raises
    ------
    ValueError
        If no configurations are given, if a per-configuration argument does not match them in number, or if several configurations are given without a label or propagated from the routines below.
    RuntimeError
        Propagated from :func:`delensing.check_setup` when a configuration is to be run and FisherLens or a compiled CLASS_delens is not usable.

    Warns
    -----
    UserWarning
        If a cache was built from a different ILC product than the one supplied alongside it, and propagated from the routines below.

    Notes
    -----
    Each configuration costs one iterative CLASS_delens run plus two per varied parameter because the delensed spectra and their derivatives depend on the noise of that configuration.
    The cache of each is therefore written as soon as it is complete rather than at the end so that a failure part way through a multi-configuration run does not discard the configurations that already succeeded.

    The satellite block is built from the first configuration's cache and added once.
    It may be built from any of them: it spans multipoles below ``lmin_ilc``, where the effective noise holds no ILC residual and is the satellite's own, and it uses the lensed spectra, which do not depend on the noise. Both were verified by perturbing each in turn and finding the block unchanged.
    """

    ilc_list = [ilc] if isinstance(ilc, (str, dict)) else list(ilc)
    if not ilc_list:
        raise ValueError('There is nothing to forecast. Give at least one ILC product, as a path '
                         'or as the dictionary get_ilc_residuals.run_ilc returns.')
    count = len(ilc_list)
    fsky_list = _per_configuration(fsky, count, 'fsky', sentinel=FSKY_AUTO)
    cache_list = _per_configuration(cache, count, 'cache', sentinel=CACHE_NONE)

    if settings is None:
        settings = load_fisher_params(paramfile)
    lmax = int(settings['lmax'])
    lmax_calc = lmax + int(settings['lbuffer'])
    lmin_ilc = int(settings['lmin_ilc'])
    vary_params = list(settings['vary_params'])

    names = []
    for index, entry in enumerate(ilc_list):
        if isinstance(entry, str):
            base = os.path.basename(entry)
            names.append(base[:-4] if base.endswith('.npy') else base)
        else:
            #a dictionary handed straight from run_ilc carries no file name of its own, so it is named by position.
            names.append('configuration_%d' % (index))
    if label is None:
        if count > 1:
            raise ValueError('A forecast combining %d configurations needs a name of its own, '
                             'since it belongs to no single ILC product and has nowhere to be '
                             'written. Give label or -label on the command line.' % (count))
        label = names[0]
    opdir = os.path.join(FORECAST_DIR, str(label))
    if opfname is None:
        opfname = os.path.join(opdir, '%s_forecast.npy' % (label))

    if verbose:
        print('\nForecast %s: %d configuration(s), lmax %d plus a buffer of %d'
              % (label, count, lmax, int(settings['lbuffer'])))
        print('  varying %s' % (', '.join(vary_params)))

    if any(item is None for item in cache_list):
        #a fresh checkout reports every missing piece of the CLASS_delens setup at once here, rather than one per attempt from inside the first run.
        paths = delensing.check_setup(verbose=verbose, **(paths or {}))

    caches, configurations = [], []
    for index, entry in enumerate(ilc_list):
        name = names[index]
        if verbose:
            print('\n[%d/%d] %s' % (index+1, count, name))
        ilc_dic = entry if isinstance(entry, dict) else load_ilc_product(entry)
        record = {'label': name,
                  'ilc_fname': entry if isinstance(entry, str) else None,
                  'fsky': get_fsky(ilc_dic, fsky_list[index]),
                  'cache_fname': None,
                  'lensing_noise_fname': None,
                  'from_cache': cache_list[index] is not None}
        if verbose:
            print('  fsky %.4f%s' % (record['fsky'],
                                     '' if fsky_list[index] is None else '  (supplied)'))

        if cache_list[index] is not None:
            this_cache = read_forecast_cache(cache_list[index])
            record['cache_fname'] = cache_list[index]
            built_from = this_cache.get('ilc_fname')
            if built_from is not None and record['ilc_fname'] is not None \
                    and os.path.basename(record['ilc_fname']) != built_from:
                warnings.warn('The cache %s was built from %s but is being used alongside %s. The '
                              'spectra, derivatives and noise would then come from one '
                              'configuration while the sky fraction and any foreground residual '
                              'come from another.'
                              % (os.path.basename(cache_list[index]), built_from,
                                 os.path.basename(record['ilc_fname'])), stacklevel=2)
            if verbose:
                print('  read %s, no CLASS run' % (os.path.basename(cache_list[index])))
        else:
            noise = delensing.build_effective_noise(
                ilc_dic, lmax_calc, lmin_ilc=lmin_ilc,
                planck_expname=str(settings.get('sat_expname', 'planck')),
                verbose=verbose)
            recon_mask = delensing.reconstruction_mask(int(settings['lmin_lensing']),
                                                       int(settings['lmax_T_lensing']),
                                                       int(settings['lmax_P_lensing']))
            root_name = '%s_c%d' % (label, index)
            powers, reconstruction = delensing.run_class(
                settings['cosmo_fid'], lmax_calc, cmb_noise=noise, recon_mask=recon_mask,
                delensing=settings.get('delensing'), root_name=root_name, paths=paths,
                verbose=verbose)
            deflection = delensing.deflection_noise_for_class(reconstruction)
            param_derivs = delensing.parameter_derivatives(
                settings['cosmo_fid'], settings['step_sizes'], vary_params, lmax_calc, noise,
                deflection, root_name='%s_deriv' % (root_name), paths=paths, verbose=verbose)
            this_cache = forecast_cache(
                powers, param_derivs, noise, deflection, reconstruction=reconstruction,
                settings=settings,
                extra={'ilc_fname': os.path.basename(entry) if isinstance(entry, str) else name,
                       'lmax_calc': lmax_calc, 'vary_params': vary_params})

        #each cache is written as it is finished rather than at the end so that a failure on a later configuration does not discard the CLASS runs already paid for.
        if write_cache:
            record['cache_fname'] = write_forecast_cache(
                this_cache, os.path.join(opdir, '%s_cache.npy' % (name)), overwrite=overwrite)
            if verbose:
                print('  wrote %s' % (record['cache_fname']))
        if write_lensing_noise:
            if this_cache.get('reconstruction') is None:
                warnings.warn('The lensing-noise curves of %s were asked for but its cache holds '
                              'no reconstruction, which is the case whenever the deflection noise '
                              'was supplied rather than reconstructed. Nothing is written for it.'
                              % (name), stacklevel=2)
            else:
                curves = delensing.lensing_noise_curves(
                    this_cache['powers'], this_cache['reconstruction'],
                    lmax=int(settings['lmax_dd']), settings=settings,
                    extra={'ilc_fname': this_cache.get('ilc_fname'), 'fsky': record['fsky']})
                record['lensing_noise_fname'] = write_lensing_noise_curves(
                    curves, os.path.join(opdir, '%s_lensing_noise.npy' % (name)),
                    overwrite=overwrite)
                if verbose:
                    print('  wrote %s' % (record['lensing_noise_fname']))

        caches.append(this_cache)
        configurations.append(record)

    check_cache_consistency(caches, labels=names)
    windows = fisher.fisher_windows(settings)

    blocks, vectors = [], []
    for index, this_cache in enumerate(caches):
        lensing_derivs, unlensed_derivs = None, None
        if int( settings.get('do_non_gaussian', 0) ):
            lensing_derivs, unlensed_derivs = obtain_ng_derivatives(
                settings['ng_derivatives'], settings['cosmo_fid'], lmax_calc, 'delensed',
                fname=os.path.join(opdir, '%s_ng_delensed.npy' % (names[index])),
                cmb_noise=this_cache['noise'], deflection_noise=this_cache['deflection_noise'],
                recon_mask=delensing.reconstruction_mask(int(settings['lmin_lensing']),
                                                         int(settings['lmax_T_lensing']),
                                                         int(settings['lmax_P_lensing'])),
                paths=paths, overwrite=overwrite, verbose=verbose)
        blocks.append( fisher.survey_fisher_block(
            this_cache['powers'], this_cache['param_derivs'], this_cache['noise'],
            this_cache['deflection_noise'], vary_params, fsky=configurations[index]['fsky'],
            windows=windows, lensing_derivs=lensing_derivs, unlensed_derivs=unlensed_derivs,
            label=names[index]) )

    satellite = None
    if int( settings.get('include_sat', 0) ):
        #built from the first cache and added once. Any of them would do: the block spans multipoles below lmin_ilc, where the effective noise carries no ILC residual, and it uses the lensed spectra, which do not depend on the noise.
        satellite = fisher.satellite_fisher_block(
            caches[0]['powers'], caches[0]['param_derivs'], caches[0]['noise'], vary_params,
            fsky=float(settings['fsky_sat']), lmin_sat=int(settings['lmin_sat']),
            lmax_sat=int(settings['lmax_sat']), label=str( settings.get('sat_expname', 'planck') ))

    total = fisher.combine_fisher(blocks, satellite=satellite, priors=settings.get('priors'),
                                  fix_params=fix_params)

    biases = None
    if int( settings.get('do_foreground_bias', 0) ):
        fractions = settings.get('fg_bias_fractions') or {}
        for index, entry in enumerate(ilc_list):
            ilc_dic = entry if isinstance(entry, dict) else load_ilc_product(entry)
            systematic = fisher.foreground_systematic(ilc_dic, fractions, lmax=lmax)
            vectors.append( fisher.bias_vector(
                caches[index]['powers'], caches[index]['param_derivs'], caches[index]['noise'],
                caches[index]['deflection_noise'], vary_params, configurations[index]['fsky'],
                systematic, windows=windows, label=names[index]) )
        biases = fisher.parameter_bias(total, vectors)

    product = forecast_product(total, configurations, settings=settings, windows=windows,
                               biases=biases, extra={'label': label, 'fg_bias_fractions':
                                                     settings.get('fg_bias_fractions') or {}})
    if save:
        opfname = write_forecast_product(product, opfname, overwrite=overwrite)

    if verbose:
        report_forecast(product)
        print('\n%s' % (opfname if save else '%s (not written; save is off)' % (opfname)))

    return product, opfname


def main(argv=None):
    """
    Parse the command line and run the forecast.

    Parameters
    ----------
    argv : list of str, optional
        Argument list. Default is ``sys.argv[1:]``.
    """

    args = parse_args(argv)
    run_forecast(args.ilc_fname, fsky=args.fsky, cache=args.cache, label=args.label,
                 opfname=args.opfname, paramfile=args.paramfile, write_cache=bool(args.write_cache),
                 write_lensing_noise=bool(args.write_lensing_noise), save=bool(args.save),
                 overwrite=bool(args.overwrite), verbose=bool(args.verbose))
    print('\nDone.')


if __name__ == '__main__':
    main()
