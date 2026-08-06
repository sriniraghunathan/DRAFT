r"""
Compute delensed CMB spectra and Fisher forecasts from ILC residuals.

This is the driver for the forecasting stage of DRAFT.
It reads the residual power spectra produced by the component separation in ``get_ilc_residuals.py``, combines them with the effective satellite noise, drives the iterative delensing of CLASS_delens through FisherLens, and reduces the result to a Fisher matrix and projected parameter uncertainties.

The command line is added in a later step. The routines below are usable directly::

    import get_fisher_forecasts

    ilc_dic = get_fisher_forecasts.load_ilc_product('products/202310xx_PBDR_config/s4wide/s4wide_ilc_galaxy1_27-39-93-145-225-278_TT-EE_galmask2_AZ_lmax6500_for7years.npy')
    fsky = get_fisher_forecasts.get_fsky(ilc_dic)

Nothing need be written to disk: ``get_ilc_residuals.run_ilc`` returns the same dictionary that :func:`load_ilc_product` reads back, so a script can run component separation and forecasting in one process::

    import get_ilc_residuals

    ilc_dic, _ = get_ilc_residuals.run_ilc('s4wide_202310xx_pbdr_config', include_gal=1, which_gal_mask=2, save=False)

Notes
-----
Effective noise assembly and the CLASS_delens calls live in :mod:`delensing`, and the Fisher algebra in :mod:`fisher`.
This module holds only the input and output handling, the parameter file and the command line, so that neither of those modules depends on a file format.
"""

import os
import sys
import warnings

import numpy as np

sys.path.append( os.path.join( os.path.dirname(os.path.abspath(__file__)), 'modules' ) )
import delensing
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

PRIOR_PREFIX = 'prior_'
"""
Prefix marking a Gaussian prior width in the parameter file.
"""


# Parameter file

def load_fisher_params(paramfile=PARAMFILE):
    r"""
    Read and validate the forecasting parameter file.

    The file uses the same flat ``key = value`` dialect as ``params.ini`` and is read with :func:`misc.get_param_dict`.
    Three groups of entries are collected by prefix rather than by section, since that dialect has none: ``fid_<name>`` gives a fiducial cosmological parameter, ``step_<name>`` its derivative step size and ``prior_<name>`` the width of a Gaussian prior on it.
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

    Warns
    -----
    UserWarning
        Propagated from :func:`delensing.validate_cosmology` and if a prior names a parameter that is not being varied, since such a prior has no effect.

    Notes
    -----
    Cosmological values and step sizes are coerced to ``float`` by :func:`delensing.validate_cosmology`.
    That matters because :func:`misc.get_param_dict` returns an ``int`` for any value written without a fractional part, so ``fid_N_eff = 3.0`` would otherwise arrive as an integer.
    Coercing here rather than asking for a particular way of writing the file keeps that detail out of the user's way.

    Multipole ranges, the lensing reconstruction cuts, the satellite block and the option flags are returned unchanged, since they are consumed by :mod:`fisher` and by :func:`delensing.run_class` rather than here.
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

    settings['cosmo_fid'] = cosmo_fid
    settings['step_sizes'] = step_sizes
    settings['priors'] = priors
    settings['vary_params'] = vary_params

    return settings


# Input

def load_ilc_product(ilc_fname):
    r"""
    Read an ILC product file.

    Parameters
    ----------
    ilc_fname : str
        Path to a ``.npy`` file written by ``get_ilc_residuals.run_ilc``.

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
        If the file can only be read with ``encoding='latin1'``, which means it was written under Python 2.

    Notes
    -----
    The products are written with :func:`numpy.save` applied to a dictionary, so they are read back with ``allow_pickle`` and unwrapped from the zero-dimensional object array that :func:`numpy.load` returns.

    The default ASCII decoding is tried first and ``encoding='latin1'`` only as a fallback, rather than passing ``latin1`` unconditionally.

    This loader is deliberately the only place that knows the product layout; :func:`delensing.build_effective_noise` and the routines in :mod:`fisher` take the dictionary itself, so that they can equally be handed the return value of ``get_ilc_residuals.run_ilc`` without anything having been written to disk.
    """

    if not os.path.exists(ilc_fname):
        raise FileNotFoundError('No ILC product at %s' % (ilc_fname))

    try:
        contents = np.load(ilc_fname, allow_pickle=True)
    except UnicodeError:
        warnings.warn('%s could not be read with the default ASCII decoding, which means it was '
                      "written under Python 2; re-reading it with encoding='latin1'." % (ilc_fname),
                      stacklevel=2)
        contents = np.load(ilc_fname, allow_pickle=True, encoding='latin1')

    if not hasattr(contents, 'item') or contents.shape != ():
        raise ValueError('%s does not contain a single dictionary; got %s'
                         % (ilc_fname, type(contents).__name__ if not hasattr(contents, 'shape')
                            else 'an array of shape %s' % (contents.shape,)))
    ilc_dic = contents.item()
    if not isinstance(ilc_dic, dict):
        raise ValueError('%s contains a %s rather than a dictionary'
                         % (ilc_fname, type(ilc_dic).__name__))

    return ilc_dic


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
        If a supplied value overrides a stored one, since the stored value is the sky fraction the foreground spectra were computed on.

    Notes
    -----
    ``get_ilc_residuals.build_output_dic`` records ``'fsky_val'`` only when ``include_gal`` is set, so a product built without galactic foregrounds carries no sky fraction and one must be supplied.
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

    Notes
    -----
    This will not read the older products under ``products/``, whose layout is described in :func:`write_lensing_noise_curves`.

    TODO: This function will be extended to also read the older products.
    """

    return load_ilc_product(fname)


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

    The noise is stored in the canonical form rather than either of the re-indexed views, since :func:`delensing.noise_for_fisher` and :func:`delensing.noise_for_class` derive those cheaply and storing one of them would invite using it in the wrong place.

    Written and read with :func:`write_forecast_cache` and :func:`read_forecast_cache`. These files are large, since they hold every derivative at every multipole, so they belong under ``results/`` rather than in the repository.
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

    return load_ilc_product(fname)


# Driver

#...

#def main(argv=None):
#    # ...
#
#
#if __name__ == '__main__':
#    main()
