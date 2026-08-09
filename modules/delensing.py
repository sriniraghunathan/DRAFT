r"""
Iterative delensing of CMB spectra with CLASS_delens, driven through FisherLens.

This module supplies the delensing stage of DRAFT.
It takes the residual noise and foreground power spectra produced by the internal linear combination in :mod:`get_ilc_residuals`, combines them with the effective Planck instrumental noise, and hands the result to CLASS_delens, which iteratively reconstructs the lensing potential and returns the delensed CMB spectra together with the lensing-reconstruction noise :math:`N_\ell^{dd}`.
Those outputs are consumed by :mod:`fisher`, which performs the Fisher-matrix forecast.

The module is grouped into five sections:

* Setup: :func:`nearest_existing_dir`, :func:`resolve_paths`, :func:`import_fisherlens`, :func:`check_setup`
* Cosmology: :func:`validate_cosmology`
* Noise: :func:`planck_noise`, :func:`check_ilc_product`, :func:`build_effective_noise`, :func:`noise_for_fisher`, :func:`noise_for_class`
* CLASS diagnostics: :func:`missing_class_outputs`, :func:`missing_class_derivatives`, :func:`unused_class_parameters`, :func:`ng_derivative_bytes`, :func:`class_run_parameters`, :func:`derivative_flags_consumed`
* CLASS runs: :func:`check_noise_for_class`, :func:`deflection_noise_for_class`, :func:`reconstruction_mask`, :func:`run_class`, :func:`lensing_derivatives`, :func:`lensing_noise_curves`, :func:`parameter_derivatives`

The CLASS runs section has three entry points, one for each kind of CLASS_delens run.
:func:`run_class` returns the spectra and the reconstruction noise of a single cosmology, :func:`lensing_derivatives` the derivative matrices of the non-Gaussian covariance, and :func:`parameter_derivatives` the derivatives with respect to each varied parameter.
The rest of that section prepares their inputs and assembles their outputs, while the CLASS diagnostics section establishes what CLASS actually read and produced.

:func:`validate_cosmology`, :func:`resolve_paths`, :func:`build_effective_noise`, :func:`reconstruction_mask`, :func:`deflection_noise_for_class`, :func:`run_class`, :func:`ng_derivative_bytes`, :func:`class_run_parameters`, :func:`lensing_derivatives`, :func:`lensing_noise_curves` and :func:`parameter_derivatives` are called by :mod:`get_fisher_forecasts`, and :func:`import_fisherlens` and :func:`noise_for_fisher` by :mod:`fisher`.
:func:`check_setup` confirms that FisherLens and a compiled CLASS_delens are usable, and is called by :mod:`get_fisher_forecasts` before the first CLASS run of a forecast so that a fresh checkout reports every missing piece at once.
The remaining routines are internal helpers.

CLASS_delens and FisherLens are not distributed with DRAFT.
FisherLens is a git submodule at the repository root and brings CLASS_delens as a submodule of its own::

    git submodule update --init --recursive
    cd FisherLens/CLASS_delens
    make class

See :doc:`installation` for the full installation instructions.
Power spectra are in units of μK² throughout, matching the convention of :mod:`ilc` and the ``cl_residual`` entries of the ILC product files.
Multipoles are held under ``el``, as elsewhere in DRAFT, in every dictionary this module builds.
FisherLens uses ``l`` instead, so the dictionaries it builds and consumes keep that name and the conversion happens in one place, :func:`noise_for_class`.
"""

# Note:
# Several routines in this module work around behavior of CLASS_delens and FisherLens that was established empirically rather than from their documentation.
# Each such workaround carries a comment naming the pinned commit it was verified against since the reasoning is not otherwise recoverable from the code.
# The pinned commits are those recorded in the respective ``.gitmodules``; if either submodule pointer is advanced, those comments are the places to re-check.

import difflib
import os
import shutil
import sys
import tempfile
import uuid
import warnings

import numpy as np

import exp_specs
import misc


# Constants

DEFAULT_FISHERLENS_DIR = os.path.join(misc.REPO_ROOT, 'FisherLens')
"""
Default location of the FisherLens submodule.
"""

DEFAULT_CLASS_SUBDIR = 'CLASS_delens'
"""
Name of the CLASS_delens submodule inside FisherLens.
"""

DEFAULT_CLASS_DATA_DIR = os.path.join(misc.REPO_ROOT, 'results', 'class_data')
"""
Default scratch location for the CLASS input and output files.

These are small for a Gaussian forecast since every parameter step overwrites the previous one, but reach tens of gigabytes once the non-Gaussian lensing derivatives are requested.
Point this at a scratch filesystem in that case.
"""

CLASS_BINARY = 'class'
"""
Name of the compiled CLASS executable, produced by ``make class``.
"""

DELENSING_MODES = ('iterative', 'yes', 'no')
"""
Values CLASS_delens accepts for its ``delensing`` input.

``'iterative'`` re-reconstructs the lensing potential from iteratively delensing until the spectra converge, ``'yes'`` delenses once with the reconstruction obtained from the lensed spectra and ``'no'`` does not delens at all.
Only the first is matched by name in ``input.c`` of CLASS_delens, which sets ``has_itr_delensing``, and any other accepted value simply leaves iteration off.
"""

#List is extracted from the source of the pinned FisherLens commit (652eaec), because that routine reads its cosmology dictionary through a chain of membership tests and silently ignores anything it does not recognize, so a misspelled name would otherwise pass unnoticed and CLASS would quietly use its own default.
#If the FisherLens submodule pointer is advanced, this list is one of the places to re-check.
CLASS_COSMO_KEYS = (
    'A_s', 'Gamma_0_nadm', 'H0', 'N_eff', 'N_idr', 'Yhe',
    'alpha_ad_bi', 'alpha_ad_cdi', 'alpha_ad_nid', 'alpha_ad_niv', 'alpha_bi', 'alpha_bi_cdi',
    'alpha_bi_nid', 'alpha_bi_niv', 'alpha_cdi', 'alpha_cdi_nid', 'alpha_cdi_niv', 'alpha_nid',
    'alpha_nid_niv', 'alpha_niv', 'alpha_s', 'bbn_alpha_sensitivity', 'c_ad_bi', 'c_ad_cdi',
    'c_ad_nid', 'c_ad_niv', 'c_bi_cdi', 'c_bi_nid', 'c_bi_niv', 'c_cdi_nid',
    'c_cdi_niv', 'c_min', 'c_nid_niv', 'eta_0', 'fEDE', 'f_bi',
    'f_cdi', 'f_idm', 'f_nid', 'f_niv', 'f_scf', 'h',
    'log10f_scf', 'log10m_scf', 'log10z_c', 'logA', 'm_scf', 'mnu',
    'n_ad_bi', 'n_ad_cdi', 'n_ad_nid', 'n_ad_niv', 'n_bi', 'n_bi_cdi',
    'n_bi_nid', 'n_bi_niv', 'n_cdi', 'n_cdi_nid', 'n_cdi_niv', 'n_nid',
    'n_nid_niv', 'n_niv', 'n_s', 'n_t', 'omega_b_h2', 'omega_c_h2',
    'omk', 'r', 'tau', 'theta_s', 'thetai_scf', 'varying_alpha',
    'varying_me', 'varying_transition_redshift',
    )
"""
Cosmological parameter names that ``classWrapTools.class_generate_data`` recognizes.
"""

GUARD_DELTA_L_MAX = 2000
"""
Multipole margin by which ``classWrapTools`` requires an externally supplied noise spectrum to exceed the requested ``lmax``.

Its length check is applied before its ``extraParams`` override takes effect, so this value is fixed however ``delta_l_max`` is set, and it happens to coincide with the ``delta_l_max`` that ``classWrapTools`` itself installs.
CLASS reads the noise only up to its ``l_lensed_max``, which equals the requested ``lmax``, so the margin is never read and exists only to satisfy that check.
"""

CLASS_COSMO_KEYS_REQUIRED = ('omega_b_h2', 'omega_c_h2', 'n_s', 'tau')
"""
Cosmological parameters that ``classWrapTools.class_generate_data`` indexes unconditionally so that omitting one raises ``KeyError`` from inside FisherLens.
"""

LENSING_DERIVATIVE_FILES = {'cl_TT': 'dClTTdCldd', 'cl_EE': 'dClEEdCldd',
                            'cl_TE': 'dClTEdCldd', 'cl_BB': 'dClBBdCldd'}
"""
File stem CLASS_delens writes each lensing-deflection derivative to, keyed as ``classWrapTools.loadLensingDerivatives`` of FisherLens keys the dictionary it returns.
The full name is ``<root_name>_<stem>_<derivative type>.dat`` under ``output/``.
"""

UNLENSED_DERIVATIVE_FILES = {'cl_TT_cl_TT': 'dClTTdClTT', 'cl_TE_cl_TE': 'dClTEdClTE',
                             'cl_EE_cl_EE': 'dClEEdClEE', 'cl_EE_cl_BB': 'dClEEdClBB',
                             'cl_BB_cl_EE': 'dClBBdClEE', 'cl_BB_cl_BB': 'dClBBdClBB'}
"""
Same as :data:`LENSING_DERIVATIVE_FILES` for the derivatives with respect to the unlensed spectra, keyed as ``classWrapTools.loadUnlensedSpectraDerivatives`` keys them.
"""

DERIVATIVE_TYPES = ['lensed', 'delensed']
"""
Derivative types CLASS_delens will produce.
There is no unlensed entry: the unlensed spectra do not depend on the lensing potential, so those derivatives vanish and the unlensed non-Gaussian covariance is the Gaussian one.
"""

NG_DERIVATIVE_MODES = ['recompute', 'save', 'load']
"""
Ways of obtaining the non-Gaussian derivative matrices.
They are recomputed and discarded, recomputed and written to a DRAFT product, or read back from an earlier run without a CLASS run at all.
No default is offered since the three trade a large file against a long calculation.
"""

DERIVATIVE_FLAGS = ['delensing derivatives', 'output_derivatives', 'derv_binedges',
                    'derivative type', 'calculate_derviaties_wrt_unlensed',
                    'unlensed derivative type']
"""
CLASS parameters that ``classWrapTools`` sets to request the derivative matrices, spelled as it spells them.
The transposition in ``calculate_derviaties_wrt_unlensed`` is reproduced deliberately: if CLASS spells that name correctly then the flag is never consumed and the unlensed derivatives are silently not requested, which is what :func:`unused_class_parameters` is for.
"""


# Setup

def nearest_existing_dir(path):
    r"""
    Walk up a path until an existing directory is found.

    This identifies where a directory tree would actually start being created since :func:`os.makedirs` creates every missing level rather than only the last one.

    Parameters
    ----------
    path : str
        Path, which need not exist.

    Returns
    -------
    str
        Closest ancestor of ``path`` that exists, which in the worst case is the filesystem root.
    """

    path = os.path.abspath(path)
    while not os.path.exists(path):
        parent = os.path.dirname(path)
        if parent == path:  #reached the root without finding anything
            break
        path = parent

    return path


def resolve_paths(fisherlens_dir=None, class_exec_dir=None, class_data_dir=None):
    r"""
    Resolve the FisherLens, CLASS_delens and scratch locations to absolute paths.

    All three are returned as absolute paths ending in a separator.
    Both properties are required rather than cosmetic.
    CLASS is launched by FisherLens as ``cd <class_exec_dir> ; ./class <class_data_dir>/input/<root>.ini``, so a relative ``class_data_dir`` is resolved once against the working directory when FisherLens writes the file and again against ``class_exec_dir`` when CLASS reads it, and the two do not coincide.
    The trailing separator is required because FisherLens forms its paths by string concatenation, as in ``class_data_dir + 'input/'``.

    Parameters
    ----------
    fisherlens_dir : str, optional
        Location of the FisherLens checkout, which must contain ``classWrapTools.py``.
        Default is ``None``, which uses :data:`DEFAULT_FISHERLENS_DIR`.
    class_exec_dir : str, optional
        Directory holding the compiled CLASS executable.
        Default is ``None``, which uses the :data:`DEFAULT_CLASS_SUBDIR` folder inside ``fisherlens_dir``.
    class_data_dir : str, optional
        Scratch directory for the CLASS input and output files.
        Default is ``None``, which uses :data:`DEFAULT_CLASS_DATA_DIR`.

    Returns
    -------
    paths : dict
        Absolute paths ending in a separator, keyed ``'fisherlens_dir'``, ``'class_exec_dir'`` and ``'class_data_dir'``.
        Nothing is created and nothing is checked for existence here; see :func:`check_setup`.

    Notes
    -----
    A relative path is resolved against the current working directory, not against the repository root.
    """

    if fisherlens_dir is None:
        fisherlens_dir = DEFAULT_FISHERLENS_DIR
    fisherlens_dir = os.path.abspath( os.path.expanduser(fisherlens_dir) )

    if class_exec_dir is None:
        class_exec_dir = os.path.join(fisherlens_dir, DEFAULT_CLASS_SUBDIR)
    class_exec_dir = os.path.abspath( os.path.expanduser(class_exec_dir) )

    if class_data_dir is None:
        class_data_dir = DEFAULT_CLASS_DATA_DIR
    class_data_dir = os.path.abspath( os.path.expanduser(class_data_dir) )

    paths = {}
    for key, path in [('fisherlens_dir', fisherlens_dir),
                      ('class_exec_dir', class_exec_dir),
                      ('class_data_dir', class_data_dir)]:
        #FisherLens concatenates, e.g. class_data_dir + 'input/', so a trailing separator is required.
        paths[key] = os.path.join(path, '')

    return paths


def import_fisherlens(fisherlens_dir=None):
    r"""
    Import the FisherLens modules from a checkout that is not on the import path.

    Parameters
    ----------
    fisherlens_dir : str, optional
        Location of the FisherLens checkout. Default is ``None``, which uses :data:`DEFAULT_FISHERLENS_DIR`.

    Returns
    -------
    class_wrap_tools : module
        The imported ``classWrapTools``.
    fisher_tools : module
        The imported ``fisherTools``.

    Raises
    ------
    ImportError
        If either module cannot be imported, with a message distinguishing an unpopulated submodule, a missing CAMB installation and any other cause.
        Also if the modules were already imported from a different checkout since Python caches modules by name and only one FisherLens can be used per process.

    Notes
    -----
    ``fisherTools`` imports ``cambWrapTools`` at module scope, which imports CAMB, so CAMB must be installed even though nothing on the CLASS code path calls it.
    On the CLASS path, ``cambWrapTools`` is reached only from ``getParamDerivsBAOandH0`` and ``camb_class_generate_data``, neither of which this module uses.
    """

    if fisherlens_dir is None:
        fisherlens_dir = DEFAULT_FISHERLENS_DIR
    fisherlens_dir = os.path.abspath( os.path.expanduser(fisherlens_dir) )

    if not os.path.exists( os.path.join(fisherlens_dir, 'classWrapTools.py') ):
        raise ImportError(
            'No classWrapTools.py under %s. If the FisherLens folder is empty the submodule '
            'has not been populated; run "git submodule update --init --recursive" from the '
            'repository root.' % (fisherlens_dir))

    if fisherlens_dir not in sys.path:
        sys.path.insert(0, fisherlens_dir)

    try:
        import classWrapTools
        import fisherTools
    except ImportError as err:
        if 'camb' in str(err):
            raise ImportError(
                'FisherLens could not be imported because CAMB is missing (%s). fisherTools '
                'imports cambWrapTools at module scope, which imports CAMB, so it is required '
                'even though the CLASS code path never calls it. Install it with '
                '"pip install camb".' % (err)) from err
        raise ImportError('FisherLens could not be imported from %s: %s' % (fisherlens_dir, err)) from err

    #Python caches modules by name, so a second call naming a different checkout would silently
    #return the first one. Fail loudly instead of forecasting against the wrong FisherLens.
    for module in (classWrapTools, fisherTools):
        loaded_from = os.path.dirname( os.path.abspath(module.__file__) )
        if loaded_from != fisherlens_dir:
            raise ImportError(
                '%s was already imported from %s, not from the requested %s. Python caches '
                'modules by name, so only one FisherLens checkout can be used per process.'
                % (module.__name__, loaded_from, fisherlens_dir))

    return classWrapTools, fisherTools


def check_setup(fisherlens_dir=None, class_exec_dir=None, class_data_dir=None, verbose=True):
    r"""
    Check that FisherLens and a compiled CLASS_delens are usable.

    Every problem found is collected before raising so that a fresh checkout reports all of its
    missing pieces at once rather than one per attempt.

    Parameters
    ----------
    fisherlens_dir : str, optional
        Location of the FisherLens checkout. Default is ``None``, which uses :data:`DEFAULT_FISHERLENS_DIR`.
    class_exec_dir : str, optional
        Directory holding the compiled CLASS executable. Default is ``None``, which uses the :data:`DEFAULT_CLASS_SUBDIR` folder inside ``fisherlens_dir``.
    class_data_dir : str, optional
        Scratch directory for the CLASS files. Default is ``None``, which uses :data:`DEFAULT_CLASS_DATA_DIR`.
        Neither it nor any of its parents need exist since FisherLens creates ``input/`` and ``output/`` beneath it with :func:`os.makedirs`, which creates every missing level. What must be writable is therefore the nearest existing ancestor, which is checked.
    verbose : bool, optional
        Print the resolved paths and the outcome of each check. Default is ``True``.

    Returns
    -------
    paths : dict
        The resolved paths, as returned by :func:`resolve_paths`.

    Raises
    ------
    RuntimeError
        If FisherLens cannot be imported, CLASS has not been compiled or the scratch directory cannot be created, listing every problem found.

    Notes
    -----
    The scratch-directory test uses :func:`os.access`, which consults the real rather than the effective user ID, so it is a cheap pre-flight check rather than a guarantee.
    """

    paths = resolve_paths(fisherlens_dir=fisherlens_dir, class_exec_dir=class_exec_dir, class_data_dir=class_data_dir)
    problems = []

    if verbose:
        print('FisherLens : %s' % (paths['fisherlens_dir']))
        print('CLASS exec : %s' % (paths['class_exec_dir']))
        print('CLASS data : %s' % (paths['class_data_dir']))

    try:
        import_fisherlens(paths['fisherlens_dir'])
        if verbose:
            print('  FisherLens imports: yes')
    except ImportError as err:
        problems.append( str(err) )
        if verbose:
            print('  FisherLens imports: NO')

    class_binary = os.path.join(paths['class_exec_dir'], CLASS_BINARY)
    if not os.path.exists(class_binary):
        problems.append(
            'No CLASS executable at %s. Compile it with "cd %s ; make class". Use the "class" '
            'target rather than "make" or "make all", which also build the "classy" Python '
            'wrapper and need Cython.' % (class_binary, paths['class_exec_dir']))
        if verbose:
            print('  CLASS compiled: NO')
    elif not os.access(class_binary, os.X_OK):
        problems.append('The CLASS executable %s is not executable; check its permissions.' % (class_binary))
        if verbose:
            print('  CLASS compiled: yes, but not executable')
    elif verbose:
        print('  CLASS compiled: yes')

    scratch_ancestor = nearest_existing_dir(paths['class_data_dir'])
    if not os.access(scratch_ancestor, os.W_OK):
        problems.append(
            'The CLASS scratch directory %s cannot be created: its nearest existing ancestor %s '
            'is not writable.' % (paths['class_data_dir'], scratch_ancestor))
        if verbose:
            print('  scratch writable: NO')
    elif verbose:
        print('  scratch writable: yes')

    if problems:
        raise RuntimeError('Setup is incomplete:\n  - %s' % ('\n  - '.join(problems)))

    if verbose:
        print('Setup looks complete.')

    return paths


# Cosmology

def validate_cosmology(cosmo_fid, step_sizes=None, vary_params=None):
    r"""
    Check a fiducial cosmology and its derivative step sizes, and normalize their types.

    Parameters
    ----------
    cosmo_fid : dict
        Fiducial parameter values, keyed by the names ``classWrapTools`` uses. Every parameter CLASS should know about belongs here, including those held fixed.
    step_sizes : dict, optional
        Step sizes for the numerical derivatives, keyed the same way. Required for every varied parameter and ignored for the rest. Default is ``None``, which is treated as an empty dictionary.
    vary_params : list of str, optional
        Parameters to differentiate. Default is ``None``, which varies every entry of ``cosmo_fid``, matching the behavior of ``fisherTools.getPowerDerivWithParams``.

    Returns
    -------
    cosmo_fid : dict
        The fiducial values, with each coerced to ``float``.
    step_sizes : dict
        The step sizes of the varied parameters only, each coerced to ``float``.
    vary_params : list of str
        The parameters to differentiate, as a list.

    Raises
    ------
    ValueError
        If ``cosmo_fid`` is empty, if it holds a name CLASS would not recognize, if one of :data:`CLASS_COSMO_KEYS_REQUIRED` is missing, if ``mnu`` is given without ``N_eff``, if a varied parameter is absent from ``cosmo_fid``, or if a varied parameter has no step size or a non-positive one.

    Warns
    -----
    UserWarning
        If neither ``A_s`` nor ``logA`` is given, or if none of ``H0``, ``h`` and ``theta_s`` is given. CLASS then silently uses its own default.
        Also if a step size would carry a positive parameter to zero or below since the two-sided derivative would then straddle an unphysical value.

    Notes
    -----
    The unrecognized-name check is important. ``classWrapTools.class_generate_data`` of FisherLens reads its cosmology through tests of the form ``if 'x' in cosmo``, so a misspelling is dropped without complaint and CLASS falls back to its own default: writing ``N_effective`` in place of ``N_eff`` yields a forecast at :math:`N_\mathrm{eff} = 3.046` rather than the requested value, with nothing said.

    ``N_eff`` is required as soon as ``mnu`` is present because ``classWrapTools.class_generate_data`` of FisherLens then computes ``N_ur`` by subtracting the contribution of one massive species from it and indexes ``N_eff`` unconditionally along that branch.

    Every supplied step size is returned, but only those of the varied parameters are checked. ``fisherTools.getPowerDerivWithParams`` indexes ``stepSizes[name]`` for each parameter it differentiates and additionally indexes ``stepSizes['mnu']`` whenever ``mnu`` appears in the cosmology to decide whether the neutrino mass needs a one-sided derivative. A step size for ``mnu`` is therefore required even when it is held fixed.

    Parameters held fixed are handled by leaving them out of ``vary_params`` while keeping them in ``cosmo_fid``. They then reach CLASS but are not differentiated.
    """

    if not isinstance(cosmo_fid, dict) or len(cosmo_fid) == 0:
        raise ValueError('cosmo_fid must be a non-empty dictionary of parameter values')
    if step_sizes is None:
        step_sizes = {}
    if vary_params is None:
        vary_params = list(cosmo_fid.keys())
    vary_params = list(vary_params)

    unrecognized = [name for name in cosmo_fid if name not in CLASS_COSMO_KEYS]
    if unrecognized:
        hints = []
        for name in sorted(unrecognized):
            close = difflib.get_close_matches(name, CLASS_COSMO_KEYS, n=2, cutoff=0.6)
            hints.append('%s%s' % (name, ' (did you mean %s?)' % (' or '.join(close)) if close else ''))
        raise ValueError('These names are not recognized by classWrapTools.class_generate_data of '
                         'FisherLens and would be silently ignored, leaving CLASS to use its own '
                         'defaults: %s' % (', '.join(hints)))

    missing = [name for name in CLASS_COSMO_KEYS_REQUIRED if name not in cosmo_fid]
    if missing:
        raise ValueError('cosmo_fid is missing %s, which classWrapTools.class_generate_data of '
                         'FisherLens indexes unconditionally' % (', '.join(missing)))
    if 'mnu' in cosmo_fid and 'N_eff' not in cosmo_fid:
        raise ValueError("cosmo_fid gives mnu but not N_eff. classWrapTools.class_generate_data of "
                         'FisherLens derives N_ur from N_eff whenever mnu is present, so omitting it'
                         'raises KeyError inside FisherLens.')

    if not any(name in cosmo_fid for name in ('A_s', 'logA')):
        warnings.warn('cosmo_fid gives neither A_s nor logA, so CLASS will use its own primordial '
                      'amplitude rather than a requested one.', stacklevel=2)
    if not any(name in cosmo_fid for name in ('H0', 'h', 'theta_s')):
        warnings.warn('cosmo_fid gives none of H0, h and theta_s, so CLASS will use its own expansion '
                      'rate rather than a requested one.', stacklevel=2)

    absent = [name for name in vary_params if name not in cosmo_fid]
    if absent:
        raise ValueError('These parameters are to be varied but have no fiducial value in cosmo_fid: '
                         '%s' % (', '.join(sorted(absent))))

    cosmo_fid = {name: float(value) for name, value in cosmo_fid.items()}
    #every supplied step size is kept, not only those of the varied parameters: FisherLens indexes
    #stepSizes['mnu'] unconditionally whenever mnu is in the cosmology, to decide whether to take a
    #one-sided derivative, so dropping it raises KeyError even when mnu is held fixed
    checked_steps = {name: float(value) for name, value in step_sizes.items()}
    for name in vary_params:
        if name not in step_sizes:
            raise ValueError('No step size for %s, which is to be varied. A two-sided derivative '
                             'needs one.' % (name))
        step = float(step_sizes[name])
        if not step > 0.:
            raise ValueError('The step size for %s must be positive, got %s. getPowerDerivWithParams '
                             'from FisherLens divides by it.' % (name, step))
        if cosmo_fid[name] > 0. and cosmo_fid[name] - step <= 0.:
            warnings.warn('The step size for %s is %g against a fiducial value of %g, so the '
                          'two-sided derivative would evaluate it at %g. If it cannot be negative, '
                          'reduce the step; getPowerDerivWithParams from the pinned FisherLens takes '
                          'a one-side derivative only for mnu and a short hardcoded list.'
                          % (name, step, cosmo_fid[name], cosmo_fid[name]-step), stacklevel=2)
        checked_steps[name] = step

    return cosmo_fid, checked_steps, vary_params


# Noise

def planck_noise(el, expname='planck', include_pol=True):
    r"""
    Effective Planck instrumental noise, inverse-variance combined over frequency bands.

    The per-band beams and white-noise levels are taken from :func:`exp_specs.get_exp_specs` and the per-band spectra are formed with :func:`misc.get_nl`.
    Low-frequency (1/f) noise is neglected for a satellite, which the stored specifications express by setting :math:`\ell_\mathrm{knee}` to ``-1`` in every band so that :math:`N_\ell = \Delta_X^2 B_\ell^{-2}` with beam function :math:`B_\ell`.
    The bands are then combined as

    .. math::

        \frac{1}{N_\ell} = \sum_i \frac{1}{N_{\ell, i}}\, .

    Parameters
    ----------
    el : array_like
        Multipole moments :math:`\ell` at which to evaluate the noise. May start at zero.
    expname : str, optional
        Configuration name passed to :func:`exp_specs.get_exp_specs`. Default is ``'planck'``.
    include_pol : bool, optional
        Include the polarization sensitivity. When ``False`` the polarization noise is infinite, i.e. the experiment contributes to temperature only. Default is ``True``.

    Returns
    -------
    dict
        Multipoles under ``'el'``, and noise power spectra in μK² under ``'cl_TT'``, ``'cl_EE'``, ``'cl_BB'`` and ``'cl_TE'``.

    Notes
    -----
    ``fisherTools.getPlanckInvVarNoise`` computes the same quantity from the same table.
    It is not used here for three reasons: it would put the band specifications outside :mod:`exp_specs`, it would make this function depend on FisherLens for a purely numerical step and it returns spectra that are not everywhere finite.
    The two routes were compared directly and effectively agree.
    """

    specs_dic = exp_specs.get_exp_specs(expname)[0]

    el = np.asarray(el)
    inv_nl = {'T': np.zeros(len(el)), 'P': np.zeros(len(el))}
    #the largest exponent np.exp can represent; beyond it a band contributes nothing
    max_exponent = np.log( np.finfo(float).max )

    for freq in sorted(specs_dic):
        beam_arcmins, delta_T, elknee_T, alphaknee_T, delta_P, elknee_P, alphaknee_P = specs_dic[freq]
        sigma = np.radians(beam_arcmins/60.) / np.sqrt(8.*np.log(2.))
        for TP, noiseval, elknee, alphaknee in [('T', delta_T, elknee_T, alphaknee_T),
                                                ('P', delta_P, elknee_P, alphaknee_P)]:
            if TP == 'P' and not include_pol:
                continue
            noiseval_radians = noiseval * np.radians(1./60.)
            #misc.get_nl forms exp( l(l+1) sigma**2 ) and only then multiplies by noiseval_radians**2.
            # Evaluate the band only where both are representable, which keeps an overflow from arising.
            headroom = max_exponent - max(0., 2.*np.log(noiseval_radians))
            safe = el * (el+1.) * sigma**2. < headroom
            if not safe.any():
                continue
            nl = misc.get_nl(noiseval, el[safe], beam_arcmins, elknee=elknee, alphaknee=alphaknee)
            contribution = np.zeros(len(el))
            contribution[safe] = np.divide(1., nl, out=np.zeros_like(nl), where=nl > 0.)
            inv_nl[TP] = inv_nl[TP] + contribution

    nl_dic = {}
    for TP in ('T', 'P'):
        nl_dic[TP] = np.divide(1., inv_nl[TP], out=np.full(len(el), np.inf), where=inv_nl[TP] > 0.)

    return {'el': el,
            'cl_TT': nl_dic['T'],
            'cl_EE': nl_dic['P'],
            'cl_BB': np.copy(nl_dic['P']),  #polarization noise taken to be mode independent
            'cl_TE': np.zeros(len(el))}


def check_ilc_product(ilc_dic, lmax_calc, lmin_ilc):
    r"""
    Validate an ILC product against the requested multipole range.

    Parameters
    ----------
    ilc_dic : dict
        ILC product, as written by :func:`get_ilc_residuals.build_output_dic` or returned by :func:`get_ilc_residuals.run_ilc`.
        Must hold ``'el'`` and ``'cl_residual'``, the latter with ``'TT'`` and ``'EE'`` entries.
    lmax_calc : int
        Highest multipole at which the effective noise is needed, which is the ``lmax`` handed to CLASS.
    lmin_ilc : int
        Lowest multipole at which the ILC residual is to be used.

    Returns
    -------
    first_nonzero : int
        Lowest multipole at which both the ``TT`` and ``EE`` residuals are positive.

    Raises
    ------
    ValueError
        If required entries are missing, if the multipoles are not consecutive integers starting at zero, if the ILC does not reach ``lmax_calc``, if ``lmin_ilc`` reaches into the region where the residual vanishes or if the residual is not positive everywhere in the range that will be used.

    Warns
    -----
    UserWarning
        If the multipole below which the residual vanishes disagrees with ``param_dict['lmin']`` recorded in the product.

    Notes
    -----
    The residual is zero for :math:`\ell \le` ``param_dict['lmin']`` because :func:`get_ilc_residuals.get_noise_spectra` zeroes the noise there, so an inverse-variance combination that reaches those multipoles would divide by zero.
    Rather than trusting the recorded ``lmin``, the vanishing region is measured from the residual itself and the two are compared.
    """

    for key in ('el', 'cl_residual'):
        if key not in ilc_dic:
            raise ValueError("The ILC product has no '%s' entry; keys are %s" % (key, sorted(ilc_dic.keys())))
    for spec in ('TT', 'EE'):
        if spec not in ilc_dic['cl_residual']:
            raise ValueError("The ILC product has no cl_residual['%s']; keys are %s"
                             % (spec, sorted(ilc_dic['cl_residual'].keys())))

    el = np.asarray(ilc_dic['el'])
    if el.ndim != 1 or len(el) == 0:
        raise ValueError('The ILC multipoles must be a non-empty one-dimensional array, got shape %s' % (el.shape,))
    if el[0] != 0 or not np.array_equal( el, np.arange(len(el)) ):
        raise ValueError('The ILC multipoles must be consecutive integers starting at zero, got %s ... %s'
                         % (el[:4], el[-2:]))

    cl_tt = np.asarray(ilc_dic['cl_residual']['TT'], dtype=float)
    cl_ee = np.asarray(ilc_dic['cl_residual']['EE'], dtype=float)
    for spec, cl in [('TT', cl_tt), ('EE', cl_ee)]:
        if len(cl) != len(el):
            raise ValueError("cl_residual['%s'] has %d entries but there are %d multipoles"
                             % (spec, len(cl), len(el)))

    if el[-1] < lmax_calc:
        raise ValueError('The ILC product reaches only ell = %d, but effective noise is needed to '
                         'ell = %d. Rerun the ILC with a higher lmax, or lower lmax and lbuffer.'
                         % (el[-1], lmax_calc))

    positive = (cl_tt > 0.) & (cl_ee > 0.)
    if not positive.any():
        raise ValueError('The ILC residual is nowhere positive in both TT and EE.')
    first_nonzero = int( el[positive][0] )

    recorded_lmin = ilc_dic.get('param_dict', {}).get('lmin')
    if recorded_lmin is not None and first_nonzero != recorded_lmin + 1:
        warnings.warn('The ILC residual first becomes positive at ell = %d, but the product records '
                      "param_dict['lmin'] = %s, which implies ell = %d. Using the measured value."
                      % (first_nonzero, recorded_lmin, recorded_lmin + 1), stacklevel=2)

    if lmin_ilc < first_nonzero:
        raise ValueError('lmin_ilc = %d reaches into the range where the ILC residual vanishes, which '
                         'would divide by zero in the inverse-variance combination. It must be at '
                         'least %d for this product.' % (lmin_ilc, first_nonzero))
    if lmin_ilc > lmax_calc:
        raise ValueError('lmin_ilc = %d exceeds lmax_calc = %d, so the ILC would never be used.'
                         % (lmin_ilc, lmax_calc))

    used = np.arange(lmin_ilc, lmax_calc + 1)
    if not positive[used].all():
        bad = used[~positive[used]]
        raise ValueError('The ILC residual is not positive at %d multipoles in the range that would be '
                         'used, starting at ell = %d.' % (len(bad), bad[0]))

    return first_nonzero


def build_effective_noise(ilc_dic, lmax_calc, lmin_ilc=30, guard_pad=GUARD_DELTA_L_MAX, include_planck=True, planck_expname='planck', verbose=False):
    r"""
    Effective post-component-separation noise to hand to CLASS_delens.

    The ILC residual supplies the effective noise above ``lmin_ilc`` and is combined with the effective Planck noise by inverse-variance weighting,

    .. math::

        N_\ell = \left( \frac{1}{C_\ell^\mathrm{res}} + \frac{1}{N_\ell^\mathrm{Planck}} \right)^{-1} ,

    while below ``lmin_ilc`` the Planck noise is used on its own.
    That cut serves two purposes: it excludes the largest angular scales, which are most affected by systematic effects in a ground-based survey, and it keeps the combination away from the low multipoles where the ILC residual vanishes identically.

    Parameters
    ----------
    ilc_dic : dict
        ILC product, validated by :func:`check_ilc_product`.
    lmax_calc : int
        Highest multipole at which real noise is required, which is the ``lmax`` handed to CLASS.
        The caller is responsible for choosing it above the multipole range used in the Fisher sum; see the Notes.
    lmin_ilc : int, optional
        Lowest multipole at which the ILC residual is used. Default is 30, following the accompanying paper.
    guard_pad : int, optional
        Number of dead multipoles appended above ``lmax_calc``. Default is :data:`GUARD_DELTA_L_MAX`; see the Notes.
    include_planck : bool, optional
        Combine in the Planck noise. When ``False`` the ILC residual is used alone above ``lmin_ilc`` and the noise below it is infinite, which excludes those multipoles rather than filling them with Planck. Default is ``True``.
    planck_expname : str, optional
        Configuration name for the satellite whose noise is combined in, passed to :func:`planck_noise`. Default is ``'planck'``.
    verbose : bool, optional
        Print a summary of the assembled noise. Default is ``False``.

    Returns
    -------
    dict
        Multipoles under ``'el'``, running from 0 to ``lmax_calc + guard_pad`` inclusive, and noise power spectra in μK² under ``'cl_TT'``, ``'cl_EE'``, ``'cl_BB'`` and ``'cl_TE'``.
        See :func:`noise_for_class` and :func:`noise_for_fisher` for the two forms its consumers need.

    Raises
    ------
    ValueError
        Propagated from :func:`check_ilc_product`, or if ``lmax_calc`` or ``guard_pad`` is negative.

    Notes
    -----
    The multipoles start at zero.
    CLASS_delens stores an externally supplied noise spectrum indexed by the *row* of the file it is read from, in ``external_temperature_noise_spectrum_init``, while the equivalent internal routine indexes by multipole and every consumer reads it by multipole.
    The array holding the multipole column of the file is written but never read, so there is no resampling step.
    Starting the multipoles at zero therefore makes the row index coincide with the multipole and the two paths agree.
    This was verified against CLASS_delens ``f27ff1b`` by injecting a step into the supplied noise and reading it back from ``_spectra_noise.dat``.
    If that submodule pin is advanced, this is one of the places to re-check.

    The padding above ``lmax_calc`` is never read.
    ``classWrapTools`` refuses noise that does not extend to ``lmax + 2000``, testing against a hardcoded value before its ``extraParams`` override is applied, so the array must be that long.
    CLASS itself reads the noise only up to ``l_lensed_max``, which equals the ``lmax`` it was given, so everything above that is dead weight and its values are immaterial.
    The last real value is repeated rather than extrapolated, to make that explicit.

    The :math:`TE` noise is zero.
    The ILC does compute a :math:`TE` residual when it solves for temperature and polarization jointly, but neither consumer can accept one: ``classWrapTools`` writes only a temperature and a polarization noise file, and ``fisherTools.getGaussianCov`` builds the :math:`TE` covariance with no :math:`N_\ell^{TE}` term.
    Zero is therefore what both codes assume.

    The :math:`BB` noise is set equal to the :math:`EE` noise, as is standard for a survey whose polarization noise is not mode dependent.

    ``lmax_calc`` is not compared against the multipole range of the Fisher sum, because that range is not known here.
    The caller should choose ``lmax_calc`` above it by a margin of at least a few hundred multipoles so that the lensing convolution is accurate at the highest multipole used.
    """

    if lmax_calc < 0:
        raise ValueError('lmax_calc must be non-negative, got %s' % (lmax_calc))
    if guard_pad < 0:
        raise ValueError('guard_pad must be non-negative, got %s' % (guard_pad))

    lmax_calc = int(lmax_calc)
    lmin_ilc = int(lmin_ilc)
    first_nonzero = check_ilc_product(ilc_dic, lmax_calc, lmin_ilc)

    #row index must equal multipole; see the Notes and CLASS_delens f27ff1b
    el = np.arange(0, lmax_calc + int(guard_pad) + 1)

    if include_planck:
        planck = planck_noise(el, expname=planck_expname)
        nl_tt = np.copy(planck['cl_TT'])
        nl_ee = np.copy(planck['cl_EE'])
    else:
        #no Planck means no information below lmin_ilc, which infinite noise expresses
        nl_tt = np.full(len(el), np.inf)
        nl_ee = np.full(len(el), np.inf)

    ilc_tt = np.asarray(ilc_dic['cl_residual']['TT'], dtype=float)
    ilc_ee = np.asarray(ilc_dic['cl_residual']['EE'], dtype=float)
    combine = (el >= lmin_ilc) & (el <= lmax_calc)
    for nl, ilc in [(nl_tt, ilc_tt), (nl_ee, ilc_ee)]:
        if include_planck:
            nl[combine] = 1. / ( 1./ilc[el[combine]] + 1./nl[combine] )
        else:
            nl[combine] = ilc[el[combine]]

    #dead padding: the values above lmax_calc are never read, so repeat rather than extrapolate
    pad = el > lmax_calc
    if pad.any():
        nl_tt[pad] = nl_tt[el == lmax_calc][0]
        nl_ee[pad] = nl_ee[el == lmax_calc][0]

    noise_dic = {'el': el,
                 'cl_TT': nl_tt,
                 'cl_EE': nl_ee,
                 'cl_BB': np.copy(nl_ee),  #polarization noise taken to be mode independent
                 'cl_TE': np.zeros(len(el))}  #zero is what both consumers assume; see the Notes

    if verbose:
        print('effective noise: ell = %d to %d, real up to %d, %d dead multipoles appended'
              % (el[0], el[-1], lmax_calc, int(pad.sum())))
        print('  ILC residual used from ell = %d (first positive at %d), Planck %s'
              % (lmin_ilc, first_nonzero, 'combined in' if include_planck else 'not included'))
        for spec, nl in [('TT', nl_tt), ('EE', nl_ee)]:
            finite = np.isfinite(nl[:lmax_calc+1])
            print('  %s: %d non-finite below lmax_calc, min %.4g, at ell=%d %.4g'
                  % (spec, int((~finite).sum()), nl[:lmax_calc+1][finite].min(),
                     lmax_calc, nl[lmax_calc]))

    return noise_dic


def noise_for_fisher(noise_dic):
    r"""
    Re-index an effective-noise dictionary for the Fisher assembly.

    Parameters
    ----------
    noise_dic : dict
        Effective noise from :func:`build_effective_noise`, whose entries begin at :math:`\ell = 0`.

    Returns
    -------
    dict
        The same entries with the first two multipoles dropped so that element zero is :math:`\ell = 2`.

    Raises
    ------
    ValueError
        If the supplied multipoles do not begin at zero since the offset would then be wrong.

    Notes
    -----
    ``fisherTools.getGaussianCov`` adds the fiducial spectra and the noise element by element, and the fiducial spectra come from the CLASS output files, which begin at :math:`\ell = 2`.
    The noise handed to the Fisher calculation must therefore begin at :math:`\ell = 2` as well, whereas the noise written out for CLASS must begin at :math:`\ell = 0` so that its row index coincides with the multipole.
    Both requirements are met by building one array from :math:`\ell = 0` and slicing it here, rather than assembling the noise twice.

    The returned entries are views on the input arrays, not copies. Nothing in FisherLens writes to a caller-supplied noise dictionary, so sharing the underlying data is safe.
    """

    el = np.asarray(noise_dic['el'])
    if len(el) < 3 or el[0] != 0:
        raise ValueError('The effective noise must start at ell = 0 and cover at least ell = 2, got '
                         '%d entries starting at %s' % (len(el), el[0] if len(el) else None))

    return {key: np.asarray(value)[2:] for key, value in noise_dic.items()}


def noise_for_class(noise_dic):
    r"""
    Re-key an effective-noise dictionary for the FisherLens CLASS wrapper.

    Parameters
    ----------
    noise_dic : dict
        Effective noise from :func:`build_effective_noise`, with multipoles under ``'el'``.

    Returns
    -------
    dict
        The same entries with the multipoles moved to ``'l'``, which is the name ``classWrapTools`` requires.

    Raises
    ------
    ValueError
        If the supplied dictionary has no ``'el'`` entry.

    Notes
    -----
    This is the single place where the DRAFT naming meets the FisherLens naming.
    ``classWrapTools.class_generate_data`` reads ``cmbNoise['l']`` in three places, to check the length of the array and to write the two noise files, and that is the only routine in FisherLens that reads the multipole entry of a noise dictionary at all; ``getGaussianCov`` takes the multipoles from the fiducial spectra instead.
    Converting here rather than naming the arrays ``'l'`` from the outset keeps every dictionary this module builds consistent with the rest of DRAFT, where multipoles are ``el``.

    The arrays are passed through by reference, not copied.
    """

    if 'el' not in noise_dic:
        raise ValueError("The effective noise has no 'el' entry; its keys are %s"
                         % (sorted(noise_dic.keys())))

    converted = {key: value for key, value in noise_dic.items() if key != 'el'}
    converted['l'] = noise_dic['el']

    return converted


# CLASS diagnostics

def missing_class_outputs(class_data_dir, root_name, expect_reconstruction):
    r"""
    List the CLASS output files that ``classWrapTools`` will read but which are absent.

    Parameters
    ----------
    class_data_dir : str
        Scratch directory handed to CLASS, ending in a separator.
    root_name : str
        Base name of the run.
    expect_reconstruction : bool
        Whether a lensing-reconstruction noise file is expected, which is the case when CLASS was asked to reconstruct iteratively.

    Returns
    -------
    list of str
        Names of the missing files, empty if all are present.
    """

    #CLASS appends an underscore to the root it is given (in ``input.c``), which is why the names carry one and why ``classWrapTools`` does not add it itself.
    suffixes = ['_cl.dat', '_cl_lensed.dat', '_cl_delensed.dat', '_spectra_noise.dat']
    if expect_reconstruction:
        suffixes.append('_lensing_noise_rcn.dat')

    return [root_name + suffix for suffix in suffixes
            if not os.path.exists( os.path.join(class_data_dir, 'output', root_name + suffix) )]


def missing_class_derivatives(class_data_dir, root_name, derivative_type, include_unlensed=True):
    r"""
    List the derivative files that ``classWrapTools`` will read but which are absent.

    Parameters
    ----------
    class_data_dir : str
        Scratch directory handed to CLASS, ending in a separator.
    root_name : str
        Base name of the run.
    derivative_type : str
        One of :data:`DERIVATIVE_TYPES`.
    include_unlensed : bool, optional
        Whether the derivatives with respect to the unlensed spectra were requested. Default is ``True``.

    Returns
    -------
    list of str
        Names of the missing files, empty if all are present.

    Notes
    -----
    Kept separate from :func:`missing_class_outputs` rather than folded into it because the two are checked at different points: a derivative run returns only the derivative matrices, so the spectrum files it also writes are never read back by it. It does write them since ``classWrapTools`` sets the same ``output`` on both kinds of run and the derivative branch does not alter it, so :func:`missing_class_outputs` applies to a derivative run unchanged.
    """

    #TODO: establish the spectrum files left behind

    stems = dict(LENSING_DERIVATIVE_FILES)
    if include_unlensed:
        stems.update(UNLENSED_DERIVATIVE_FILES)
    names = ['%s_%s_%s.dat' % (root_name, stem, derivative_type)
             for stem in sorted( stems.values() )]

    return [name for name in names
            if not os.path.exists( os.path.join(class_data_dir, 'output', name) )]


def unused_class_parameters(class_data_dir, root_name, among=None):
    r"""
    List the parameters CLASS read but did not consume.

    Parameters
    ----------
    class_data_dir : str
        Scratch directory handed to CLASS, ending in a separator.
    root_name : str
        Base name of the run.
    among : iterable of str, optional
        Restrict the result to these names. Default is ``None``, which returns all of them.

    Returns
    -------
    list of str
        Names CLASS ignored, sorted and without repeats.

    Raises
    ------
    FileNotFoundError
        If the file is absent, which means CLASS was not asked to record its parameters.

    Notes
    -----
    ``classWrapTools`` sets ``write parameters``, so CLASS records the deck it was given in ``<root_name>_parameters.ini`` and everything it did not recognize in ``<root_name>_unused_parameters``.
    That second file is the only way to tell a misspelled CLASS parameter from a working one: CLASS ignores what it does not recognize without complaint, so a flag whose name is wrong has no effect and no error.
    """

    path = os.path.join(class_data_dir, 'output', root_name + '_unused_parameters')
    if not os.path.exists(path):
        raise FileNotFoundError('No %s. CLASS records the parameters it ignored there when it is '
                                'asked to write its parameters, which classWrapTools of FisherLens '
                                'does, so its absence means the run did not reach that point.'
                                % (path))

    names = []
    with open(path) as handle:
        for line in handle:
            line = line.split('#')[0].strip()
            if '=' in line:
                names.append( line.split('=')[0].strip() )
    if among is not None:
        among = set(among)
        names = [name for name in names if name in among]

    return sorted( set(names) )


def ng_derivative_bytes(lmax_calc, include_unlensed=True, derivative_types=1):
    r"""
    Size in memory of a set of non-Gaussian derivative matrices.

    Parameters
    ----------
    lmax_calc : int
        Highest multipole the matrices cover.
    include_unlensed : bool, optional
        Whether the derivatives with respect to the unlensed spectra are counted. Default is ``True``.
    derivative_types : int, optional
        How many of :data:`DERIVATIVE_TYPES` are held at once. Default is 1.

    Returns
    -------
    int
        Size in bytes.

    Notes
    -----
    Each matrix is square in the multipole, of side ``lmax_calc + 1`` once the two rows and columns below the first multipole CLASS returns are padded on and is held as ``float64``.
    The size is therefore quadratic in ``lmax_calc``.
    """

    matrices = len(LENSING_DERIVATIVE_FILES)
    if include_unlensed:
        matrices += len(UNLENSED_DERIVATIVE_FILES)

    return (int(lmax_calc)+1)**2 * 8 * matrices * int(derivative_types)


def class_run_parameters(class_data_dir, root_name):
    r"""
    Read back the input deck CLASS recorded for a run.

    Parameters
    ----------
    class_data_dir : str
        Scratch directory handed to CLASS, ending in a separator.
    root_name : str
        Base name of the run.

    Returns
    -------
    dict
        Every parameter CLASS read, as written, values left as strings.

    Raises
    ------
    FileNotFoundError
        If the file is absent, which means the run did not reach the point of writing it.

    Notes
    -----
    The cosmological parameters appear under the names CLASS uses, not the ones DRAFT uses: ``omega_c_h2`` is recorded as ``omega_cdm``, ``tau`` as ``tau_reio``, ``theta_s`` as ``100*theta_s`` and ``N_eff`` as a combination of ``N_ur``, ``N_ncdm``, ``T_ncdm`` and ``m_ncdm``. Recovering a DRAFT cosmology from this file therefore means reproducing the translation inside ``classWrapTools``, which is why it is not used to check whether a stored product matches a requested one; a manifest written alongside the product is used for that instead. This file is read only for the parameters whose names DRAFT sets directly.
    """

    #Note: ``classWrapTools`` of FisherLens sets ``write parameters``, so CLASS records what it read in ``<root_name>_parameters.ini`` and what it did not in ``<root_name>_unused_parameters``. Between them the two files account for the whole deck that was sent.
    #``root`` is recorded as an absolute path, so it differs between machines and after any move of the working directory, and is never a meaningful thing to compare.

    path = os.path.join(class_data_dir, 'output', root_name + '_parameters.ini')
    if not os.path.exists(path):
        raise FileNotFoundError('No %s. CLASS writes the deck it read there, so its absence means '
                                'the run did not get that far.' % (path))

    recorded = {}
    with open(path) as handle:
        for line in handle:
            line = line.split('#')[0].strip()
            if '=' in line:
                key, value = line.split('=', 1)
                recorded[ key.strip() ] = value.strip()

    return recorded


def derivative_flags_consumed(class_data_dir, root_name, derivative_type, include_unlensed=True):
    r"""
    Check that CLASS both received and consumed the flags requesting the derivative matrices.

    Parameters
    ----------
    class_data_dir : str
        Scratch directory handed to CLASS, ending in a separator.
    root_name : str
        Base name of the run.
    derivative_type : str
        One of :data:`DERIVATIVE_TYPES`, checked against what CLASS recorded.
    include_unlensed : bool, optional
        Whether the derivatives with respect to the unlensed spectra were requested. Default is ``True``.

    Returns
    -------
    list of str
        Complaints, one per flag that was ignored, absent or recorded with an unexpected value. Empty when everything landed.

    Notes
    -----
    A CLASS parameter whose name it does not recognize is ignored without complaint, so a misspelled flag has no effect and raises nothing. ``calculate_derviaties_wrt_unlensed`` carries a transposition in ``classWrapTools``.
    """

    complaints = []
    try:
        ignored = set( unused_class_parameters(class_data_dir, root_name, among=DERIVATIVE_FLAGS) )
    except FileNotFoundError as err:
        return ['%s' % (err)]
    try:
        recorded = class_run_parameters(class_data_dir, root_name)
    except FileNotFoundError as err:
        return ['%s' % (err)]

    expected = {'delensing derivatives': 'yes', 'output_derivatives': 'yes',
                'derivative type': str(derivative_type)}
    if include_unlensed:
        expected['calculate_derviaties_wrt_unlensed'] = 'yes'
        expected['unlensed derivative type'] = str(derivative_type)

    for flag in sorted(expected):
        if flag in ignored:
            complaints.append('%s was passed but CLASS did not recognize it' % (flag))
        elif flag not in recorded:
            complaints.append('%s reached neither the read nor the ignored list' % (flag))
        elif recorded[flag] != expected[flag]:
            complaints.append('%s was read as %r, expected %r'
                              % (flag, recorded[flag], expected[flag]))

    return complaints


# CLASS runs

def check_noise_for_class(cmb_noise, lmax_calc):
    r"""
    Check that an effective-noise dictionary is acceptable to the FisherLens CLASS wrapper.

    Parameters
    ----------
    cmb_noise : dict
        Effective noise from :func:`build_effective_noise`, with multipoles under ``'el'``.
    lmax_calc : int
        Highest multipole CLASS will be asked for.

    Returns
    -------
    int
        The multipole the noise was required to reach, for reporting.

    Raises
    ------
    ValueError
        If the multipoles do not start at zero or do not extend far enough for the length check inside ``classWrapTools``.
    """

    #Note: Both conditions are tested here rather than left to ``classWrapTools``, whose own length check reports failure by evaluating a bare name that does not exist and so raises ``NameError`` with no explanation.
    #That check also runs before its ``extraParams`` override is applied, so it always demands multipoles up to ``lmax + 2000`` however ``delta_l_max`` is set; see :data:`GUARD_DELTA_L_MAX`.

    el_noise = np.asarray(cmb_noise['el'])
    if len(el_noise) == 0 or el_noise[0] != 0:
        raise ValueError('The supplied noise must start at multipole zero so that its row index '
                         'matches the multipole once CLASS reads it back; see '
                         'build_effective_noise.')
    required = int(lmax_calc) + GUARD_DELTA_L_MAX
    if el_noise[-1] < required:
        raise ValueError('The supplied noise reaches only ell = %d, but classWrapTools requires '
                         'ell >= %d for lmax_calc = %d. That bound is hardcoded and is tested '
                         'before any delta_l_max override takes effect, so pad the noise to it; '
                         'build_effective_noise does this with its guard_pad argument.'
                         % (el_noise[-1], required, lmax_calc))

    return required


def deflection_noise_for_class(reconstruction, estimator='MV', from_multipole_zero=True):
    r"""
    Re-index a reconstruction-noise curve for supply back to CLASS_delens.

    Parameters
    ----------
    reconstruction : dict or array_like
        Either the reconstruction output of :func:`run_class`, from which ``estimator`` is taken or the noise array itself indexed by multipole from two upward, as CLASS returns it.
    estimator : str, optional
        Which estimator to take when a dictionary is given. Default is ``'MV'``; see the note under :func:`lensing_noise_curves` about what CLASS puts there.
    from_multipole_zero : bool, optional
        Prepend two entries so that the array is indexed by multipole from zero rather than from two. Default is ``True``; see the Notes for why that is not merely cosmetic and pass ``False`` to reproduce the behavior of the earlier pipeline.

    Returns
    -------
    ndarray
        Noise array in the requested convention. Two entries longer than the input when ``from_multipole_zero`` is set.

    Raises
    ------
    ValueError
        If ``estimator`` is absent from the supplied dictionary or if the array is not one dimensional and non-empty.

    Warns
    -----
    UserWarning
        If the supplied curve is not everywhere positive and finite since CLASS divides by it.
    """

    #Note: CLASS_delens reads an externally supplied reconstruction noise by *row* rather than by multipole and converts it from the deflection to the potential convention using that row counter in place of the multipole:
    #    ple->pk_rcn_ext[index_l] = plrn[index_l]/ll/(1.+ll);     /* ll = (double)index_l */
    #The multipole column read from the file is stored in ``l_rcn_ext`` and never used.
    #The result is then consumed as ``pk_rcn_ext[l]``, indexed by multipole and starting at ``l = 0``.
    #Since ``classWrapTools`` writes the file starting at multipole two, an array supplied in the convention CLASS itself returns is displaced by two multipoles and both the conversion factor and the lookup are then taken at the wrong multipole.
    #Indexing from zero instead makes the row counter coincide with the multipole, so the conversion and the lookup are both correct.
    #The first entry is read but cannot be made meaningful: the conversion divides by ``ll*(1.+ll)`` with ``ll`` zero, giving infinity for a positive input and an indeterminate value for a zero one. The prepended entries are therefore set to the lowest supplied value rather than to zero so that the result is an infinity rather than a nan. CLASS's own source marks this region as needing improvement.
    #Verified against CLASS_delens ``f27ff1b``. If that submodule pointer is advanced, this is one of the places to re-check.

    if isinstance(reconstruction, dict):
        if estimator not in reconstruction:
            raise ValueError("The reconstruction output has no '%s' entry; its estimators are %s"
                             % (estimator, sorted(key for key in reconstruction if key != 'l')))
        noise = np.asarray(reconstruction[estimator], dtype=float)
    else:
        noise = np.asarray(reconstruction, dtype=float)

    if noise.ndim != 1 or len(noise) == 0:
        raise ValueError('The reconstruction noise must be a non-empty one-dimensional array, got '
                         'shape %s' % (noise.shape,))
    usable = np.isfinite(noise) & (noise > 0.)
    if not usable.all():
        warnings.warn('The reconstruction noise is not everywhere positive and finite: %d of %d '
                      'entries are not. CLASS divides by it.' % ((~usable).sum(), len(noise)),
                      stacklevel=2)
    if not from_multipole_zero:
        return noise

    #a positive value here yields an infinity rather than a nan from the division by zero at ll = 0
    filler = noise[usable].min() if usable.any() else 1.

    return np.concatenate([np.repeat(filler, 2), noise])


def reconstruction_mask(lmin_lensing=2, lmax_T_lensing=3000, lmax_P_lensing=5000):
    r"""
    Assemble the multipole mask for the lensing reconstruction.

    Parameters
    ----------
    lmin_lensing : int, optional
        Lowest multipole entering the reconstruction, in every spectrum. Default is 2.
    lmax_T_lensing : int, optional
        Highest temperature multipole entering the reconstruction. Default is 3000, which mitigates reconstruction biases from non-Gaussian foregrounds.
    lmax_P_lensing : int, optional
        Highest polarization multipole entering the reconstruction, applied to both :math:`E` and :math:`B`. Default is 5000.

    Returns
    -------
    dict
        Mask in the form ``classWrapTools`` expects, with the six entries ``lmin_T``, ``lmax_T``, ``lmin_E``, ``lmax_E``, ``lmin_B`` and ``lmax_B``.

    Raises
    ------
    ValueError
        If a lower bound is below 2 or is not below its matching upper bound.

    Notes
    -----
        The mask is read by ``classWrapTools.class_generate_data`` of FisherLens only on the iterative branch, that is only when no external deflection noise is supplied since otherwise CLASS performs no reconstruction of its own.
    """

    if lmin_lensing < 2:
        raise ValueError('lmin_lensing must be at least 2, got %s' % (lmin_lensing))
    for label, lmax in [('lmax_T_lensing', lmax_T_lensing), ('lmax_P_lensing', lmax_P_lensing)]:
        if lmax <= lmin_lensing:
            raise ValueError('%s must exceed lmin_lensing = %s, got %s' % (label, lmin_lensing, lmax))

    #A fresh dictionary is returned on every call because ``fisherTools.createEllsToUseDict``, which handles the analogous structure for the Fisher sum, edits the dictionary it is given in place; keeping to fresh copies avoids the same pattern biting here.

    return {'lmin_T': int(lmin_lensing), 'lmax_T': int(lmax_T_lensing),
            'lmin_E': int(lmin_lensing), 'lmax_E': int(lmax_P_lensing),
            'lmin_B': int(lmin_lensing), 'lmax_B': int(lmax_P_lensing)}


def run_class(cosmo_fid,
              lmax_calc,
              cmb_noise=None,
              deflection_noise=None,
              recon_mask=None,
              delensing=None,
              extra_params=None,
              root_name=None,
              paths=None,
              use_temp_scratch=False,
              class_wrap_tools=None,
              verbose=True
              ):
    r"""
    Run CLASS_delens once and return the CMB spectra and the lensing-reconstruction noise.

    Parameters
    ----------
    cosmo_fid : dict
        Cosmological parameters, validated by :func:`validate_cosmology`.
    lmax_calc : int
        Highest multipole CLASS is asked for, passed through as its ``l_max_scalars``.
        Choose it above the multipole range of the Fisher sum so that the lensing convolution is accurate there.
    cmb_noise : dict, optional
        Effective noise from :func:`build_effective_noise`, with multipoles under ``'el'`` starting at zero. It is re-keyed for FisherLens by :func:`noise_for_class` on the way in.
        Default is ``None``, which makes CLASS use its own idealized white-noise model instead.
    deflection_noise : array_like, optional
        Lensing-reconstruction noise :math:`N_\ell^{dd}`, indexed by multipole from zero, as :func:`deflection_noise_for_class` returns.
        Default is ``None``, which is the case of interest: CLASS then reconstructs the lensing potential iteratively and returns the noise it achieved.
        Supplying an array instead switches off the iteration and delenses with the given noise, on a code path whose indexing is discussed under :func:`deflection_noise_for_class`.
    recon_mask : dict, optional
        Reconstruction mask from :func:`reconstruction_mask`. Default is ``None``, which lets CLASS use its own bounds and so places no cut on the temperature multipoles.
    delensing : str, optional
        One of :data:`DELENSING_MODES`, passed through to CLASS unchanged. Default is ``None``, which leaves whatever ``classWrapTools`` chose, namely iterative delensing when no external deflection noise is supplied and a single pass when one is.
    extra_params : dict, optional
        Additional CLASS settings, overriding those ``classWrapTools`` sets. Default is ``None``.
    root_name : str, optional
        Base name for the CLASS input and output files. Default is ``None``, which generates a unique one.
    paths : dict, optional
        Resolved paths from :func:`resolve_paths`. Default is ``None``, which resolves them.
    use_temp_scratch : bool, optional
        Place the CLASS files in a temporary directory and delete it afterward, leaving nothing behind. Default is ``False``, which uses the scratch directory in ``paths``.
    class_wrap_tools : module, optional
        An already-imported ``classWrapTools``. Default is ``None``, which imports it.
    verbose : bool, optional
        Report the settings and where the files went. Default is ``True``.

    Returns
    -------
    powers : dict
        Spectra keyed ``'unlensed'``, ``'lensed'``, ``'delensed'`` and ``'lensing'``, as returned by ``classWrapTools.class_generate_data`` of FisherLens.
        The first three hold ``'l'`` and ``'cl_XY'`` in μK², the last holds ``'cl_phiphi'``, ``'cl_dd'`` and ``'cl_kk'``.
    reconstruction : dict or None
        Lensing-reconstruction noise :math:`N_\ell^{dd}` for every estimator, keyed ``'l'``, ``'MV'``, ``'TT'``, ``'TE'``, ``'EE'``, ``'BB'``, ``'EB'`` and ``'TB'``.
        ``None`` when ``deflection_noise`` was supplied since CLASS then performs no reconstruction.

    Raises
    ------
    ValueError
        If ``lmax_calc`` is not positive, or if ``cmb_noise`` does not start at multipole zero or does not extend far enough for the length check inside ``classWrapTools``.
        Also if ``delensing`` is not one of :data:`DELENSING_MODES`, if iterative delensing is requested alongside an external deflection noise or if the mode is given both directly and through ``extra_params``.
    RuntimeError
        If CLASS did not produce the output files that ``classWrapTools`` goes on to read.
    """

    #Note: Three code choices of ``classWrapTools`` are worked around here.
    #(i) Its own length check on the supplied noise runs before its ``extraParams`` override is applied, so it always demands multipoles up to ``lmax + 2000`` however ``delta_l_max`` is set, and reports the failure by evaluating a bare name that does not exist, giving a ``NameError`` rather than an explanation. The same condition is therefore tested first here, with a message that says what to do.
    #(ii) It launches CLASS with :func:`os.system` and discards the return code, so a failed run surfaces later as an unrelated error from reading a file that was never written. The output files CLASS is expected to have produced are therefore checked before ``classWrapTools.class_generate_data`` of FisherLens is allowed to read them, which is possible because a unique ``root_name`` guarantees that anything found was written by this call.
    #(iii) It sets ``overwrite_root``, so runs sharing a ``root_name`` and scratch directory overwrite the output of one another. Unique names are generated by default for that reason and ``use_temp_scratch`` isolates a run completely. CLASS prints a good deal to the terminal, which cannot be suppressed from here because it is a separate process writing to the inherited output stream.
    #(iv) It couples the delensing mode to the source of the reconstruction noise, choosing iterative delensing with an internal reconstruction when no deflection noise is supplied and a single pass with an external one otherwise, so an internal reconstruction followed by a single delensing pass is not reachable through its interface. It is reachable here because its ``extraParams`` override is applied after that choice is made and nothing in between depends on it and because CLASS writes the reconstruction noise whenever the reconstruction is internal, independently of how many delensing passes follow.

    if lmax_calc <= 0:
        raise ValueError('lmax_calc must be positive, got %s' % (lmax_calc))
    lmax_calc = int(lmax_calc)

    class_settings = {} if extra_params is None else dict(extra_params)
    if delensing is not None:
        #the parameter file returns a numpy string, which compares equal but prints awkwardly
        delensing = str(delensing)
        if delensing not in DELENSING_MODES:
            raise ValueError('delensing must be one of %s, got %r'
                             % (', '.join(DELENSING_MODES), delensing))
        if delensing == 'iterative' and deflection_noise is not None:
            raise ValueError('Iterative delensing needs an internal reconstruction to iterate on, '
                             'but an external deflection noise was supplied. Leave '
                             'deflection_noise as None to iterate or ask for a single pass.')
        if 'delensing' in class_settings:
            raise ValueError('The delensing mode was given both as an argument (%r) and through '
                             'extra_params (%r). Give it once.'
                             % (delensing, class_settings['delensing']))
        class_settings['delensing'] = delensing

    if class_wrap_tools is None:
        class_wrap_tools = import_fisherlens(None if paths is None else paths['fisherlens_dir'])[0]
    if paths is None:
        paths = resolve_paths()
    if root_name is None:
        root_name = 'draft_%s' % (uuid.uuid4().hex[:10])

    if cmb_noise is not None:
        check_noise_for_class(cmb_noise, lmax_calc)

    scratch = None
    if use_temp_scratch:
        scratch = tempfile.mkdtemp(prefix='draft_class_')
        class_data_dir = os.path.join(scratch, '')
    else:
        class_data_dir = paths['class_data_dir']
        os.makedirs(class_data_dir, exist_ok=True)

    if verbose:
        print('CLASS run %s' % (root_name))
        print('  lmax_calc %d, scratch %s' % (lmax_calc, class_data_dir))
        print('  noise %s, deflection noise %s, reconstruction mask %s'
              % ('supplied' if cmb_noise is not None else 'idealized',
                 'supplied' if deflection_noise is not None else 'reconstructed internally',
                 'supplied' if recon_mask is not None else 'CLASS defaults'))
        print('  delensing %s' % (class_settings['delensing'] if 'delensing' in class_settings
                                  else '%s (chosen by classWrapTools)'
                                  % ('yes' if deflection_noise is not None else 'iterative')))

    try:
        try:
            output = class_wrap_tools.class_generate_data(
                cosmo_fid,
                rootName=root_name,
                cmbNoise=None if cmb_noise is None else noise_for_class(cmb_noise),
                deflectionNoise=deflection_noise,
                reconstructionMask=recon_mask,
                lmax=lmax_calc,
                classExecDir=paths['class_exec_dir'],
                classDataDir=class_data_dir,
                extraParams=class_settings,
                outputAllReconstructions=deflection_noise is None
                )
        except Exception as err:
            missing = missing_class_outputs(class_data_dir, root_name, deflection_noise is None)
            if missing:
                raise RuntimeError('CLASS did not write %s under %soutput/. Its own messages are '
                                   'above; a failure there is not reported by classWrapTools of '
                                   'FisherLens, which discards the return code of the command it '
                                   'launches.' % (', '.join(missing), class_data_dir)) from err
            raise

        powers, reconstruction = output
        if deflection_noise is not None:
            reconstruction = None
    finally:
        if scratch is not None:
            shutil.rmtree(scratch, ignore_errors=True)
            if verbose:
                print('  temporary scratch removed')

    return powers, reconstruction


def lensing_derivatives(cosmo_fid,
                        lmax_calc,
                        derivative_type='delensed',
                        cmb_noise=None,
                        deflection_noise=None,
                        recon_mask=None,
                        include_unlensed=True,
                        extra_params=None,
                        root_name=None,
                        paths=None,
                        class_wrap_tools=None,
                        verbose=True
                        ):
    r"""
    Run CLASS_delens once for the derivative matrices of the non-Gaussian covariance.

    Parameters
    ----------
    cosmo_fid : dict
        Cosmological parameters, validated by :func:`validate_cosmology`.
    lmax_calc : int
        Highest multipole CLASS is asked for.
    derivative_type : str, optional
        One of :data:`DERIVATIVE_TYPES`. Default is ``'delensed'``.
    cmb_noise : dict, optional
        Effective noise from :func:`build_effective_noise`. Required for the delensed derivatives and not used for the lensed ones. Default is ``None``.
    deflection_noise : array_like, optional
        Lensing-reconstruction noise. Required for the delensed derivatives. Default is ``None``.
    recon_mask : dict, optional
        Reconstruction cuts from :func:`reconstruction_mask`. Default is ``None``, which leaves the CLASS defaults.
    include_unlensed : bool, optional
        Whether to ask for the derivatives with respect to the unlensed spectra as well. Default is ``True``.
    extra_params : dict, optional
        Further CLASS parameters. Default is ``None``.
    root_name : str, optional
        Base name of the run. Default is ``None``, which generates a unique one.
    paths : dict, optional
        Resolved paths from :func:`resolve_paths`. Default is ``None``.
    class_wrap_tools : module, optional
        An already imported ``classWrapTools``. Default is ``None``.
    verbose : bool, optional
        Whether to report what is being run. Default is ``True``.

    Returns
    -------
    lensing_derivs : dict
        ``dCldCLd``, the derivatives with respect to the lensing-deflection power, keyed as :data:`LENSING_DERIVATIVE_FILES` is.
    unlensed_derivs : dict or None
        ``dCldCLu``, the derivatives with respect to the unlensed spectra, or ``None`` when ``include_unlensed`` is false.

    Raises
    ------
    ValueError
        If ``derivative_type`` is not one of :data:`DERIVATIVE_TYPES` or if the delensed derivatives are asked for without both noise arrays.
    RuntimeError
        If CLASS did not write the derivative files or wrote them without consuming the flags that request them.

    Warns
    -----
    UserWarning
        If noise is supplied for the lensed derivatives, which do not depend on it.

    Notes
    -----
    This is a separate CLASS run from the one :func:`run_class` performs and much heavier: it produces a matrix in the spectrum multipole against the deflection multipole for each spectrum rather than one number per multipole. A run returns only the derivative matrices and discards the spectra it computed on the way, so it cannot double as a fiducial run.

    The lensed derivatives do not depend on the noise, so one run serves every configuration and its result may be shared. The delensed derivatives do depend on it, through both the effective noise and the reconstruction noise, so they are per configuration. That asymmetry is what makes caching worthwhile for the first and not for the second.

    ``calculate_derviaties_wrt_unlensed`` is checked against the parameters CLASS reports as unused since a name CLASS does not recognize is ignored silently. See :func:`unused_class_parameters`.
    """

    derivative_type = str(derivative_type)
    if derivative_type not in DERIVATIVE_TYPES:
        raise ValueError('derivative_type must be one of %s, got %r. There is no unlensed '
                         'derivative: the unlensed spectra do not depend on the lensing potential.'
                         % (', '.join(DERIVATIVE_TYPES), derivative_type))
    if lmax_calc <= 0:
        raise ValueError('lmax_calc must be positive, got %s' % (lmax_calc))
    lmax_calc = int(lmax_calc)

    if derivative_type == 'delensed' and (cmb_noise is None or deflection_noise is None):
        raise ValueError('The delensed derivatives depend on both the effective noise and the '
                         'reconstruction noise, so cmb_noise and deflection_noise are both needed. '
                         'Without them, the result would not describe any configuration.')
    if derivative_type == 'lensed' and (cmb_noise is not None or deflection_noise is not None):
        warnings.warn('The lensed derivatives do not depend on the noise, so supplying it makes '
                      'this run configuration specific for no gain. One run without it serves '
                      'every configuration.', stacklevel=2)

    if cmb_noise is not None:
        check_noise_for_class(cmb_noise, lmax_calc)
    if class_wrap_tools is None:
        class_wrap_tools = import_fisherlens(None if paths is None else paths['fisherlens_dir'])[0]
    if paths is None:
        paths = resolve_paths()
    if root_name is None:
        root_name = 'draft_derivs_%s' % (uuid.uuid4().hex[:10])
    class_data_dir = paths['class_data_dir']
    os.makedirs(class_data_dir, exist_ok=True)

    if verbose:
        print('CLASS derivative run %s' % (root_name))
        print('  %s derivatives, lmax_calc %d, scratch %s'
              % (derivative_type, lmax_calc, class_data_dir))
        print('  %d matrices of (%d, %d), %.1f GB in memory'
              % (len(LENSING_DERIVATIVE_FILES)
                 + (len(UNLENSED_DERIVATIVE_FILES) if include_unlensed else 0),
                 lmax_calc+1, lmax_calc+1,
                 ng_derivative_bytes(lmax_calc, include_unlensed)/1e9))

    try:
        output = class_wrap_tools.class_generate_data(
            cosmo_fid,
            rootName=root_name,
            cmbNoise=None if cmb_noise is None else noise_for_class(cmb_noise),
            deflectionNoise=deflection_noise,
            reconstructionMask=recon_mask,
            lmax=lmax_calc,
            classExecDir=paths['class_exec_dir'],
            classDataDir=class_data_dir,
            extraParams={} if extra_params is None else dict(extra_params),
            calculateDerivatives=derivative_type,
            includeUnlensedSpectraDerivatives=include_unlensed
            )
    except Exception as err:
        missing = missing_class_derivatives(class_data_dir, root_name, derivative_type,
                                            include_unlensed)
        dropped = []
        try:
            dropped = unused_class_parameters(class_data_dir, root_name, among=DERIVATIVE_FLAGS)
        except FileNotFoundError:
            pass
        if dropped:
            raise RuntimeError('CLASS ignored %s, so it was never asked for these derivatives and '
                               'did not write %s. The name is spelled by classWrapTools of '
                               'FisherLens as CLASS is expected to spell it, so a mismatch here '
                               'means the submodule pin and the CLASS_delens version disagree.'
                               % (', '.join(dropped), ', '.join(missing) or 'them')) from err
        if missing:
            raise RuntimeError('CLASS did not write %s under %soutput/. Its own messages are '
                               'above; a failure there is not reported by classWrapTools of '
                               'FisherLens, which discards the return code of the command it '
                               'launches.' % (', '.join(missing), class_data_dir)) from err
        raise

    complaints = derivative_flags_consumed(class_data_dir, root_name, derivative_type,
                                           include_unlensed)
    if complaints:
        warnings.warn('CLASS did not take the derivative request as given: %s. The matrices below '
                      'may not be what was asked for.' % ('; '.join(complaints)), stacklevel=2)

    #The return is a pair when the unlensed derivatives are requested and a bare dictionary otherwise (FisherLens 652eaec).
    if include_unlensed:
        lensing_derivs, unlensed_derivs = output
    else:
        lensing_derivs, unlensed_derivs = output, None

    return lensing_derivs, unlensed_derivs


def lensing_noise_curves(powers, reconstruction, lmax=None, settings=None, extra=None):
    r"""
    Assemble the lensing-reconstruction noise product.

    Parameters
    ----------
    powers : dict
        Spectra from :func:`run_class`, whose ``'lensing'`` entry supplies the fiducial lensing spectra.
    reconstruction : dict
        Reconstruction noise from :func:`run_class`, keyed by estimator.
    lmax : int, optional
        Highest multipole to keep. Default is ``None``, which keeps everything CLASS returned.
    settings : dict, optional
        Forecasting settings, from which the reconstruction cuts and the multipole ranges are recorded alongside the curves so the product is self describing. Default is ``None``.
    extra : dict, optional
        Further entries to record, such as the configuration name and the effective noise the reconstruction was performed with. Default is ``None``.

    Returns
    -------
    dict
        Multipoles under ``'el'``; the fiducial lensing spectra under ``'cl_phiphi'``, ``'cl_dd'`` and ``'cl_kk'``; the reconstruction noise per estimator under ``'nl_dd'``, keyed ``'MV'``, ``'TT'``, ``'TE'``, ``'EE'``, ``'BB'``, ``'EB'`` and ``'TB'``; the estimator that ``'MV'`` actually holds under ``'mv_estimator'``; and the recorded settings.

    Raises
    ------
    ValueError
        If ``reconstruction`` is ``None``, which is the case when CLASS was given an external deflection noise and so performed no reconstruction of its own.

    Notes
    -----
    Everything is in deflection convention, :math:`C_\ell^{dd} = \ell(\ell+1) C_\ell^{\phi\phi}`, which is what CLASS_delens returns and what the Fisher sum uses.
    Convert to the convergence convention with

    .. math::

        N_\ell^{\kappa\kappa} = \frac{\ell(\ell+1)}{4} N_\ell^{dd}\, ,

    and likewise for the signal; ``cl_kk`` is included so that a signal-to-noise ratio can be formed without repeating the conversion.

    Two properties of the reconstruction noise are worth knowing before using it:

    * The :math:`BB` estimator is not a number at any multipole when there are no primordial tensor modes, because the unlensed :math:`B`-mode power then vanishes and the estimator has no response.
      This does not corrupt the minimum-variance entry: CLASS_delens combines inverse variances and an estimator with no response contributes nothing to such a sum, so its infinite noise drops out before the value written to the file acquires its indeterminate form.
      No promise is made here that every entry of ``'nl_dd'`` is finite.

    * CLASS_delens builds its minimum-variance entry on one of three branches, depending on whether the full estimator covariance, only its diagonal or an iterative :math:`EB` reconstruction is requested, and on the last of those it sets the entry equal to a single estimator rather than combining any.
      The label alone therefore does not say what the entry holds, so it is compared against each individual estimator and the result recorded under ``'mv_estimator'``, which reads ``'MV'`` when the entry is a genuine combination.
      ``classWrapTools`` of FisherLens requests the diagonal form and the entry was checked to reproduce the inverse-variance sum of the five estimators.
    """

    #Note: This format departs from the lensing-noise products shipped in ``products/``, which were written by an earlier pipeline.
    #Those store the noise as ``Nl_TT``, ``Nl_EB`` and so on at the top level, in the convergence convention, as ``complex128`` arrays with a zero imaginary part, and use ``Nl_ET`` where CLASS_delens says ``TE``.
    #They carry one indeterminate entry, at multipole zero and are finite everywhere above it including below their own recorded ``lmin``, so nothing in the file marks where the reconstruction actually began.
    #Here the arrays are real, free of ``nan``, grouped under a single ``'nl_dd'`` entry and named as CLASS_delens names them.
    #Readers of the older products therefore need updating rather than reusing.

    if reconstruction is None:
        raise ValueError('There is no reconstruction noise to record. CLASS only reconstructs the '
                         'lensing potential when it is not given an external deflection noise, so '
                         'run_class must be called with deflection_noise left as None.')

    el = np.asarray(powers['lensing']['l'])
    keep = slice(None) if lmax is None else slice(None, int(np.searchsorted(el, lmax, side='right')))

    product = {'el': el[keep]}
    for spec in ('cl_phiphi', 'cl_dd', 'cl_kk'):
        product[spec] = np.asarray(powers['lensing'][spec], dtype=float)[keep]

    estimators = [key for key in reconstruction if key != 'l']
    product['nl_dd'] = {key: np.asarray(reconstruction[key], dtype=float)[keep]
                        for key in estimators}
    product['estimators'] = sorted(estimators)

    #CLASS labels one entry minimum variance but, where the reconstruction is iterative, sets it
    #equal to a single estimator instead of combining them. Record which one it turned out to be
    #rather than trusting the label; see the Notes.
    product['mv_estimator'] = 'MV'
    if 'MV' in product['nl_dd']:
        for name in sorted(key for key in estimators if key != 'MV'):
            if np.array_equal( product['nl_dd']['MV'], product['nl_dd'][name] ):
                product['mv_estimator'] = name
                break

    recorded = ('lmin_lensing', 'lmax_T_lensing', 'lmax_P_lensing', 'lmin', 'lmax', 'lmax_TT',
                'lmin_dd', 'lmax_dd', 'lbuffer', 'lmin_ilc', 'sat_expname', 'include_sat')
    if settings is not None:
        for key in recorded:
            if key in settings:
                product[key] = settings[key]
    if extra is not None:
        product.update(extra)

    return product


def parameter_derivatives(cosmo_fid,
                          step_sizes,
                          vary_params,
                          lmax_calc,
                          cmb_noise,
                          deflection_noise,
                          pol_combs=None,
                          spectrum_types=None,
                          extra_params=None,
                          root_name=None,
                          paths=None,
                          use_temp_scratch=False,
                          fisher_tools=None,
                          verbose=True
                          ):
    r"""
    Differentiate the CMB and lensing spectra with respect to each varied parameter.

    Two CLASS runs are performed per varied parameter, at the fiducial value displaced by plus and minus its step size, and the central difference taken.
    The reconstruction noise is held fixed at the value the fiducial run achieved, so these runs do not repeat the iterative reconstruction and are much cheaper than it.

    Parameters
    ----------
    cosmo_fid : dict
        Fiducial parameter values, validated by :func:`validate_cosmology`.
    step_sizes : dict
        Step sizes, validated by :func:`validate_cosmology`. Must include ``mnu`` whenever the cosmology does, even when it is held fixed; see the Notes there.
    vary_params : list of str
        Parameters to differentiate. Anything in ``cosmo_fid`` but absent here reaches CLASS and is held fixed.
    lmax_calc : int
        Highest multipole CLASS is asked for, as in :func:`run_class`.
    cmb_noise : dict
        Effective noise from :func:`build_effective_noise`, with multipoles under ``'el'`` starting at zero.
    deflection_noise : array_like
        Reconstruction noise from :func:`deflection_noise_for_class`, indexed by multipole from zero.
    pol_combs : list of str, optional
        Spectra to differentiate. Default is ``None``, which uses ``['cl_TT', 'cl_TE', 'cl_EE', 'cl_dd']``, matching the multipole windows of the Fisher sum.
    spectrum_types : list of str, optional
        Which spectra to differentiate. Default is ``None``, which takes ``['unlensed', 'lensed', 'delensed', 'lensing']``, matching the default of ``getPowerDerivWithParams`` itself.
        ``'lensing'`` is required whenever ``'cl_dd'`` appears in ``pol_combs`` since that is where the lensing derivative is taken from.
    extra_params : dict, optional
        Additional CLASS settings. Default is ``None``.
    root_name : str, optional
        Base name for the CLASS files. Default is ``None``, which generates a unique one. See the Notes on why uniqueness matters more here than elsewhere.
    paths : dict, optional
        Resolved paths from :func:`resolve_paths`. Default is ``None``, which resolves them.
    use_temp_scratch : bool, optional
        Place the CLASS files in a temporary directory and delete it afterward. Default is ``False``.
    fisher_tools : module, optional
        An already-imported ``fisherTools``. Default is ``None``, which imports it.
    verbose : bool, optional
        Report the settings and the parameters being differentiated. Default is ``True``.

    Returns
    -------
    dict
        Derivatives indexed by parameter, then by spectrum type, then by spectrum, as ``fisherTools.getPowerDerivWithParams`` returns them. The lensing derivative is under ``[parameter]['lensing']['cl_dd']``.

    Raises
    ------
    ValueError
        If ``vary_params`` is empty, if a varied parameter has no step size, if ``'cl_dd'`` is requested without ``'lensing'`` among the spectrum types, if the returned derivatives are missing an entry that was asked for or propagated from :func:`check_noise_for_class` and :func:`deflection_noise_for_class`.

    Notes
    -----
    ``fisherTools.getPowerDerivWithParams`` gives every CLASS run it launches the same root name, and ``classWrapTools`` sets ``overwrite_root``, so the runs overwrite one another's output in turn.
    That is harmless in sequence since each is read before the next begins, but two forecasts sharing a root name and scratch directory would destroy each other's intermediate files.
    A unique root name is therefore generated by default and ``use_temp_scratch`` isolates a forecast completely.

    There is deliberately no delensing-mode argument here, unlike :func:`run_class`. These runs are given a fixed reconstruction noise, so ``classWrapTools`` selects a single delensing pass for them and there is nothing to iterate on; the mode is settled by the fiducial run that produced that noise.
    That means a forecast built on an iterative fiducial run mixes conventions, holding the iteratively obtained reconstruction noise fixed while delensing each parameter step in a single pass. Holding it fixed is the usual approximation, but the reconstruction noise in principle responds to a parameter step, so the derivatives of the delensed spectra are approximate in that respect.

    A one-sided derivative is taken only for a short hardcoded list of parameters and for ``mnu`` when its step size exceeds its fiducial value. Any other parameter whose step size would carry it past a physical boundary is differentiated two-sidedly across that boundary, which is why :func:`validate_cosmology` warns about it.
    """

    if pol_combs is None:
        pol_combs = ['cl_TT', 'cl_TE', 'cl_EE', 'cl_dd']
    if spectrum_types is None:
        spectrum_types = ['unlensed', 'lensed', 'delensed', 'lensing']
    pol_combs = list(pol_combs)
    spectrum_types = list(spectrum_types)
    if 'cl_dd' in pol_combs and 'lensing' not in spectrum_types:
        raise ValueError("'cl_dd' is to be differentiated but 'lensing' is absent from "
                         'spectrum_types. getPowerDerivWithParams takes the lensing derivative from '
                         "the 'lensing' spectra, and if that entry is missing it prints a message and "
                         'returns without the derivative rather than failing, so the omission would '
                         'only surface later as a missing key.')
    vary_params = list(vary_params)
    if not vary_params:
        raise ValueError('vary_params is empty, so there is nothing to differentiate')
    missing = [name for name in vary_params if name not in step_sizes]
    if missing:
        raise ValueError('No step size for %s, which are to be differentiated'
                         % (', '.join(sorted(missing))))

    lmax_calc = int(lmax_calc)
    check_noise_for_class(cmb_noise, lmax_calc)
    deflection_noise = np.asarray(deflection_noise, dtype=float)
    if deflection_noise.ndim != 1 or len(deflection_noise) == 0:
        raise ValueError('deflection_noise must be a non-empty one-dimensional array, got shape %s'
                         % (deflection_noise.shape,))

    if fisher_tools is None:
        fisher_tools = import_fisherlens(None if paths is None else paths['fisherlens_dir'])[1]
    if paths is None:
        paths = resolve_paths()
    if root_name is None:
        root_name = 'draft_deriv_%s' % (uuid.uuid4().hex[:10])

    scratch = None
    if use_temp_scratch:
        scratch = tempfile.mkdtemp(prefix='draft_class_')
        class_data_dir = os.path.join(scratch, '')
    else:
        class_data_dir = paths['class_data_dir']
        os.makedirs(class_data_dir, exist_ok=True)

    if verbose:
        print('derivatives for %d parameter(s): %s' % (len(vary_params), ', '.join(vary_params)))
        print('  %d CLASS runs, root %s, scratch %s'
              % (2 * len(vary_params), root_name, class_data_dir))
        print('  spectra %s, types %s' % (', '.join(pol_combs), ', '.join(spectrum_types)))

    #``useClass`` is passed explicitly because it defaults to ``False`` in ``getPowerDerivWithParams``, which would route the calculation through CAMB rather than CLASS_delens and so produce no delensed spectra at all.
    try:
        param_derivs = fisher_tools.getPowerDerivWithParams(
            cosmoFid=cosmo_fid,
            stepSizes=step_sizes,
            polCombs=pol_combs,
            cmbNoiseSpectraK=noise_for_class(cmb_noise),
            deflectionNoisesK=deflection_noise,
            spectrumTypes=spectrum_types,
            paramsToDifferentiate=vary_params,
            useClass=True,
            lmax=lmax_calc,
            fileNameBase=root_name,
            classExecDir=paths['class_exec_dir'],
            classDataDir=class_data_dir,
            extraParams={} if extra_params is None else dict(extra_params)
            )
    finally:
        if scratch is not None:
            shutil.rmtree(scratch, ignore_errors=True)
            if verbose:
                print('  temporary scratch removed')

    #getPowerDerivWithParams reports a missing spectrum type by printing and carrying on, so check
    #that everything asked for came back rather than letting the gap surface as a missing key later
    for name in vary_params:
        if name not in param_derivs:
            raise ValueError('No derivative was returned for %s' % (name))
        for spectrum in spectrum_types:
            expected = ['cl_dd'] if spectrum == 'lensing' else [c for c in pol_combs if c != 'cl_dd']
            absent = [c for c in expected if c not in param_derivs[name].get(spectrum, {})]
            if absent:
                raise ValueError('The derivatives for %s are missing %s under %r. Look above for a '
                                 'message from getPowerDerivWithParams, which reports such gaps by '
                                 'printing rather than by failing.'
                                 % (name, ', '.join(absent), spectrum))

    return param_derivs
