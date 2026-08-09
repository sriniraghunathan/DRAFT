r"""
Internal linear combination (ILC) of multi-frequency CMB data supporting the residual calculation in :mod:`get_ilc_residuals`.

The multi-frequency CMB observations are combined into a single estimate of the CMB signal using a harmonic-space internal linear combination.
The ILC exploits the known frequency independence of the CMB signal (in thermodynamic temperature units), together with the distinct frequency dependence of the foregrounds, to construct an optimal estimate by linearly combining the maps from all available frequency channels, while down-weighting modes that are dominated by foreground emission or instrumental noise.

The module is grouped into four sections:

* Covariance: :func:`get_analytic_covariance`
* Frequency response: :func:`get_cib_freq_dep`, :func:`get_radio_freq_dep`, :func:`get_acap`
* Covariance assembly: :func:`get_teb_spec_combination`, :func:`create_clmat`, :func:`get_clinv`, :func:`corr_from_cov`
* ILC residuals: :func:`residual_power`

:func:`get_analytic_covariance`, :func:`residual_power` and :func:`create_clmat` are called by :mod:`get_ilc_residuals`.
:func:`corr_from_cov` is a standalone utility that is currently not used elsewhere in this repository.
The remaining routines are internal helpers.

In harmonic space, the ILC estimate of the CMB signal at multipole :math:`\ell` is

.. math::

    S_\ell = \sum_{i=1}^{N_f} w_\ell^i\, M_\ell^i\, ,

where :math:`M_\ell^i` is the harmonic-space representation of the observed map in frequency channel :math:`i` and the sum runs over all :math:`N_f` frequency channels.
The multipole-dependent weights :math:`w_\ell^i` are determined by minimizing the total variance of :math:`S_\ell` subject to the constraint that the CMB signal is preserved, with the option to additionally null or suppress the frequency response of one or more other sky components.
The spectral energy distributions of all components included in the solution are collected into a matrix :math:`\mathcal{F} = [A_s, B_{f_1}, C_{f_2}, \ldots]_{N_f \times N_c}`, where :math:`N_c` is the number of components.
The first column :math:`A_s` is the target signal and :math:`B_{f_1}, C_{f_2}, \ldots` are the :math:`N_c - 1` undesired (foreground) components.
For the CMB, the target signal is frequency independent, so :math:`A_s = A_\mathrm{CMB} = [1, 1, \ldots, 1]`.
The constrained-ILC weights then take the form

.. math::

    w_\ell = \mathbf{C}_\ell^{-1} \mathcal{F} \left( \mathcal{F}^\mathrm{T} \mathbf{C}_\ell^{-1} \mathcal{F} \right)^{\!-1} U\, ,

where the constraint vector :math:`U = [1, 0, \ldots, 0]` has length :math:`N_c` and selects :math:`A_\mathrm{CMB}` as the preserved component.
The :math:`N_f \times N_f` covariance matrix :math:`\mathbf{C}_\ell` of the noise and foreground power at multipole :math:`\ell` is constructed from the auto- and cross-frequency power spectra of the instrumental noise and of the galactic and extragalactic foreground components, after computing them on the masked sky.

In the limit where no additional components are nulled, :math:`\mathcal{F}` reduces to the :math:`N_f \times 1` frequency response vector :math:`A_\mathrm{CMB}` and the weights simplify to the standard minimum-variance form,

.. math::

    w_\ell^\mathrm{MV} = \frac{\mathbf{C}_\ell^{-1} A_\mathrm{CMB}}{A_\mathrm{CMB}^\mathrm{T} \mathbf{C}_\ell^{-1} A_\mathrm{CMB}}\, .

The corresponding residual power spectrum after component separation is

.. math::

    C_\ell^\mathrm{res} = \left( A_\mathrm{CMB}^\mathrm{T} \mathbf{C}_\ell^{-1} A_\mathrm{CMB} \right)^{-1}\, ,

which represents the effective noise that includes residuals from both the experimental noise and the foreground signals in the component-separated CMB map.

The framework defines four forms of the ILC:

* Standard minimum-variance ILC minimizes the variance subject to unit response for the required component. This is the form used for the forecasts in the accompanying paper (arXiv:2608.XXXXX).
* Constrained ILC additionally nulls the frequency response of one or more other sky components, selected with the ``null_comp`` argument of :func:`residual_power`.
* Partial ILC suppresses the targeted component by artificially rescaling the covariance of that component rather than imposing an explicit nulling constraint, selected with the ``cl_multiplier_dic`` argument of :func:`get_analytic_covariance`. This is useful when the frequency response of an unwanted component is not known a priori or cannot be captured by a single spectral energy distribution, as is the case for the cosmic infrared background.
* Cross-ILC constructs two separate constrained-ILC maps that null different sky components and estimates the CMB power spectrum from their cross-correlation. *[To be ported.]*

:func:`get_analytic_covariance` builds :math:`\mathbf{C}_\ell`, :func:`get_acap` builds :math:`A_\mathrm{CMB}` and :func:`residual_power` combines them.
Power spectra are in units of μK² and frequencies are in GHz.
Band ordering follows the ``freqarr`` supplied by the caller: :func:`get_acap`, :func:`create_clmat` and the returned weights are all indexed in that order, so ``freqarr``, ``freqcalib_fac`` and the covariance must be consistent.
Component names follow the spelling used in the returned foreground dictionary, namely ``cib`` for the cosmic infrared background and ``galdust``/``galsync`` for the galactic terms.
Various parameter values and data locations are taken from ``param_dict``, which is populated from ``params.ini``.
"""

import re
import warnings

import numpy as np

import foregrounds as fg
import misc


# Covariance

def get_analytic_covariance(
        param_dict,
        freqarr,
        el=None,
        nl_dic=None,
        bl_dic=None,
        ignore_fg=[],
        which_spec='TT',
        pol_frac_per_cent_dust=0.02,
        pol_frac_per_cent_radio=0.03,
        pol_frac_per_cent_tsz=0.,
        pol_frac_per_cent_ksz=0.,
        include_gal=0,
        cib_corr_coeffs=None,
        null_highfreq_radio=1,
        reduce_radio_power_150=None,
        reduce_tsz_power=None,
        reduce_cib_power=None,
        remove_cib_decorr=0,
        cl_multiplier_dic=None,
        return_fg_spectra=True,
        force_cl_dic=None
        ):
    r"""
    Signal and noise covariance :math:`\mathbf{C}_\ell` across a set of frequency bands.

    The covariance combines the CMB, the extragalactic and galactic foregrounds, and the instrumental noise:

    .. math::

        C_{\ell, ij} = C_\ell^\mathrm{CMB} + C_\ell^\mathrm{kSZ} + C_{\ell, ij}^\mathrm{tSZ}
        + C_{\ell, ij}^\mathrm{radio} + C_{\ell, ij}^\mathrm{CIB}
        + C_{\ell, ij}^{\mathrm{tSZ} \times \mathrm{CIB}} + C_{\ell, ij}^\mathrm{gal} + N_{\ell, ij}\, ,

    where the CMB and kSZ terms are frequency independent in thermodynamic units.
    Individual terms can be excluded with ``ignore_fg``, which is how the component being recovered is kept out of the covariance, or rescaled with ``cl_multiplier_dic`` for a partial ILC.


    Parameters
    ----------
    param_dict : dict
        Parameter dictionary, as returned by :func:`misc.get_param_dict`.
        Must contain ``freq0``, ``fg_model``, ``spec_index_rg``, ``spec_index_dg_po``, ``spec_index_dg_clus`` and ``Tcib``, and the galactic file names when ``include_gal`` is set.
    freqarr : array_like
        Frequency bands in GHz. The returned dictionary is keyed by pairs drawn from this list and the ordering is the one used throughout the ILC.
    el : array_like, optional
        Multipole moments :math:`\ell`, which must start at zero.
        Default is ``None``, in which case the multipoles of the SPT foreground templates are used.
    nl_dic : dict, optional
        Noise power spectra, keyed either by band or by band pair. A band pair that is absent is taken to have zero cross-band noise.
        Each entry must cover at least the requested multipoles. Default is ``None``, i.e. a noise-free covariance.
    bl_dic : dict, optional
        Beam transfer functions :math:`B_\ell` keyed by band, used only to correct the galactic templates, which carry the S4 beams. Default is ``None``.
    ignore_fg : list, optional
        Components to leave out of the covariance, drawn from ``cmb``, ``tsz``, ``y``, ``ksz``, ``radio``, ``cib``, ``noise`` and ``tsz_cib``.
        Requesting ``cmb`` also removes the kSZ. Default is an empty list.
    which_spec : str, optional
        Spectrum to compute, one of ``'TT'``, ``'EE'`` or ``'TE'``. Default is ``'TT'``.
    pol_frac_per_cent_dust : float, optional
        Polarization fraction of the CIB. Default is 0.02.
    pol_frac_per_cent_radio : float, optional
        Polarization fraction of the radio sources. Default is 0.03.
    pol_frac_per_cent_tsz : float, optional
        Polarization fraction of the tSZ. Default is 0.
    pol_frac_per_cent_ksz : float, optional
        Polarization fraction of the kSZ. Default is 0.
    include_gal : bool, optional
        Include galactic dust and synchrotron from the PySM simulations. Default is 0.
    cib_corr_coeffs : dict, optional
        CIB decorrelation coefficients keyed by band pair, applied to the cross-band CIB power.
        Every off-diagonal pair must be present if the dictionary is supplied at all.
        Default is ``None``, i.e. fully correlated.
    null_highfreq_radio : bool, optional
        Null the radio power above the high-frequency threshold of :func:`foregrounds.get_cl_radio`. Default is 1.
    reduce_radio_power_150 : float, optional
        Rescale the radio power at 150 GHz by this factor. Default is ``None``.
    reduce_tsz_power : float, optional
        Rescale the tSZ power by this factor. Default is ``None``.
    reduce_cib_power : float, optional
        Rescale the CIB power by this factor. Default is ``None``.
    remove_cib_decorr : bool, optional
        Replace the cross-band CIB with the geometric mean of the two auto-spectra, which removes the decorrelation built into the template. Overrides ``cib_corr_coeffs``. Default is 0.
    cl_multiplier_dic : dict, optional
        Factors by which individual components are rescaled in the covariance for a partial ILC.
        Keys are drawn from ``cmb``, ``ksz``, ``tsz``, ``radio``, ``cib``, ``tsz_cib``, ``galdust``, ``galsync`` and ``noise``. Default is ``None``.
    return_fg_spectra : bool, optional
        Also return the individual component spectra. Default is ``True``.
    force_cl_dic : dict, optional
        Spectra supplied by the caller which override the internally computed ones.
        Recognized keys are ``cmb``, ``ksz``, ``y``, ``tsz``, ``cib``, ``tsz_cib`` or ``cib_tsz``, and ``rad`` or ``radio``. Default is ``None``.

    Returns
    -------
    el : ndarray
        Multipole moments over which the covariance is defined.
    cl_dic : dict
        Total covariance, keyed by band pair.
    fg_cl_dic : dict
        Spectrum of each component, keyed by component name and then by band pair.
        Only returned when ``return_fg_spectra`` is set.

    Raises
    ------
    ValueError
         If ``freqarr`` is empty or contains a repeated band, if ``el`` is empty or does not start at zero, if an entry of ``nl_dic`` is shorter than ``el``, if ``ignore_fg`` or ``cl_multiplier_dic`` contains an unrecognized component or if ``cib_corr_coeffs`` is supplied without an entry for every band pair.

    Notes
    -----
    ``cib`` denotes the cosmic infrared background throughout, whereas ``galdust`` denotes galactic dust.
    Note that :func:`foregrounds.get_cl_galactic` uses ``dust`` for the latter.

    Noise that has been amplified by beam deconvolution is clipped at :math:`5 \times 10^4`\ μK² to limit the dynamic range of the covariance.
    The galactic terms are currently not clipped, so at multipoles where the beam correction of a large-beam channel diverges the covariance can still become badly conditioned.
    """

    #ignore_fg = foreground terms that must be ignored
    #debug=False  #unused
    #'dust' renamed to 'cib'
    #possible_ignore_fg = ['cmb', 'tsz', 'y', 'ksz', 'radio', 'dust', 'noise', 'tsz_cib']
    possible_ignore_fg = ['cmb', 'tsz', 'y', 'ksz', 'radio', 'cib', 'noise', 'tsz_cib']
    if len(ignore_fg) > 0:
        ignore_fg = list(ignore_fg)
        if 'cmb' in ignore_fg: ignore_fg.append('ksz')
        if not all( [ currfg in possible_ignore_fg for currfg in ignore_fg] ):
            bad_fg = [currfg for currfg in ignore_fg if currfg not in possible_ignore_fg]
            raise ValueError('ignore_fg entries must be one of %s, got %s' % (possible_ignore_fg, bad_fg))
        ignore_fg = np.unique(ignore_fg)

    el_, cl_cmb = fg.get_foreground_power_spt('CMB', freq1=param_dict['freq0'], freq2=param_dict['freq0'])
    if el is None: el = np.copy(el_)
    if len(freqarr) == 0:
        raise ValueError('freqarr must not be empty')
    if len(set(freqarr)) != len(freqarr):
        raise ValueError('freqarr must not contain repeated bands, got %s' % (list(freqarr)))
    if len(el) == 0:
        raise ValueError('el must not be empty')
    if min(el) != 0:
        raise ValueError('el must start at 0, got el[0] = %s' % (min(el)))
    if nl_dic is not None:
        for nl_key in nl_dic:
            if len(nl_dic[nl_key]) < len(el):
                raise ValueError('nl_dic[%s] has %s entries, fewer than the %s multipoles in el' % (nl_key, len(nl_dic[nl_key]), len(el)))
    cl_cmb = np.interp(el, el_, cl_cmb)
    el_, cl_ksz = fg.get_foreground_power_spt('kSZ', freq1=param_dict['freq0'], freq2=param_dict['freq0'])
    cl_ksz = np.interp(el, el_, cl_ksz)

    if which_spec == 'EE':
        cl_ksz = cl_ksz * pol_frac_per_cent_ksz**2.
    if which_spec == 'TE':
        cl_ksz = cl_ksz * 0.

    if cl_multiplier_dic is not None:
        if 'cmb' in cl_multiplier_dic:
            cl_cmb = np.copy(cl_cmb) * cl_multiplier_dic['cmb']
        if 'ksz' in cl_multiplier_dic:
            cl_ksz = np.copy(cl_ksz) * cl_multiplier_dic['ksz']

    #20220428 - force cl if force_cl_dic is supplied
    if force_cl_dic is None: force_cl_dic = {}
    if 'cmb' in force_cl_dic:
        cl_cmb = np.interp(el, np.arange(len(force_cl_dic['cmb'])), force_cl_dic['cmb'])
    if 'ksz' in force_cl_dic:
        cl_ksz = np.interp(el, np.arange(len(force_cl_dic['ksz'])), force_cl_dic['ksz'])

    if cl_multiplier_dic is not None:
        possible_multipliers = ['cmb', 'ksz', 'tsz', 'radio', 'cib', 'tsz_cib', 'galdust', 'galsync', 'noise']
        bad_mult = [kk for kk in cl_multiplier_dic if kk not in possible_multipliers]
        if len(bad_mult) > 0:
            raise ValueError('cl_multiplier_dic keys must be one of %s, got %s' % (possible_multipliers, bad_mult))

    cl_dic = {}
    cl_ori = np.zeros(len(el))

    if return_fg_spectra:
        fg_cl_dic = {}
        for curr_comp in ['cmb', 'ksz', 'tsz', 'radio', 'cib', 'tsz_cib', 'noise']:
            fg_cl_dic[curr_comp] = {}
        if include_gal:
            fg_cl_dic['galdust'] = {}
            fg_cl_dic['galsync'] = {}

    for freq1 in freqarr:
        for freq2 in freqarr:
            #get tsz
            el_, cl_tsz = fg.get_cl_tsz(freq1, freq2, freq0=param_dict['freq0'], fg_model=param_dict['fg_model'], reduce_tsz_power=reduce_tsz_power)
            if which_spec == 'EE':
                cl_tsz = cl_tsz * pol_frac_per_cent_tsz**2.
            elif which_spec == 'TE':
                cl_tsz = cl_tsz * 0.
            cl_tsz = np.interp(el, el_, cl_tsz)

            #get radio
            el_, cl_radio = fg.get_cl_radio(freq1, freq2, freq0=param_dict['freq0'], fg_model=param_dict['fg_model'], spec_index_rg=param_dict['spec_index_rg'], null_highfreq_radio=null_highfreq_radio, reduce_radio_power_150=reduce_radio_power_150)
            if which_spec == 'EE':
                cl_radio = cl_radio * pol_frac_per_cent_radio**2.
            elif which_spec == 'TE':
                cl_radio = cl_radio * 0.
            cl_radio = np.interp(el, el_, cl_radio)

            #get CIB
            el_, cl_dg_po, cl_dg_clus = fg.get_cl_dust_cib(freq1, freq2, freq0=param_dict['freq0'], fg_model=param_dict['fg_model'], spec_index_dg_po=param_dict['spec_index_dg_po'], spec_index_dg_clus=param_dict['spec_index_dg_clus'], Tcib=param_dict['Tcib'], reduce_cib_power=reduce_cib_power)
            cl_dust = cl_dg_po + cl_dg_clus
            if remove_cib_decorr:
                el_, cl_dg_po1, cl_dg_clus1 = fg.get_cl_dust_cib(freq1, freq1, freq0=param_dict['freq0'], fg_model=param_dict['fg_model'], spec_index_dg_po=param_dict['spec_index_dg_po'], spec_index_dg_clus=param_dict['spec_index_dg_clus'], Tcib=param_dict['Tcib'], reduce_cib_power=reduce_cib_power)
                el_, cl_dg_po2, cl_dg_clus2 = fg.get_cl_dust_cib(freq2, freq2, freq0=param_dict['freq0'], fg_model=param_dict['fg_model'], spec_index_dg_po=param_dict['spec_index_dg_po'], spec_index_dg_clus=param_dict['spec_index_dg_clus'], Tcib=param_dict['Tcib'], reduce_cib_power=reduce_cib_power)
                cl_dust1 = cl_dg_po1 + cl_dg_clus1
                cl_dust2 = cl_dg_po2 + cl_dg_clus2
                cl_dust = np.sqrt( cl_dust1 * cl_dust2 )
                cib_corr_coeffs = None

            cl_dust[np.isnan(cl_dust)] = 0.
            if which_spec == 'EE':
                cl_dust = cl_dust * pol_frac_per_cent_dust**2.
            elif which_spec == 'TE':
                cl_dust = cl_dust * 0.

            #include CIB decorrelation if available
            if cib_corr_coeffs is not None:
                if freq1 == freq2:
                    corr_coeff = 1.
                else:
                    if (freq1, freq2) in cib_corr_coeffs:
                        corr_coeff = cib_corr_coeffs[(freq1, freq2)]
                    elif (freq2, freq1) in cib_corr_coeffs:
                        corr_coeff = cib_corr_coeffs[(freq2, freq1)]
                    else:
                        raise ValueError('cib_corr_coeffs has no entry for the (%s, %s) band pair' % (freq1, freq2))
                cl_dust *= corr_coeff
            cl_dust = np.interp(el, el_, cl_dust)

            #get tSZ x CIB
            el_, cl_tsz_cib = fg.get_cl_tsz_cib(freq1, freq2, freq0=param_dict['freq0'], fg_model=param_dict['fg_model'], spec_index_dg_po=param_dict['spec_index_dg_po'], spec_index_dg_clus=param_dict['spec_index_dg_clus'], Tcib=param_dict['Tcib'], reduce_tsz_power=reduce_tsz_power)
            if which_spec == 'EE' or which_spec == 'TE':
                cl_tsz_cib = cl_tsz_cib * 0.
            cl_tsz_cib = np.interp(el, el_, cl_tsz_cib)

            #galaxy
            if include_gal:# and not pol: #get galactic dust and sync
                el_, cl_gal_dust = fg.get_cl_galactic(param_dict, 'dust', freq1, freq2, which_spec, el=el, bl_dic=bl_dic)
                el_, cl_gal_sync = fg.get_cl_galactic(param_dict, 'sync', freq1, freq2, which_spec, el=el, bl_dic=bl_dic)

            cl = np.copy( cl_ori )

            #20220428 - force cl if force_cl_dic is supplied
            if 'y' in force_cl_dic:
                cl_y_force = force_cl_dic['y']
                tsz_fac_freq1 = fg.compton_y_to_delta_Tcmb(freq1)
                tsz_fac_freq2 = fg.compton_y_to_delta_Tcmb(freq2)
                tsz_fac = tsz_fac_freq1 * tsz_fac_freq2
                cl_tsz_force = cl_y_force * tsz_fac
                cl_tsz = np.interp(el, np.arange(len(cl_tsz_force)), cl_tsz_force)
            if 'tsz' in force_cl_dic:
                cl_tsz_force = force_cl_dic['tsz'][(freq1, freq2)]
                cl_tsz = np.interp(el, np.arange(len(cl_tsz_force)), cl_tsz_force)
            if 'cib' in force_cl_dic:
                if (freq1, freq2) in force_cl_dic['cib']:
                    cl_dust_force = force_cl_dic['cib'][(freq1, freq2)]
                else:
                    cl_dust_force = force_cl_dic['cib'][(freq2, freq1)]
                cl_dust = np.interp(el, np.arange(len(cl_dust_force)), cl_dust_force)
            if 'tsz_cib' in force_cl_dic or 'cib_tsz' in force_cl_dic:
                if 'tsz_cib' in force_cl_dic:
                    cl_tsz_cib_force = force_cl_dic['tsz_cib'][(freq1, freq2)]
                else:
                    cl_tsz_cib_force = force_cl_dic['cib_tsz'][(freq1, freq2)]
                cl_tsz_cib = np.interp(el, np.arange(len(cl_tsz_cib_force)), cl_tsz_cib_force)
            if 'rad' in force_cl_dic or 'radio' in force_cl_dic:
                if 'rad' in force_cl_dic:
                    cl_radio_force = force_cl_dic['rad'][(freq1, freq2)]
                else:
                    cl_radio_force = force_cl_dic['radio'][(freq1, freq2)]
                cl_radio = np.interp(el, np.arange(len(cl_radio_force)), cl_radio_force)

            if cl_multiplier_dic is not None:
                if 'tsz' in cl_multiplier_dic:
                    cl_tsz = np.copy(cl_tsz) * cl_multiplier_dic['tsz']
                if 'radio' in cl_multiplier_dic:
                    cl_radio = np.copy(cl_radio) * cl_multiplier_dic['radio']
                #if 'dust' in cl_multiplier_dic:
                if 'cib' in cl_multiplier_dic:
                    cl_dust = np.copy(cl_dust) * cl_multiplier_dic['cib']
                if 'tsz_cib' in cl_multiplier_dic:
                    cl_tsz_cib = np.copy(cl_tsz_cib) * cl_multiplier_dic['tsz_cib']
                if include_gal:
                    if 'galdust' in cl_multiplier_dic:
                        cl_gal_dust = np.copy(cl_gal_dust) * cl_multiplier_dic['galdust']
                    if 'galsync' in cl_multiplier_dic:
                        cl_gal_sync = np.copy(cl_gal_sync) * cl_multiplier_dic['galsync']

            if 'cmb' not in ignore_fg:
                if len(cl_cmb) < len(cl):
                    cl_cmb = np.interp(el, np.arange(len(cl_cmb)), cl_cmb)
                cl = cl + np.copy(cl_cmb[el])
            if 'ksz' not in ignore_fg:
                if len(cl_ksz) < len(cl):
                    cl_ksz = np.interp(el, np.arange(len(cl_ksz)), cl_ksz)
                cl = cl + cl_ksz[el]
            if 'tsz' not in ignore_fg and 'y' not in ignore_fg:
                if len(cl_tsz) < len(cl):
                    cl_tsz = np.interp(el, np.arange(len(cl_tsz)), cl_tsz)
                cl = cl + cl_tsz[el]
            if 'radio' not in ignore_fg:
                if len(cl_radio) < len(cl):
                    cl_radio = np.interp(el, np.arange(len(cl_radio)), cl_radio)
                cl = cl + cl_radio[el]
            #if 'dust' not in ignore_fg:
            if 'cib' not in ignore_fg:
                if len(cl_dust) < len(cl):
                    cl_dust = np.interp(el, np.arange(len(cl_dust)), cl_dust)
                cl = cl + cl_dust[el]

            #20220503 - add tszxcib if either tsz or cib is included.
            add_cl_tsz_cib = True
            #if ('cib' in ignore_fg and 'tsz' in ignore_fg) or 'tsz_cib' in ignore_fg or 'cib_tsz' in ignore_fg:
            if 'tsz_cib' in ignore_fg or 'cib_tsz' in ignore_fg or 'cib_y' in ignore_fg or 'y_cib' in ignore_fg:
                add_cl_tsz_cib = False
            if add_cl_tsz_cib: #'cib' not in ignore_fg and 'tsz' not in ignore_fg and 'tsz_cib' not in ignore_fg:
                if len(cl_tsz_cib) < len(cl):
                    cl_tsz_cib = np.interp(el, np.arange(len(cl_tsz_cib)), cl_tsz_cib)
                cl = cl + cl_tsz_cib[el]

            if include_gal:# and not pol: #get galactic dust and sync
                cl = cl + cl_gal_dust
                cl = cl + cl_gal_sync

            #el is required to start at 0 above and get_foreground_power_spt already pads to l=0, so commented
            ##make sure cl start from el=0 rather than el=10 which is the default for SPT G15 results
            #lmin = min(el)
            #cl = np.concatenate( (np.zeros(lmin), cl) )

            #noise auto power spectrum
            if nl_dic is not None:
                if (freq1, freq2) in nl_dic:
                    nl = nl_dic[(freq1, freq2)]
                else:
                    nl = nl_dic[freq1]
                    if freq1 != freq2:
                        nl = np.copy(nl) * 0.

                if len(cl) > len(nl):
                    cl = cl[:len(nl)]
                elif len(cl) < len(nl):
                    nl = nl[:len(cl)]

                #remove very large numbers because of beam deconvolution
                ini_nl = np.median(nl[:100])
                end_nl = np.median(nl[-100:])
                if end_nl > ini_nl: #this implies beam deconvolution has made end nl pretty large
                    max_nl_value = 5e4 #some large number
                    #having end_nl pretty large introduces covariance inversion issues
                    badinds = np.where(nl >= max_nl_value)[0]
                    nl = np.copy(nl)
                    nl[badinds] = max_nl_value
                    #print(ini_nl, end_nl)

                el = np.arange(len(cl))

            else:
                nl = np.zeros(len(cl))

            if cl_multiplier_dic is not None:
                if 'noise' in cl_multiplier_dic:
                    nl = np.copy(nl) * cl_multiplier_dic['noise']

            if 'noise' not in ignore_fg:
                if which_spec != 'TE':
                    cl = cl + np.copy(nl)

            if return_fg_spectra:
                #these are symmetric under an exchange of the two bands
                curr_fg_cl_arr = [('cmb', cl_cmb), ('ksz', cl_ksz), ('tsz', cl_tsz), ('radio', cl_radio), ('cib', cl_dust), ('noise', nl)]
                if include_gal:
                    curr_fg_cl_arr = curr_fg_cl_arr + [('galdust', cl_gal_dust), ('galsync', cl_gal_sync)]
                for curr_comp, curr_cl in curr_fg_cl_arr:
                    fg_cl_dic[curr_comp][(freq1, freq2)] = fg_cl_dic[curr_comp][(freq2, freq1)] = curr_cl
                #cl_tsz_cib is not symmetric under an exchange of the two bands, so only one ordering
                fg_cl_dic['tsz_cib'][(freq1, freq2)] = cl_tsz_cib

            cl[np.isnan(cl)] = 0.
            cl[np.isinf(cl)] = 0.

            #Commented since False anyway
            ##20200516 - adjusting Nl when beam is too large (for 30/40 GHz bands)
            #adjust_for_large_beams = False
            #if adjust_for_large_beams:
            #    beam_tol_for_ilc = 1000. #some large number
            #    bl = bl_dic[freq1]
            #    if 'effective' in bl_dic:
            #        bl_eff = bl_dic['effective']
            #    else:
            #        bl_eff = bl_dic[145]
            #    bad_inds = np.where( (bl_eff/bl > beam_tol_for_ilc) )[0]
            #    print(freq1, freq2, bad_inds)
            #    cl[bad_inds] = 0.

            cl_dic[(freq1, freq2)] = cl

    if return_fg_spectra:
        return el, cl_dic, fg_cl_dic
    else:
        return el, cl_dic


# Frequency response

def get_cib_freq_dep(nu, Tcib=20., beta=1.505):  ##Tcmb=2.7255,  #Tcmb was unused
    r"""
    Frequency dependence of the CIB as a modified blackbody:

    .. math::

        f(\nu) = \nu^\beta \, \frac{B_\nu(T_\mathrm{CIB})}{\mathrm{d}B_\nu/\mathrm{d}T}\, ,

    where the division by :math:`\mathrm{d}B_\nu/\mathrm{d}T` converts from flux to thermodynamic temperature units.

    Parameters
    ----------
    nu : float
        Frequency in GHz.
    Tcib : float, optional
        Dust temperature :math:`T_\mathrm{CIB}` in K. Default is 20.
    beta : float, optional
        Emissivity index :math:`\beta`. Default is 1.505, the Poisson CIB value; the clustered component uses 2.505.

    Returns
    -------
    value : float
        Frequency dependence in thermodynamic temperature units, unnormalized.

    Raises
    ------
    ValueError
        If ``nu`` is not positive, or is at or above :math:`10^4` and so is presumably in Hz rather than GHz.
    """

    misc.check_freqs_in_ghz(nu)
    bnu1 = fg.get_BnuT(nu, temp=Tcib)
    dbdt = fg.get_dB_dT(nu)
    value = ((nu * 1e9)**beta) * bnu1 / dbdt  #the frequency power law is evaluated in Hz

    return value


def get_radio_freq_dep(nu, nu0=150., spec_index_rg=-0.76, null_highfreq_radio=True, highfreq_radio_threshold=230):
    r"""
    Frequency dependence of radio sources as a power law:

    .. math::

        f(\nu) = \frac{\mathrm{d}B_{\nu_0}/\mathrm{d}T}{\mathrm{d}B_\nu/\mathrm{d}T} \left( \frac{\nu}{\nu_0} \right)^{\!\alpha} ,

    with the ratio of blackbody derivatives converting from flux to thermodynamic temperature units.

    Parameters
    ----------
    nu : float
        Frequency in GHz.
    nu0 : float, optional
        Reference frequency in GHz. Default is 150.
    spec_index_rg : float, optional
        Radio spectral index :math:`\alpha`. Default is -0.76.
    null_highfreq_radio : bool, optional
        Return zero above ``highfreq_radio_threshold``, where the radio contribution is negligible. Default is ``True``.
    highfreq_radio_threshold : float, optional
        Threshold in GHz for the above. Default is 230.

    Returns
    -------
    value : float
        Frequency dependence in thermodynamic temperature units, normalized to unity at ``nu0`` up to the unit conversion.
    Raises
    ------
    ValueError
        If ``nu`` or ``nu0`` is not positive, or is at or above :math:`10^4` and so is presumably in Hz rather than GHz.
    """

    misc.check_freqs_in_ghz(nu, nu0, highfreq_radio_threshold)

    nr = fg.get_dB_dT(nu0)
    dr = fg.get_dB_dT(nu)
    epsilon_nu1_nu0 = nr/dr
    scaling = (nu/nu0)**spec_index_rg
    value = epsilon_nu1_nu0 * scaling

    if null_highfreq_radio and (nu > highfreq_radio_threshold):
        value = 0.

    return value


def get_acap(freqarr, final_comp='cmb', freqcalib_fac=None, nspecs=1, spec_index_rg=-0.76):
    r"""
    Frequency response :math:`A_s` of a sky component across a set of bands.

    The response is the factor by which the component appears in each band relative to the reference, in the same thermodynamic units as the maps, multiplied by any calibration factor:

    .. math::

        A_i = f(\nu_i) \, g_i\, ,

    with :math:`f` the spectral dependence of the component, :math:`\nu_i` the frequency of band :math:`i` and :math:`g_i` the calibration factor.

    Parameters
    ----------
    freqarr : array_like
        Frequency bands in GHz. The response is returned in this order, which must match the band ordering used to build the covariance.
    final_comp : str, optional
        Component whose response is required. Default is ``'cmb'``. One of

        * ``'cmb'``, unit response in every band,
        * ``'tsz'`` or ``'y'``, the thermal SZ spectral function from :func:`foregrounds.compton_y_to_delta_Tcmb`,
        * ``'cib'`` or ``'cibpo'``, the Poisson CIB from :func:`get_cib_freq_dep`,
        * ``'cibclus'``, the clustered CIB, i.e. :func:`get_cib_freq_dep` with :math:`\beta = 2.505`,
        * ``'misc_cib_tcibXX_betaYY'``, the CIB with dust temperature ``XX`` and emissivity index ``YY``,
        * ``'radio'``, radio sources from :func:`get_radio_freq_dep`.

    freqcalib_fac : array_like, optional
        Multiplicative gain of each band, in the same order as ``freqarr``, by which the assumed response is scaled.
        Default is ``None``, which is equivalent to unity in every band.
    nspecs : int, optional
        Number of map components entering the ILC jointly, as returned by :func:`get_teb_spec_combination`.
        Default is 1, i.e. temperature or polarization alone.
    spec_index_rg : float, optional
        Radio spectral index, used only when ``final_comp`` is ``'radio'``. Default is -0.76.

    Returns
    -------
    acap : ndarray
        Frequency response, of shape ``(nspecs * len(freqarr), nspecs)``.
        For a joint temperature and polarization ILC only the CMB has a nonzero polarization response.

    Raises
    ------
    ValueError
        If no frequency response is implemented for ``final_comp``.
        This includes ``ksz``, ``noise`` and ``tsz_cib``, which :func:`residual_power` accepts as components but which have no response here.
         Also if ``freqcalib_fac`` does not have exactly one entry per band.
    """

    #``'misc_radio_alphaYY'``, radio with spectral index ```YY`` is not implemented, so removed from the docstring
    nc = len(freqarr)

    if freqcalib_fac is None:
        freqcalib_fac = np.ones(nc)
    elif np.shape(freqcalib_fac) != (nc,):
        raise ValueError('freqcalib_fac must have one entry per band, got shape %s for %s bands' % (np.shape(freqcalib_fac), nc))

    if final_comp.lower() == 'cmb':
        freqscale_fac = np.ones(nc)

    elif final_comp.lower() == 'tsz' or final_comp.lower() == 'y':
        freqscale_fac = []
        #for freq in sorted( freqarr ):
        for freq in freqarr:
            freqscale_fac.append( fg.compton_y_to_delta_Tcmb(freq) )

        freqscale_fac = np.asarray( freqscale_fac )

    #The tSZ-CIB cross-term does not have its own SED, so has no single-frequency response. It used to be caught here and left freqscale_fac unset; it now falls through to the else below.
    #elif final_comp.lower() == 'tsz_cib' or final_comp.lower() == 'cib_tsz':
    #    pass

    elif final_comp.lower() == 'cib' or final_comp.lower() == 'cibpo':
        freqscale_fac = []
        #for freq in sorted( freqarr ):
        for freq in freqarr:
            freqscale_fac.append( get_cib_freq_dep(freq) )

        freqscale_fac = np.asarray( freqscale_fac )
        freqscale_fac /= np.max(freqscale_fac)

    elif final_comp.lower() == 'cibclus':
        freqscale_fac = []
        #for freq in sorted( freqarr ):
        for freq in freqarr:
            freqscale_fac.append( get_cib_freq_dep(freq, beta=2.505) )

        freqscale_fac = np.asarray( freqscale_fac )
        freqscale_fac /= np.max(freqscale_fac)

    elif final_comp.lower().find('misc_cib') > -1:
        #default values
        misc_tcib = 20.
        misc_beta = 1.505
        tcib_tmp = re.findall(r'tcib\d*\.?\d+', final_comp.lower())
        if len(tcib_tmp) > 0:
            tcib_tmp = tcib_tmp[0]
            misc_tcib = float(tcib_tmp.replace('tcib', ''))

        beta_tmp = re.findall(r'beta\d*\.?\d+', final_comp.lower())
        if len(beta_tmp) > 0:
            beta_tmp = beta_tmp[0]
            misc_beta = float(beta_tmp.replace('beta', ''))

        #freqarr = [30, 44, 70, 100, 150, 217, 353, 545]
        freqscale_fac = []
        #for freq in sorted( freqarr ):
        for freq in freqarr:
            freqscale_fac.append( get_cib_freq_dep(freq, Tcib=misc_tcib, beta=misc_beta) )

        freqscale_fac = np.asarray( freqscale_fac )
        freqscale_fac /= np.max(freqscale_fac)

    elif final_comp.lower() == 'radio':
        freqscale_fac = []
        #for freq in sorted( freqarr ):
        for freq in freqarr:
            freqscale_fac.append( get_radio_freq_dep(freq, spec_index_rg=spec_index_rg) )

        freqscale_fac = np.asarray( freqscale_fac )

    else:
        #'ksz', 'noise' and 'tsz_cib' land here: residual_power accepts them as components, but no frequency response is implemented for them
        raise ValueError('No frequency response is implemented for final_comp = %s' % (final_comp))

    acap = np.zeros(nc) + (freqscale_fac * freqcalib_fac) #assuming CMB is the same and calibrations factors are same for all channel

    if nspecs > 1:
        acap_full = np.zeros( (nspecs, len(acap) * nspecs) )
        acap_full[0, :len(acap)] = acap
        if final_comp.lower() == 'cmb':
            acap_full[1, len(acap):] = acap
        else: #polarization weights are zero for other foregrounds
            acap_full[1, len(acap):] = 0.

        #acap_full = np.asmatrix(acap_full).T
        acap_full = acap_full.T  #should be nspecs*nc x nspecs
        acap = acap_full
    else:
        #acap = np.asmatrix(acap).T
        acap = acap.reshape(-1, 1)  #should be nspecs*nc x nspecs, .T alone would leave a 1d array

    return acap


# Covariance assembly

def get_teb_spec_combination(cl_dic):
    """
    Determine whether the ILC is to be performed on one map component or on several jointly.

    Parameters
    ----------
    cl_dic : dict
        Auto- and cross-spectra of different frequency channels keyed by spectrum name, i.e. ``'TT'``, ``'EE'``, ``'TE'``, ``'BB'``, ``'TB'``, ``'EB'``.

    Returns
    -------
    nspecs : int
        Number of map components entering the ILC, namely 1 for temperature or polarization alone, 2 for temperature and E-modes jointly, and 3 when B-modes are included.
    specs : list of str or None
        The spectrum names in sorted order or ``None`` if the combination of keys is not one of those recognized.

    Notes
    -----
    When :func:`create_clmat` is uaws, only 1 and 2 are usable since B-modes are absent.
    """

    specs = sorted( list( cl_dic.keys() ) )

    if specs == ['TT'] or specs == ['EE']: #only TT is supplied
        nspecs = 1
    elif specs == sorted(['TT', 'EE']) or specs == sorted(['TT', 'EE', 'TE']): #TT/EE/TE are supplied
        nspecs = 2
    elif specs == sorted(['TT', 'EE', 'BB']) or specs == sorted(['TT', 'EE', 'BB', 'TE', 'TB', 'EB']): #TT/EE/BB are supplied
        nspecs = 3
    else:
        nspecs = 1
        specs = None

    return nspecs, specs


def create_clmat(freqarr, elcnt, cl_dic):
    r"""
    Band-band covariance matrix at one multipole.

    For a single map component, the matrix is :math:`\mathbf{C}_\ell` over the requested bands.
    For a joint temperature and polarization ILC, it is the block matrix

    .. math::

        \begin{pmatrix} C^{TT}_{ij} & C^{TE}_{ij} \\ C^{TE}_{ij} & C^{EE}_{ij} \end{pmatrix}

    of size ``nspecs * len(freqarr)`` on a side.


    Parameters
    ----------
    freqarr : array_like
        Frequency bands in GHz. Rows and columns follow this order.
    elcnt : int
        Index into the multipole axis of the spectra, not the multipole itself.
    cl_dic : dict
        Auto- and cross-spectra at all :math:`\ell` keyed by spectrum name and then by band pair.

    Returns
    -------
    clmat : ndarray
        Covariance matrix at the requested multipole index.

    Raises
    ------
    ValueError
        If the keys of ``cl_dic`` are not a recognized combination of spectra.
    """

    nc = len(freqarr)
    nspecs, pspec_arr = get_teb_spec_combination(cl_dic)
    if pspec_arr is None:
        raise ValueError('cl_dic keys %s are not a supported spectrum combination' % (sorted(cl_dic.keys())))
    clmat = np.zeros( (nspecs * nc, nspecs * nc) )

    for pspecind, pspec in enumerate( pspec_arr ):
        curr_cl_dic = cl_dic[pspec]
        if nspecs == 1: #clmat for TT or EE or BB
            for ncnt1, freq1 in enumerate(freqarr):
                for ncnt2, freq2 in enumerate(freqarr):
                    j, i = ncnt2, ncnt1
                    clmat[j, i] = curr_cl_dic[(freq1, freq2)][elcnt]
        else: #joint or separate TT/EE constraints #fix me: include BB for joint constraints.
            if pspec == 'TT':
                for ncnt1, freq1 in enumerate(freqarr):
                    for ncnt2, freq2 in enumerate(freqarr):
                        j, i = ncnt2, ncnt1
                        clmat[j, i] = curr_cl_dic[(freq1, freq2)][elcnt]
            elif pspec == 'EE':
                for ncnt1, freq1 in enumerate(freqarr):
                    for ncnt2, freq2 in enumerate(freqarr):
                        j, i = ncnt2 + nc, ncnt1 + nc
                        clmat[j, i] = curr_cl_dic[(freq1, freq2)][elcnt]
            elif pspec == 'TE':
                for ncnt1, freq1 in enumerate(freqarr):
                    for ncnt2, freq2 in enumerate(freqarr):
                        j, i = ncnt2 + nc, ncnt1
                        clmat[j, i] = curr_cl_dic[(freq1, freq2)][elcnt]
                        clmat[i, j] = curr_cl_dic[(freq1, freq2)][elcnt]

    return clmat


def get_clinv(freqarr, elcnt, cl_dic, return_clmat=False):
    r"""
    Inverse of the band-band covariance matrix at one multipole :math:`\ell`.

    Parameters
    ----------
    freqarr : array_like
        Frequency bands in GHz.
    elcnt : int
        Index into the multipole axis of the spectra, not the multipole :math:`\ell` itself.
    cl_dic : dict
        Auto- and cross-spectra at all :math:`\ell` keyed by spectrum name and then by band pair.
    return_clmat : bool, optional
        Also return the covariance itself, which is convenient for debugging. Default is ``False``.

    Returns
    -------
    clinv : ndarray
        Inverse covariance at the requested multipole index.
    clmat : ndarray
        Covariance at the same index. Only returned when ``return_clmat`` is set.

    Warns
    -----
    UserWarning
        If the covariance contains non-finite entries. :func:`numpy.linalg.pinv` returns zeros for an infinite entry without raising, which would otherwise silently zero the weights at that multipole, and raises :exc:`numpy.linalg.LinAlgError` for a NaN one.

    Notes
    -----
    The Moore-Penrose pseudo-inverse is used, so a rank-deficient covariance, for instance one in which a band carries no information, is handled without raising.
    """

    #clmat = np.asmatrix( create_clmat(freqarr, elcnt, cl_dic) )
    clmat = create_clmat(freqarr, elcnt, cl_dic)
    if not np.all( np.isfinite(clmat) ):
        #np.linalg.pinv silently returns zeros for inf and raises LinAlgError for nan
        warnings.warn('The covariance contains non-finite values. Infinities give zero ILC weights at those multipoles and NaNs make the inversion fail.', stacklevel=2)
    clinv = np.linalg.pinv(clmat)

    if return_clmat:
        return clinv, clmat
    else:
        return clinv


def corr_from_cov(covmat):
    r"""
    Correlation matrix corresponding to a covariance matrix:

    .. math::

        R_{ij} = \frac{C_{ij}}{\sqrt{C_{ii} C_{jj}}}\, .

    Parameters
    ----------
    covmat : array_like
        Covariance matrix, assumed square.

    Returns
    -------
    corrmat : ndarray
        Correlation matrix, of the same shape as ``covmat``.

    Raises
    ------
    ValueError
        If any diagonal entry of ``covmat`` is not positive, which would otherwise give a non-finite correlation.
    """

    diags = np.diag(covmat)
    if np.any(diags <= 0.):
        raise ValueError('covmat must have a positive diagonal, got %s' % (diags))
    corrmat = np.zeros_like(covmat)
    for i in range(covmat.shape[0]):
        for j in range(covmat.shape[0]):
            corrmat[i, j] = covmat[i, j] / np.sqrt( diags[i] * diags[j] )
    #corrmat = covmat / np.outer(np.sqrt(diags), np.sqrt(diags))

    return corrmat


# ILC residuals

def residual_power(
        param_dict,  #unused argument
        freqarr,
        el,
        cl_dic,
        final_comp='cmb',
        null_comp=None,
        spec_index_rg=-0.76,
        freqcalib_fac=None,
        lmin=10,
        return_weights=True
        ):
    r"""
    Residual power and weights of an ILC.

    For a standard ILC, the weights and the residual power follow from the covariance :math:`\mathbf{C}_\ell` and the frequency response :math:`A_s` of the component being recovered:

    .. math::

        w_\ell = \frac{\mathbf{C}_\ell^{-1} A_s}{A_s^\mathrm{T} \mathbf{C}_\ell^{-1} A_s} \, , \qquad C_\ell^\mathrm{res} = \left( A_s^\mathrm{T} \mathbf{C}_\ell^{-1} A_s \right)^{\!-1} .

    When ``null_comp`` is given, a constrained ILC is performed instead.
    Here, :math:`\mathcal{F}` is the matrix whose columns are the response of the component to be recovered followed by those of the components to be nulled, and :math:`U = [1, 0, \ldots, 0]` is the constraint vector,

    .. math::

        w_\ell = \mathbf{C}_\ell^{-1} \mathcal{F} \left( \mathcal{F}^\mathrm{T} \mathbf{C}_\ell^{-1} \mathcal{F} \right)^{\!-1} U \, ,

    which returns the required component with unit response while nulling the others.

    Parameters
    ----------
    param_dict : dict
        Parameter dictionary. Currently unused, but retained so that the call signature matches the other routines.
    freqarr : array_like
        Frequency bands in GHz, in the same order as the covariance and the weights.
    el : array_like
        Multipole moments :math:`\ell`.
    cl_dic : dict
        Auto- and cross-spectra at all :math:`\ell` keyed by spectrum name and then by band pair, as returned by :func:`get_analytic_covariance`.
    final_comp : str, optional
        Component to be recovered. Default is ``'cmb'``.
        Options are 'cmb', 'tsz', 'y', 'ksz', 'radio', 'cib', 'noise' and 'tsz_cib'.
    null_comp : str or list of str, optional
        Component or components to be nulled, which selects a constrained ILC. Default is ``None``.
        Options are the same as for ``final_comp``.
    spec_index_rg : float, optional
        Radio spectral index, passed to :func:`get_acap` and used only for radio components. Default is -0.76.
    freqcalib_fac : array_like, optional
        Multiplicative gain of each band, in the same order as ``freqarr``, by which the assumed response is scaled.
        Default is ``None``, which is equivalent to unity in every band.
    lmin : int, optional
        Multipoles at or below this value are skipped and their residual left at zero. Default is 10.
    return_weights : bool, optional
        Also return the weights. Default is ``True``.

    Returns
    -------
    cl_residual : ndarray
        Residual power, of shape ``(3, len(el))``.
        The rows are the temperature, polarization and cross-terms, with only the first being filled for a single map component.
    weightsarr : ndarray
        ILC weights, of shape ``(nspecs * len(freqarr), nspecs, len(el))``.
        Only returned when ``return_weights`` is set.

    Raises
    ------
    ValueError
        If ``freqarr`` is empty or contains a repeated band, if ``final_comp`` or any entry of ``null_comp`` is not a recognized component, if ``freqcalib_fac`` does not have exactly one entry per band, or if the response of ``final_comp`` lies in the span of the nulled components and so cannot be preserved and nulled at once.

    Notes
    -----
    Non-finite residuals, which arise where the covariance cannot be inverted, are set to zero.

    A constrained ILC is supported only for a single map component.
    Combining ``null_comp`` with a joint temperature and polarization covariance raises ``IndexError`` since the joint constrained formalism is not implemented.
    """

    #assert final_comp in ['cmb', 'tsz', 'y', 'ksz', 'radio', 'dust', 'noise', 'tsz_cib']
    #possible_comps = ['cmb', 'tsz', 'y', 'ksz', 'radio', 'cib', 'noise', 'tsz_cib']
    #Note: 'ksz', 'noise' and 'tsz_cib' were previously accepted here and then failed inside get_acap; 'ksz' should be straight-forward to implement, though
    possible_comps = ['cmb', 'tsz', 'y', 'radio', 'cib']
    if final_comp not in possible_comps:
        raise ValueError('final_comp must be one of %s, got %s' % (possible_comps, final_comp))
    #assert null_comp in [None, 'cmb', 'tsz', 'y', 'ksz', 'radio', 'dust', 'noise', 'tsz_cib']
    if null_comp is not None:
        for curr_null_comp in ([null_comp] if np.ndim(null_comp) == 0 else null_comp):
            if curr_null_comp not in possible_comps:
                raise ValueError('null_comp must be one of %s, got %s' % (possible_comps, curr_null_comp))
    #if freqcalib_fac is not None:
    #    assert len(freqcalib_fac) == len(freqarr)
    if len(freqarr) == 0:
        raise ValueError('freqarr must not be empty')
    if len(set(freqarr)) != len(freqarr):
        raise ValueError('freqarr must not contain repeated bands, got %s' % (list(freqarr)))
    if freqcalib_fac is not None and np.shape(freqcalib_fac) != (len(freqarr),):
        raise ValueError('freqcalib_fac must have one entry per band, got shape %s for %s bands' % (np.shape(freqcalib_fac), len(freqarr)))

    nspecs, pspec_arr = get_teb_spec_combination(cl_dic) #20200527 - teb
    acap = get_acap(freqarr, final_comp=final_comp, freqcalib_fac=freqcalib_fac, nspecs=nspecs, spec_index_rg=spec_index_rg)

    if null_comp is not None:
        total_comp_to_null = 0
        if np.ndim(null_comp) == 0:
            bcap = get_acap(freqarr, final_comp=null_comp, freqcalib_fac=freqcalib_fac, nspecs=nspecs, spec_index_rg=spec_index_rg)
            total_comp_to_null += 1
        else:
            bcap = None
            for curr_null_comp in null_comp:
                curr_bcap = get_acap(freqarr, final_comp=curr_null_comp, freqcalib_fac=freqcalib_fac, nspecs=nspecs, spec_index_rg=spec_index_rg)
                if bcap is None:
                    bcap = np.copy( curr_bcap )
                else:
                    bcap = np.column_stack( (bcap, curr_bcap) )
                total_comp_to_null += 1
            #bcap = np.asmatrix(bcap)

        #A null_comp whose response is collinear with that of final_comp asks for the component to be preserved and nulled at once.
        #Note that 'tsz' and 'y' are aliases, so comparing the names alone would not catch every case.
        bcap_2d = bcap if np.ndim(bcap) > 1 else np.reshape(bcap, (-1, 1))
        if np.linalg.matrix_rank(np.column_stack( (acap, bcap_2d) )) == np.linalg.matrix_rank(bcap_2d):
            raise ValueError('The response of final_comp = %s lies in the span of the nulled components %s, so it cannot be preserved and nulled at the same time' % (final_comp, null_comp))

    nc = len(freqarr)
    weightsarr = np.zeros( (nspecs * nc, nspecs, len( el ) ) )
    cl_residual = np.zeros( (3, len(el)) )

    #cl_residual_tmp = []  #unused

    for elcnt, currel in enumerate(el):
        if currel <= lmin: continue ## or el>=lmax: continue
        #clinv = get_clinv( freqarr, elcnt, cl_dic )
        clinv, clmat = get_clinv( freqarr, elcnt, cl_dic, return_clmat=True )

        if null_comp is None:
            nr = np.dot(clinv, acap)
            dr = np.dot( acap.T, np.dot(clinv, acap) )
            drinv = np.linalg.pinv(dr)
            weight = np.dot(nr, drinv)
        else:
            G = np.column_stack( (acap, bcap) )
            #G = np.asmatrix(G)

            total_comps = G.shape[1]
            ncap = np.zeros( total_comps )#total_comp_to_null + 1 )
            ncap[0] = 1.
            #ncap = np.asmatrix( ncap ).T
            ncap = ncap.reshape(-1, 1)

            nr = np.dot(clinv, G)
            dr = np.dot( G.T, np.dot(clinv, G) )
            drinv = np.dot( np.linalg.pinv(dr), ncap)
            weight = np.dot(nr, drinv)

        #ILC residuals
        if nspecs > 1:
            cl_residual_tt, cl_residual_ee, cl_residual_te = drinv[0, 0], drinv[1, 1], drinv[0, 1]
            cl_residual[:, elcnt] = cl_residual_tt, cl_residual_ee, cl_residual_te
        else:
            cl_residual[0, elcnt] = drinv[0, 0]

        weightsarr[:, :, elcnt] = weight

        #cl_residual_tmp.append( drinv )  #unused

    weightsarr = np.asarray( weightsarr )
    cl_residual = np.asarray( cl_residual )

    cl_residual[np.isinf(cl_residual)] = 0.
    cl_residual[np.isnan(cl_residual)] = 0.

    if return_weights:
        return cl_residual, weightsarr
    else:
        return cl_residual


# Superseded code, kept for now, could be removed later

#coth and compton_y_to_delta_Tcmb are duplicates of the foregrounds.py versions.
#The copies below reference h, k_B and Tcmb, none of which are defined in this module, so every call failed.
#Now instead using fg.compton_y_to_delta_Tcmb.
#def coth(x):
#    """
#    Coth function for tSZ frequency response
#
#    Parameters
#    ----------
#    x: float
#        x = h*nu/(k_B * T_CMB)
#
#    Returns
#    -------
#    coth_x: float
#        Hyperbolic cotangent of x.
#    """
#
#    coth_x = (np.exp(x) + np.exp(-x)) / (np.exp(x) - np.exp(-x))
#
#    return coth_x
#
#
#def compton_y_to_delta_Tcmb(nu):
#    """
#    Get the tSZ frequency response
#
#    Parameters
#    ----------
#    nu: float
#        Frequency in Hz.
#
#    Returns
#    -------
#    g_nu_with_tcmb: float
#        tSZ frequency response in Tcmb units.
#    """
#
#    nu *= 1e9
#    x = h * nu / k_B / Tcmb
#    g_nu = x * coth(x/2.) - 4.
#    g_nu_with_tcmb = Tcmb * np.mean(g_nu)
#
#    return g_nu_with_tcmb
