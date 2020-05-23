"""
simulations/python/foregrounds.py

Functions for setting up simulations of gaussian-
and poisson-distributed foregrounds.

"""
import os
import numpy as np
import healpy as hp
import pickle as pk

# Planck constant in m2 kg / s in G3units
h = 6.62607004e-34
# Boltzmann constant in m2 kg s-2 K-1 in G3units
k_B = 1.38064852e-23
# velocity of light in m/s in G3units
c = 2.997997e8
cmb_temp = 2.725
data_folder = '/Users/sraghunathan/Research/SPTPol/analysis/git/DRAFT/data/'


def band_to_freq_to_query_idl_files(band):
    """
    Function to convert band names to integer values to query idl sav files

    Parameters
    ----------
    band: str
        can be 95GHz, 150GHz, 220GHz

    Returns
    -------
    freq: int
        frequency corresponding to band. 
        This is not the effective frequency.
        To get the effective frequency look into get_spt_effective_frequencies.
    """
    if band not in ['90GHz', '95GHz', '150GHz', '220GHz']:
        raise ValueError("band must be one of 90/95, 150, 220GHz")

    band_dic = {'90GHz': 95, '95GHz': 95, '150GHz': 150, '220GHz': 220}
    return band_dic[band]

def get_spt_effective_frequencies(which_spt, band, component):
    """
    Function to get the SPT-SZ/SPT-SZ+SPTpol/SPT-3G effective freqeuncies for 
    different foreground component.

    Parameters
    ----------
    which_spt: str
        spt experiment for which we need the effective frequencies.
        can be one of: 'george', sptsz', 'reichardt', 'sptszpol', 'stp3g', '3g'
        george/sptsz - Sec 5.7 of https://arxiv.org/pdf/1408.3161.pdf
        reichardt/sptszpol - Sec 5.1 of https://arxiv.org/pdf/2002.06197.pdf
        spt3g/3g - The effective band centres are calculated by Zhaodi 
        (https://pole.uchicago.edu/spt3g/images/Effective_band_center_for_3G.pdf).
    band: str
        95GHz or 150GHz or 220GHz for which we need the effective SPT3G band centre
    component : str
        The foreground component to use. Must be one of
        'tSZ', 'DG-Cl', 'DG-Po', 'DG', 'RG'
    Returns
    -------
    spt_band_centre: float in g3units for the desired foreground component.
    """
    which_spt = which_spt.lower()
    component = component.lower()
    assert which_spt in ['george', 'sptsz', 'reichardt', 'sptszpol', 'spt3g', '3g']
    assert component in ['tsz', 'dg-cl', 'dg-po', 'dg', 'rg']

    if band == '90GHz':
        band = '95GHz'

    if band not in ['95GHz', '150GHz', '220GHz']:
        raise ValueError("band must be one of 90/95, 150, 220GHz")

    if component == 'dg-cl' or component == 'dg-po' or component == 'dg':
        component = 'dg'

    if which_spt == 'george' or which_spt == 'sptsz':
        spt_eff_frequencies = {
            'dg': {'95GHz': 96.9, '150GHz': 153.4, '220GHz': 221.6},
            'rg': {'95GHz': 93.5, '150GHz': 149.5, '220GHz': 215.8},
            'tsz': {'95GHz': 96.6, '150GHz': 152.3, '220GHz': 220.1},
        }
    elif which_spt == 'reichardt' or which_spt == 'sptszpol':
        spt_eff_frequencies = {
            'dg': {'95GHz': 96.9, '150GHz': 153.4, '220GHz': 221.6},
            'rg': {'95GHz': 93.5, '150GHz': 149.5, '220GHz': 215.8},
            'tsz': {'95GHz': 96.6, '150GHz': 152.3, '220GHz': 220.1},
        }
    elif which_spt == 'spt3g' or which_spt == '3g':
        spt_eff_frequencies = {
            'dg': {'95GHz': 95.96, '150GHz': 150.01, '220GHz': 222.76},
            'rg': {'95GHz': 93.52, '150GHz': 145.92, '220GHz': 213.34},
            'tsz': {'95GHz': 95.69, '150GHz': 148.85, '220GHz': 220.15},
        }

    return spt_eff_frequencies[component][band]


def get_spt_foreground_amplitudes_templates(
    component, band='150GHz', fg_model='reichardt', dg_clus_template_id=0, lmax=12000,
):

    """
    Function to get the SPT-SZ/SPTpol (G15/R20) foreground 
    amplitudes of various components and 
    the respective \ell dependences of D_{\ell}.

    Parameters
    ----------
    band: str
        95GHz or 150GHz or 220GHz for which we need the effective SPT3G band centre
    component : str
        The foreground component to use. Must be one of
        'tsz', 'dg-cl', 'dg-po', 'ksz', 'rg'
    fg_model: str
        Must be george (SPT-SZ G15: https://arxiv.org/abs/1408.3161) 
        or 
        reichardt (SPT-SZ/SPTpol R20: https://arxiv.org/abs/2002.06197 )
    dg_clus_template_id: int
        template for dg-cl
        0: contribution is split into 1- and 2-halo terms. Default choice.
        1: D_{\ell} \propto \ell^0.8.
    lmax: int
        lmax upto which template must be calculated. Default is 12000.
    Returns
    -------
    SPT-SZ/SPTpol pivot \ell, foreground amplitudes [uK^2] at pivot \ell, 
    foreground template filenames. 
    """
    component = component.lower()
    if fg_model not in ['george', 'reichardt']:
        raise NotImplementedError(fg_model)
    assert component in ['tsz', 'dg-cl', 'dg-po', 'dg', 'rg', 'ksz']

    if band == '90GHz':
        band = '95GHz'

    if band not in ['95GHz', '150GHz', '220GHz']:
        raise ValueError("band must be one of 90/95, 150, 220GHz")

    clustered_component_template_dic = {}

    if fg_model == 'george':
        el_norm = 3000
        tsz_freq0 = 143

        # D_{\ell}[\ell = el_norm] in uk^2
        tsz_power = 4.38  # sec 6.1.3 / table 6 of https://arxiv.org/pdf/1408.3161.pdf
        fg_amplitudes = {
            # sec 6.1.1 of https://arxiv.org/pdf/1408.3161.pdf
            'dg-po': {'95GHz': 1.37, '150GHz': 9.16, '220GHz': 66.1},
            # sec 6.1.1 of https://arxiv.org/pdf/1408.3161.pdf
            # dg-cl: dictionary containing amplitude for
            # 0: 1- and 2-halo terms
            # 1: D_{\ell} \propto \ell^0.8
            'dg-cl': {
                '95GHz': {0: [0.128, 0.122], 1: 0.208},
                '150GHz': {0: [1.99, 1.90], 1: 3.46},
                '220GHz': {0: [28.3, 26.5], 1: 50.6},
            },
            # sec 6.1.2 of https://arxiv.org/pdf/1408.3161.pdf
            'rg': {'95GHz': 7.81, '150GHz': 1.06, '220GHz': None},
            # sec 8.1 of https://arxiv.org/pdf/1408.3161.pdf
            # bispectrum constraints as sec 6.1.3 only quotes an upper limit
            'ksz': {'95GHz': 2.9, '150GHz': 2.9, '220GHz': 2.9},
        }

        # released as part of likelihood code tarball
        # with R20/G15 data release.
        # downloaded from https://pole.uchicago.edu/public/data/reichardt20/
        # Srini restricted these files to \ell_max = 12000 on 20200513
        # for SPT-3G simulation purposes.
        # Look into the above link to recover values upto \ell < 25000.
        fg_template_files = {
            'dg-cl': ['dl_cib_1halo_norm1_12000.txt', 'dl_cib_2halo_norm1_12000.txt'],
            'ksz': 'dl_ksz_CSFplusPATCHY_13sep2011_norm1_fake12000.txt',
            'tsz': 'dl_shaw_tsz_s10_153ghz_norm1_fake12000.txt',
        }

    elif fg_model == 'reichardt':
        el_norm = 3000
        tsz_freq0 = 143

        # D_{\ell}[\ell = el_norm] in uk^2
        tsz_power = 3.42  # sec 6.1.3 of https://arxiv.org/pdf/2002.06197.pdf
        fg_amplitudes = {
            # sec 6.1.1 of https://arxiv.org/pdf/2002.06197.pdf
            'dg-po': {'95GHz': None, '150GHz': 7.24, '220GHz': 61.4},
            # sec 6.1.1 of https://arxiv.org/pdf/2002.06197.pdf
            # dg-cl: dictionary containing amplitude for
            # 0: 1- and 2-halo terms
            # 1: D_{\ell} \propto \ell^0.8
            'dg-cl': {
                '95GHz': None,
                '150GHz': {0: [2.21, 1.82], 1: None},  # 4.03},
                '220GHz': {0: [32.4, 27.5], 1: None},  # 59.9},
            },
            # sec 6.1.2 of https://arxiv.org/pdf/2002.06197.pdf
            'rg': {'95GHz': None, '150GHz': 1.01, '220GHz': None},
            # sec 6.1.3 of https://arxiv.org/pdf/2002.06197.pdf
            'ksz': {'95GHz': 3.0, '150GHz': 3.0, '220GHz': 3.0},
        }

        # released as part of likelihood code tarball
        # with R20/G15 data release.
        # downloaded from https://pole.uchicago.edu/public/data/reichardt20/
        # Srini restricted these files to \ell_max = 12000 on 20200513
        # for SPT-3G simulation purposes.
        # Look into the above link to recover values upto \ell < 25000.
        fg_template_files = {
            'dg-cl': ['dl_cib_1halo_norm1_12000.txt', 'dl_cib_2halo_norm1_12000.txt'],
            'ksz': 'dl_ksz_CSFplusPATCHY_13sep2011_norm1_fake12000.txt',
            'tsz': 'dl_shaw_tsz_s10_153ghz_norm1_fake12000.txt',
        }

    # scaling the tsz power to 95/150/220 from tsz_freq0 = 143GHz used in G15/R20.
    # uses the effective frequencies of each band
    tsz_amp_dic = {}
    y_power = (
        tsz_power / (compton_y_to_delta_cmb_temp(tsz_freq0)) ** 2.0
    )
    for tmpband in ['95GHz', '150GHz', '220GHz']:
        tmpband_eff = get_spt_effective_frequencies(fg_model, tmpband, 'tsz')
        tsz_amp_dic[tmpband] = (
            y_power * (compton_y_to_delta_cmb_temp(tmpband_eff)) ** 2.0,
        )
    fg_amplitudes['tsz'] = tsz_amp_dic

    ####
    # get the foreground templates
    el = np.arange(lmax+1)

    if component == 'dg-po' or component == 'rg':
        # Poisson component - Dl \propto \ell^2.
        reqd_template = (el * 1.0 / el_norm) ** 2.0

    elif component == 'dg-cl':

        if dg_clus_template_id == 0:  # 1- and 2-halo terms
            template_fname_arr = fg_template_files['dg-cl']
            reqd_template = []
            for template_fname in template_fname_arr:

                template_fname = '%s/%s' %(data_folder, template_fname)

                template_el, template_dl = np.loadtxt(
                    template_fname, usecols=[0, 1], unpack=1
                )[:lmax]
                # get it at required el
                template = np.interp(el, template_el, template_dl, left=0.0, right=0.0)
                reqd_template.append(template)

            reqd_template = np.asarray(reqd_template)

        elif dg_clus_template_id == 1:  # Clustered component - Dl \propto \ell^0.8.
            reqd_template = (el * 1.0 / el_norm) ** 0.8

    else:  # tsz/ksz

        template_fname = fg_template_files[component]

        template_fname = '%s/%s' %(data_folder, template_fname)

        template_el, template_dl = np.loadtxt(template_fname, usecols=[0, 1], unpack=1)[:lmax]

        # get it at required el
        reqd_template = np.interp(el, template_el, template_dl, left=0.0, right=0.0)

    ####

    # get foreground amplitudes
    if component == 'dg-cl':
        amp_val = np.asarray(fg_amplitudes[component][band][dg_clus_template_id])
    else:
        amp_val = fg_amplitudes[component][band]

    if amp_val is None:
        raise ValueError(
            'Amplitude of %s is not defined in %sGHz band. Try other bands.'
            % (component, band)
        )

    return el, amp_val, reqd_template


def get_cl_dusty_galaxies(
    freq1,
    freq2=None,
    fg_model='reichardt',
    dg_clus_template_id=0,
    band0='150GHz',
    spec_index_dg_po=1.48,
    spec_index_dg_clus=2.23,
    cib_temp=25.0,
    return_ells=False,
):

    """
    Auto/Cross power spectra of dust galaxies
    (both Possion and clustering components) by scaling the SPT-SZ/SPTpol results.
    (i.e:) spectra will be scaled from SPT (freq0 x freq0) to 
    desired (freq1 x freq2) frequencies using the spectral indices.

    Parameters
    ----------
    freq1 : float
        Observing frequency in g3units.
    freq2 : float
        Observing frequency in g3units. If None, we will compute the autospectra (freq1 x freq1).
    fg_model: str
        Must be george (SPT-SZ G15: https://arxiv.org/abs/1408.3161) 
        or 
        reichardt (SPT-SZ/SPTpol R20: https://arxiv.org/abs/2002.06197 )
    dg_clus_template_id: int
        template for dg-cl
        0: contribution is split into 1- and 2-halo terms. Default choice.
        1: D_{\ell} \propto \ell^0.8.
    band0: str
        reference frequency string. Must be one the three SPT-SZ bands [90/95, 150, 220]GHz.
        This is the band at which the G15 or R20 models were calculated. 
        150GHz is the default value. 
        we will get the spt effective frequency for this band as freq0 and 
        the spectra will be scaled from SPT (freq0 x freq0) 
        to the desired (freq1 x freq2) frequencies using the specified spectral indices.
    spec_index_dg_po: float
        Spectral index for dusty galaxies (DG) Poisson component 
        (c.f. Eq. 13 of G15 https://arxiv.org/pdf/1408.3161.pdf)
        Default value is 1.48 from R20 (sec 6.1.1 https://arxiv.org/pdf/2002.06197.pdf).
        For G15 value check Sec 6.1.1 of G15 https://arxiv.org/pdf/1408.3161.pdf.
    spec_index_dg_clus: float
        Spectral index for dusty galaxies (DG) clustering component 
        (c.f. Eq. 13 of G15 https://arxiv.org/pdf/1408.3161.pdf).
        Default value is 2.23 from R20 (sec 6.1.1 https://arxiv.org/pdf/2002.06197.pdf).
        For G15 value check Sec 6.1.1 of G15 https://arxiv.org/pdf/1408.3161.pdf.
    cib_temp: float
        Dust Blackbody temperature in g3units.
        Default value is R20 value cib_temp = 25 K
        Cannot find the reference for this in https://arxiv.org/pdf/2002.06197.pdf.
        For G15 it is 20 K
        (Sec 5.4.3 of G15 https://arxiv.org/pdf/1408.3161.pdf).
    return_ells: bool
        if True, return ells else only the power spectra

    Returns
    -------
    if return_ells: el, cl_dg_po, cl_dg_clus: Multipoles, freq1 x freq2 power spectra of DG-Po and DG-Cl.
    else: cl_dg_po, cl_dg_clus: freq1 x freq2 power spectra of DG-Po and DG-Cl.
    """

    # check if freq2 is specified for cross spectra (freq1 x freq2),
    # else get auto spectra (freq1 x freq1)
    if freq2 is None:
        freq2 = freq1

    # ensure band0 is one of the SPT bands
    if band0 not in ['90GHz', '150GHz', '220GHz']:
        raise ValueError('band0 must be one of 90, 150, 220 GHz')

    # first get the SPT spectra at the reference frequency band0
    if fg_model not in ['george', 'reichardt']:
        sys.exit()

    (
        el,
        dl_dg_po_amp,
        dl_dg_po_template,
    ) = get_spt_foreground_amplitudes_templates(
        'DG-Po', band0, fg_model=fg_model,
    )

    (
        el,
        dl_dg_clus_amp,
        dl_dg_clus_template,
    ) = get_spt_foreground_amplitudes_templates(
        'DG-Cl',
        band0,
        fg_model=fg_model,
        dg_clus_template_id=dg_clus_template_id,
    )

    dl_fac = el * (el + 1) / 2 / np.pi


    # ref. frequency for reference band in the ref. SPT experiment
    freq0 = get_spt_effective_frequencies(fg_model, band0, 'DG')

    # Eq. 11 of G15 https://arxiv.org/pdf/1408.3161.pdf
    nr = (get_planck_db_dt(freq0)) ** 2.0
    dr = get_planck_db_dt(freq1) * get_planck_db_dt(freq2)
    epsilon_nu1_nu2 = nr / dr

    # Eq. 13 of G15 https://arxiv.org/pdf/1408.3161.pdf
    bnuT_freq1 = get_planck_bnu_of_t(freq1, temp=cib_temp)
    bnuT_freq2 = get_planck_bnu_of_t(freq2, temp=cib_temp)
    bnuT_freq0 = get_planck_bnu_of_t(freq0, temp=cib_temp)

    # ensure freq is in Hz for subsequent calculations
    freq0 *= 1e9
    freq1 *= 1e9
    freq2 *= 1e9

    # Eq. 13 of G15 https://arxiv.org/pdf/1408.3161.pdf for DG-Po
    etanu1_dg_po = ((1.0 * freq1) ** spec_index_dg_po) * bnuT_freq1
    etanu2_dg_po = ((1.0 * freq2) ** spec_index_dg_po) * bnuT_freq2
    etanu0_dg_po = ((1.0 * freq0) ** spec_index_dg_po) * bnuT_freq0

    # Eq. 13 of G15 https://arxiv.org/pdf/1408.3161.pdf for DG-Cl
    etanu1_dg_clus = ((1.0 * freq1) ** spec_index_dg_clus) * bnuT_freq1
    etanu2_dg_clus = ((1.0 * freq2) ** spec_index_dg_clus) * bnuT_freq2
    etanu0_dg_clus = ((1.0 * freq0) ** spec_index_dg_clus) * bnuT_freq0

    # convert spectra at (freq0 x freq0) to (freq1 x freq2)
    # Eq. 12 of G15 https://arxiv.org/pdf/1408.3161.pdf
    # DG-Po
    eps_eta_term_dg_po = epsilon_nu1_nu2 * (
        1.0 * etanu1_dg_po * etanu2_dg_po / etanu0_dg_po / etanu0_dg_po
    )
    print(freq0/1e9, freq1/1e9, freq2/1e9, eps_eta_term_dg_po)
    dl_dg_po = dl_dg_po_amp * eps_eta_term_dg_po * dl_dg_po_template

    # DG-Cl
    eps_eta_term_dg_clus = epsilon_nu1_nu2 * (
        1.0 * etanu1_dg_clus * etanu2_dg_clus / etanu0_dg_clus / etanu0_dg_clus
    )
    if np.ndim(dl_dg_clus_template) > 1:  # 1-/2-halo terms
        dl_dg_clus = (
            dl_dg_clus_amp[:, None] * eps_eta_term_dg_clus * dl_dg_clus_template
        )
        dl_dg_clus = np.sum(dl_dg_clus, axis=0)
    else:  # D_{\ell} \propto \ell^0.8
        dl_dg_clus = dl_dg_clus_amp * eps_eta_term_dg_clus * dl_dg_clus_template

    # convert dl to cl
    cl_dg_po = dl_dg_po / dl_fac
    cl_dg_clus = dl_dg_clus / dl_fac

    # remove nan and inf
    cl_dg_po[np.isnan(cl_dg_po)] = 0.0
    cl_dg_po[np.isinf(cl_dg_po)] = 0.0
    cl_dg_clus[np.isnan(cl_dg_clus)] = 0.0
    cl_dg_clus[np.isinf(cl_dg_clus)] = 0.0

    if not return_ells:
        return cl_dg_po, cl_dg_clus
    else:
        return el, cl_dg_po, cl_dg_clus


def get_cl_tsz(
    freq1,
    freq2=None,
    fg_model='reichardt',
    band0='150GHz',
    return_ells=False,
):

    """
    Auto/Cross tSZ power spectra by scaling the SPT-SZ/SPTpol results.
    (i.e:) spectra will be scaled from SPT (freq0 x freq0) to 
    the desired (freq1 x freq2) frequencies using the spectral indices.

    Parameters
    ----------
    freq1 : float
        Observing frequency in g3units.
    freq2 : float
        Observing frequency in g3units. If None, we will compute the autospectra (freq1 x freq1).
    fg_model: str
        Must be george (SPT-SZ G15: https://arxiv.org/abs/1408.3161) 
        or reichardt (SPT-SZ/SPTpol R20: https://arxiv.org/abs/2002.06197 )
    band0: str
        reference frequency string. Must be one the three SPT-SZ bands [90/95, 150, 220]GHz.
        This is the band at which the G15 or R20 models were calculated. 
        150GHz is the default value. 
        we will get the spt effective frequency for this band as freq0 and 
        the spectra will be scaled from SPT (freq0 x freq0) 
        to the desired (freq1 x freq2) frequencies using the specified spectral indices.
    return_ells: bool
        if True, return ells else only the power spectra

    Returns
    -------
    if return_ells: el, cl_tsz - Multipoles, freq1 x freq2 tSZ power spectra.
    else: cl_tsz - freq1 x freq2 tSZ power spectra.
    """

    # check if freq2 is specified for cross spectra (freq1 x freq2),
    # else get auto spectra (freq1 x freq1)
    if freq2 is None:
        freq2 = freq1

    # ensure band0 is one of the SPT bands
    if band0 not in ['90GHz', '150GHz', '220GHz']:
        raise ValueError('band0 must be one of 90, 150, 220 GHz')

    # first get the SPT spectra at the reference frequency band0
    if fg_model not in ['george', 'reichardt']:
        raise NotImplementedError(model)

    el, dl_tsz_amp, dl_tsz_template = get_spt_foreground_amplitudes_templates(
        'tSZ', band0, fg_model=fg_model
    )
    dl_fac = el * (el + 1) / 2 / np.pi
    dl_tsz_freq0 = dl_tsz_amp * dl_tsz_template
    cl_tsz_freq0 = dl_tsz_freq0 / dl_fac

    # ref. frequency for reference band in the ref. SPT experiment
    freq0 = get_spt_effective_frequencies(fg_model, band0, 'tSZ')

    # get the Compton-y to \DeltaTcmb conversion factors for all frequencies
    tsz_fac_freq0 = compton_y_to_delta_cmb_temp(freq0)
    tsz_fac_freq1 = compton_y_to_delta_cmb_temp(freq1)
    tsz_fac_freq2 = compton_y_to_delta_cmb_temp(freq2)

    # get the scaling factor to convert spectra at (freq0 x freq0) to (freq1 x freq2)
    scalefac = tsz_fac_freq1 * tsz_fac_freq2 / (tsz_fac_freq0 ** 2.0)

    cl_tsz = cl_tsz_freq0 * scalefac

    cl_tsz[np.isnan(cl_tsz)] = 0.0
    cl_tsz[np.isinf(cl_tsz)] = 0.0

    if not return_ells:
        return cl_tsz
    else:
        return el, cl_tsz


def get_cl_ksz(
    fg_model='reichardt', return_ells=False,
):

    """
    get kSZ power spectra from R20/G15.

    Parameters
    ----------
    fg_model: str
        Must be george (SPT-SZ G15: https://arxiv.org/abs/1408.3161) 
        or reichardt (SPT-SZ/SPTpol R20: https://arxiv.org/abs/2002.06197 )
    return_ells: bool
        if True, return ells else only the power spectra

    Returns
    -------
    if return_ells: el, cl_ksz - kSZ power spectra.
    else: cl_ksz - kSZ power spectra.
    """

    if fg_model not in ['george', 'reichardt']:
        raise NotImplementedError(model)

    el, dl_ksz_amp, dl_ksz_template = get_spt_foreground_amplitudes_templates(
        'kSZ', fg_model=fg_model
    )
    dl_fac = el * (el + 1) / 2 / np.pi
    dl_ksz = dl_ksz_amp * dl_ksz_template
    cl_ksz = dl_ksz / dl_fac

    cl_ksz[np.isnan(cl_ksz)] = 0.0
    cl_ksz[np.isinf(cl_ksz)] = 0.0


    if not return_ells:
        return cl_ksz
    else:
        return el, cl_ksz


def get_cl_radio(
    freq1,
    freq2=None,
    band0='150GHz',
    fg_model='reichardt',
    spec_index_rg=-0.76,
    return_ells=False,
):

    """
    Auto/Cross power spectra of radio galaxies (Poisson component) by scaling the SPT-SZ/SPTpol results.
    (i.e:) spectra will be scaled from SPT (freq0 x freq0) 
    to the desired (freq1 x freq2) frequencies using the spectral indices.

    Parameters
    ----------
    freq1 : float
        Observing frequency in g3units.
    freq2 : float
        Observing frequency in g3units. If None, we will compute the autospectra (freq1 x freq1).
    fg_model: str
        Must be george (SPT-SZ G15: https://arxiv.org/abs/1408.3161) 
        or reichardt (SPT-SZ/SPTpol R20: https://arxiv.org/abs/2002.06197 )
    band0: str
        reference frequency string. Must be one the three SPT-SZ bands [90/95, 150, 220]GHz.
        This is the band at which the G15 or R20 models were calculated. 
        150GHz is the default value. 
        we will get the spt effective frequency for this band as freq0 and 
        the spectra will be scaled from SPT (freq0 x freq0) 
        to the desired (freq1 x freq2) frequencies using the specified spectral indices.
    spec_index_rg: float
        Spectral index for radio galaxies (RG) Poisson component 
        (c.f. Eq. 16 of G15 https://arxiv.org/pdf/1408.3161.pdf)
        Default value is -0.76 from R20 (sec 6.1.2 https://arxiv.org/pdf/2002.06197.pdf).
        For G15 value check Sec 6.1.2 of G15 https://arxiv.org/pdf/1408.3161.pdf.
    return_ells: bool
        if True, return ells else only the power spectra

    Returns
    -------
    if return_ells: el, cl_rg - Multipoles, freq1 x freq2 power spectra of radio galaxies.
    else: cl_rg - freq1 x freq2 power spectra of radio galaxies.
    """

    # check if freq2 is specified for cross spectra (freq1 x freq2),
    # else get auto spectra (freq1 x freq1)
    if freq2 is None:
        freq2 = freq1

    # ensure band0 is one of the SPT bands
    if band0 not in ['90GHz', '150GHz', '220GHz']:
        raise ValueError('band0 must be one of 90, 150, 220 GHz')

    # first get the SPT spectra at the reference frequency band0
    if fg_model not in ['george', 'reichardt']:
        raise NotImplementedError(model)

    el, dl_rg_amp, dl_rg_template = get_spt_foreground_amplitudes_templates(
        'rg', band0, fg_model=fg_model
    )
    dl_fac = el * (el + 1) / 2 / np.pi

    # ref. frequency for reference band in the ref. SPT experiment
    freq0 = get_spt_effective_frequencies(fg_model, band0, 'RG')

    # Eq. 11 of G15 https://arxiv.org/pdf/1408.3161.pdf
    nr = (get_planck_db_dt(freq0)) ** 2.0
    dr = get_planck_db_dt(freq1) * get_planck_db_dt(freq2)
    epsilon_nu1_nu2 = nr / dr

    # convert spectra at (freq0 x freq0) to (freq1 x freq2)
    # Eq. 16 of G15 https://arxiv.org/pdf/1408.3161.pdf
    dl_rg = (
        dl_rg_amp
        * epsilon_nu1_nu2
        * (1.0 * freq1 * freq2 / freq0 / freq0) ** spec_index_rg
        * dl_rg_template
    )

    cl_rg = dl_rg / dl_fac

    cl_rg[np.isnan(cl_rg)] = 0.0

    if not return_ells:
        return cl_rg
    else:
        return el, cl_rg


def get_planck_bnu_of_t(freq, temp=2.725):

    """
    get Planck function B_nu(T)

    Parameters
    ----------
    freq : float
        Observing frequency in g3units.
    temp: float
        Blackbody tempature in g3units.

    Returns
    -------
    Planck BnuT : float
        Sepctral radiance of given frequency at temperature 
    """

    freq *= 1e9
    
    x = h * freq / (k_B * temp)

    t1 = 2 * h * freq ** 3.0 / c ** 2.0
    t2 = 1.0 / (np.exp(x) - 1.0)

    return t1 * t2


def get_planck_db_dt(freq, freq0=None, temp=2.725):

    """
    get Planck dB/dT

    Parameters
    ----------
    freq : float
        Observing frequency in g3units.
    freq0 : float, optional
        reference frequency in g3units.
    temp: float, optional
        Blackbody tempature in g3units.

    Returns
    -------
    Planck dB/dT : float
    """

    freq *= 1e9

    x = h * freq / (k_B * temp)
    dBdT = x ** 4.0 * np.exp(x) / (np.exp(x) - 1) ** 2.0

    if freq0 is not None:
        freq0 *= 1e9
        x0 = h * freq0 / (k_B * temp)
        dBdT0 = x0 ** 4 * np.exp(x0) / (np.exp(x0) - 1) ** 2.0
        return dBdT / dbdT0
    else:
        return dBdT


def coth(x):
    return (np.exp(x) + np.exp(-x)) / (np.exp(x) - np.exp(-x))


def compton_y_to_delta_cmb_temp(freq):

    """
    Compton-$y$ to $\Delta$cmb_temp

    Parameters
    ----------
    freq : float
        Observing frequency in g3units.

    Returns
    -------
    conv_factor : float
        Conversion factor to go from Compton-y to 
        \Delta cmb_temp at the specified frequency
    """

    freq *= 1e9
    x = (h * freq) / (k_B * cmb_temp)
    g_nu = x * coth(x / 2.0) - 4.0

    return cmb_temp * g_nu
