import numpy as np, sys, os, scipy as sc, foregrounds as fg, misc, re, flatsky#, healpy as H
from pylab import *
################################################################################################################
def get_analytic_covariance(param_dict, freqarr, el = None, nl_dic = None, bl_dic = None, ignore_fg = [], which_spec = 'TT', pol_frac_per_cent_dust = 0.02, pol_frac_per_cent_radio = 0.03, pol_frac_per_cent_tsz = 0., pol_frac_per_cent_ksz = 0., include_gal = 0, max_nl_value = 5000., beam_tol_for_ilc = 1000., cib_corr_coeffs = None, use_websky_cib = 0, scale_spt_using_sptspire = 0, use_sptspire_for_hfbands = 0, minval_for_hfbands=500, use_mdpl2_cib = 0, null_highfreq_radio = 1, reduce_radio_power_150 = None, reduce_tsz_power = None, reduce_cib_power = None, remove_cib_decorr = 0, cib_flux_threshold = 1.5, mdpl2_cib_version = 'v0p3', use_mdpl2_tsz = 0, cl_multiplier_dic = None, return_fg_spectra = True, force_cl_dic = None):

    #ignore_fg = foreground terms that must be ignored
    debug=False
    possible_ignore_fg = ['cmb', 'tsz', 'y', 'ksz', 'radio', 'dust', 'noise', 'tsz_cib']
    if len(ignore_fg)>0:
        if 'cmb' in ignore_fg: ignore_fg.append('ksz')
        if not all( [ currfg in possible_ignore_fg for currfg in ignore_fg] ):
            print( '\n\t Alert: Elements of ignore_fg should be one of the following: %s\n\n' %(np.array2string(np.asarray(possible_ignore_fg))) )
            sys.exit()


    el_, cl_cmb = fg.get_foreground_power_spt('CMB', freq1 = param_dict['freq0'], freq2 = param_dict['freq0'])
    if el is None: el = np.copy(el_)    
    el_, cl_ksz = fg.get_foreground_power_spt('kSZ', freq1 = param_dict['freq0'], freq2 = param_dict['freq0'])

    if (1): #make cmb, ksz length the same as el
        cl_cmb = np.interp(el, np.arange(len(cl_cmb)), cl_cmb)
        cl_ksz = np.interp(el, np.arange(len(cl_ksz)), cl_ksz)

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
    #20220428 - force cl if force_cl_dic is supplied

    cl_dic = {}
    cl_ori = np.zeros(len(el))
    if use_sptspire_for_hfbands or scale_spt_using_sptspire:
        param_dict['reduce_radio_power_150'] = reduce_radio_power_150
        comps_to_subtract_from_spt_spire = ['CMB', 'kSZ', 'tSZ', 'radio']
        spt_spire_freq_crosses_dic = fg.get_spt_spire_bandpower(el_for_interp = el, comps_to_subtract = comps_to_subtract_from_spt_spire, param_dict = param_dict)
    else:
        spt_spire_freq_crosses_dic = None

    if return_fg_spectra:
        fg_cl_dic = {}

    for freq1 in freqarr:
        for freq2 in freqarr:

            #20220429: removing this as tsz x cib is not symmetric.
            '''
            if (freq2, freq1) in cl_dic:
                cl_dic[(freq1, freq2)] = cl_dic[(freq2, freq1)]
                continue
            '''

            #get tsz
            if use_mdpl2_tsz and use_mdpl2_cib: #20201120  - use MDPL2 tsz spectra for SPT3G/ Planck bands
                el_, cl_tsz = fg.get_cl_tsz_tszcib_mdpl2_v0p3(freq1, freq2, el = el, which_spec = 'tsz', reduce_tsz_power = reduce_tsz_power)
            else:
                el_, cl_tsz = fg.get_cl_tsz(freq1, freq2, freq0 = param_dict['freq0'], fg_model = param_dict['fg_model'], reduce_tsz_power = reduce_tsz_power)

            if (0):
                el_, cl_tsz_G15 = fg.get_cl_tsz(freq1, freq2, freq0 = param_dict['freq0'], fg_model = param_dict['fg_model'], reduce_tsz_power = reduce_tsz_power)
                dl_fac = (el * (el+1))/2./np.pi
                ax=subplot(111,yscale='log'); plot(el, dl_fac * cl_tsz_G15); plot(el, dl_fac * cl_tsz ); title('tSZ: %s,%s: %s' %(freq1, freq2, which_spec))
                ylim(0.1, 1e4); xlim(2000, 1e4)
                if freq1 == freq2:
                    freq0 = 150
                    el_, cl_tsz_freq0 = fg.get_cl_tsz_tszcib_mdpl2_v0p3(freq0, freq0, el = el, which_spec = 'tsz', reduce_tsz_power = reduce_tsz_power)
                    tsz_fac_freq0 = compton_y_to_delta_Tcmb(freq0*1e9)
                    tsz_fac_freq1 = compton_y_to_delta_Tcmb(freq1*1e9)
                    scalefac = tsz_fac_freq1 * tsz_fac_freq1/ (tsz_fac_freq0**2.)
                    cl_tsz_v2 = cl_tsz_freq0 * scalefac
                    plot(el, dl_fac * cl_tsz_v2 );
                show(); #sys.exit()

            if which_spec == 'EE':
                cl_tsz = cl_tsz * pol_frac_per_cent_tsz**2.
            elif which_spec == 'TE':
                cl_tsz = cl_tsz * 0.

            #get radio
            el_, cl_radio = fg.get_cl_radio(freq1, freq2, freq0 = param_dict['freq0'], fg_model = param_dict['fg_model'], spec_index_rg = param_dict['spec_index_rg'], null_highfreq_radio = null_highfreq_radio, reduce_radio_power_150 = reduce_radio_power_150)
            if which_spec == 'EE':
                cl_radio = cl_radio * pol_frac_per_cent_radio**2.
            elif which_spec == 'TE':
                cl_radio = cl_radio * 0.

            #get dust
            #sys.exit()
            el_,  cl_dg_po, cl_dg_clus = fg.get_cl_dust(freq1, freq2, freq0 = param_dict['freq0'], fg_model = param_dict['fg_model'], spec_index_dg_po = param_dict['spec_index_dg_po'], spec_index_dg_clus = param_dict['spec_index_dg_clus'], Tcib = param_dict['Tcib'], reduce_cib_power = reduce_cib_power)
            cl_dust = cl_dg_po + cl_dg_clus
            if remove_cib_decorr:
                el_,  cl_dg_po1, cl_dg_clus1 = fg.get_cl_dust(freq1, freq1, freq0 = param_dict['freq0'], fg_model = param_dict['fg_model'], spec_index_dg_po = param_dict['spec_index_dg_po'], spec_index_dg_clus = param_dict['spec_index_dg_clus'], Tcib = param_dict['Tcib'], reduce_cib_power = reduce_cib_power)
                el_,  cl_dg_po2, cl_dg_clus2 = fg.get_cl_dust(freq2, freq2, freq0 = param_dict['freq0'], fg_model = param_dict['fg_model'], spec_index_dg_po = param_dict['spec_index_dg_po'], spec_index_dg_clus = param_dict['spec_index_dg_clus'], Tcib = param_dict['Tcib'], reduce_cib_power = reduce_cib_power)
                cl_dust1 = cl_dg_po1 + cl_dg_clus1
                cl_dust2 = cl_dg_po2 + cl_dg_clus2
                cl_dust = np.sqrt( cl_dust1 * cl_dust2 )
                cib_corr_coeffs = None

            cl_dust[np.isnan(cl_dust)] = 0.
            cl_dust = np.interp(el, el_, cl_dust)
            cl_dust_ori = np.copy(cl_dust)
            tit = 'G15/R20 CIB'
            if use_websky_cib and (force_cl_dic is None or len(force_cl_dic) == 0):
                el_,  cl_dust = fg.get_cl_cib_websky(freq1, freq2, el = el, remove_cib_decorr = remove_cib_decorr, threshold_mjy = cib_flux_threshold)
                cib_corr_coeffs = None #do not use this as websky already takes it into account
                tit = 'Websky CIB'
            if use_mdpl2_cib and (force_cl_dic is None or len(force_cl_dic) == 0):
                #el_,  cl_dust = fg.get_cl_cib_mdpl2(freq1, freq2, el = el)
                el_,  cl_dust = fg.get_cl_cib_mdpl2_v0p3(freq1, freq2, el = el, remove_cib_decorr = remove_cib_decorr, flux_threshold = cib_flux_threshold, v0p3_or_v0p5 = mdpl2_cib_version)
                cib_corr_coeffs = None #do not use this as websky already takes it into account
                tit = 'MDPL2 CIB'                
            if use_sptspire_for_hfbands or scale_spt_using_sptspire:
                if freq1>minval_for_hfbands or freq2>minval_for_hfbands:
                    #el, cl_dust = fg.get_spt_spire_bandpower(freq1, freq2, el_for_interp = el)
                    el_, cl_spt_spire = spt_spire_freq_crosses_dic[(freq1, freq2)]
                    if use_sptspire_for_hfbands:
                        cl_dust = np.copy( cl_spt_spire )
                    elif scale_spt_using_sptspire: #20220311 - scale G15/R21 based on SPTxSPIRE
                        el_for_scaling = 3000
                        spt_spire_val = np.interp( el_for_scaling, el_, cl_spt_spire )
                        curr_val = np.interp( el_for_scaling, el, cl_dust )
                        scaling_val = spt_spire_val/curr_val
                        cl_dust = cl_dust * scaling_val

                    #20220311 - actually what is written below is not true, I think.
                    #SPTxSPIRE:
                        #cross spectra obviously does not include noise.
                        #auto spectra: also should not include noise because there are obtained by crossing different halves.
                    """
                    if which_spec == 'TT' and ( (freq1, freq2) in nl_dic or (freq2, freq1) in nl_dic ) and freq1>minval_for_hfbands or freq2>minval_for_hfbands: #null nl_TT as SPTxSPIRE bandpowers already inlcudes them.
                        nl_dic[(freq1, freq2)] = nl_dic[(freq2, freq1)] = np.zeros( len(nl_dic[(freq1, freq2)]) )
                    """
                    cib_corr_coeffs = None #do not use this as websky already takes it into account
                tit = 'SPTxSPIRE CIB'

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
                cl_dust *= corr_coeff

            #get tSZ x CIB
            if use_mdpl2_tsz and use_mdpl2_cib:#20201120  - use MDPL2 tsz x cib spectra for SPT3G/ Planck bands
                el_, cl_tsz_cib = fg.get_cl_tsz_tszcib_mdpl2_v0p3(freq1, freq2, el = el, which_spec = 'tsz_cib')
                if (0):
                    el_, cl_tsz_cib_G15 = fg.get_cl_tsz_cib(freq1, freq2, freq0 = param_dict['freq0'], fg_model = param_dict['fg_model'], spec_index_dg_po = param_dict['spec_index_dg_po'], spec_index_dg_clus = param_dict['spec_index_dg_clus'], Tcib = param_dict['Tcib'], use_websky_cib = use_websky_cib, use_sptspire_for_hfbands = use_sptspire_for_hfbands, use_mdpl2_cib = use_mdpl2_cib, cl_cib_dic = spt_spire_freq_crosses_dic, reduce_tsz_power = reduce_tsz_power, cib_flux_threshold = cib_flux_threshold, mdpl2_cib_version = mdpl2_cib_version)
                    dl_fac = (el * (el+1))/2./np.pi
                    ax=subplot(111,yscale='log'); plot(el, dl_fac * cl_tsz_cib_G15, label = 'G15'); plot(el, dl_fac * cl_tsz_cib, label = 'MDPL2'); title('%s,%s: %s' %(freq1, freq2, which_spec))
                    legend(loc = 1)
                    ylim(0.1, 1e4); xlim(2000, 1e4)
                    show()
            else: 
                el_, cl_tsz_cib = fg.get_cl_tsz_cib(freq1, freq2, freq0 = param_dict['freq0'], fg_model = param_dict['fg_model'], spec_index_dg_po = param_dict['spec_index_dg_po'], spec_index_dg_clus = param_dict['spec_index_dg_clus'], Tcib = param_dict['Tcib'], use_websky_cib = use_websky_cib, use_sptspire_for_hfbands = use_sptspire_for_hfbands, minval_for_hfbands = minval_for_hfbands,  use_mdpl2_cib = use_mdpl2_cib, cl_cib_dic = spt_spire_freq_crosses_dic, reduce_tsz_power = reduce_tsz_power, cib_flux_threshold = cib_flux_threshold, mdpl2_cib_version = mdpl2_cib_version)
            if which_spec == 'EE' or which_spec == 'TE':
                cl_tsz_cib = cl_tsz_cib * 0.

            if (0):
                freq_combs = [ (90, 90), (90,150), (90, 220), (150, 150), (150, 220), (220, 220)]
                colorarr = ['navy', 'blue', 'royalblue', 'green', 'lime', 'darkred']
                ax = subplot(111, yscale = 'log');
                for cntr, f1f2 in enumerate( freq_combs ):
                    freq1, freq2 = f1f2
                    el_, cl_tsz_cib = fg.get_cl_tsz_cib(freq1, freq2, freq0 = param_dict['freq0'], fg_model = param_dict['fg_model'], spec_index_dg_po = param_dict['spec_index_dg_po'], spec_index_dg_clus = param_dict['spec_index_dg_clus'], Tcib = param_dict['Tcib'])
                    el_, cl_tsz_cib_v2 = fg.get_foreground_power_spt('tSZ-CIB', freq1 = freq1, freq2 = freq2) 
                    plot(el, cl_tsz_cib, color = colorarr[cntr], lw = 1.); plot(el, -cl_tsz_cib_v2, color = colorarr[cntr], lw = 2., ls = '--'); 
                xlim(100, 5000)
                show()
                sys.exit()

            #galaxy
            if include_gal:# and not pol: #get galactic dust and sync
                el_, cl_gal_dust = fg.get_cl_galactic(param_dict, 'dust', freq1, freq2, which_spec, el = el, bl_dic = bl_dic)
                el_, cl_gal_sync = fg.get_cl_galactic(param_dict, 'sync', freq1, freq2, which_spec, el = el, bl_dic = bl_dic)

            cl = np.copy( cl_ori )

            #20220428 - force cl if force_cl_dic is supplied
            tsz_forced, cib_forced = False, False
            if 'y' in force_cl_dic:
                tsz_forced = True
                cl_y_force = force_cl_dic['y']
                tsz_fac_freq1 = compton_y_to_delta_Tcmb(freq1*1e9)
                tsz_fac_freq2 = compton_y_to_delta_Tcmb(freq2*1e9)
                tsz_fac = tsz_fac_freq1 * tsz_fac_freq2
                cl_tsz_force = cl_y_force * tsz_fac
                #loglog(cl_tsz, 'k'); loglog(cl_y_force, 'r'); loglog(cl_tsz_force, 'g'); show(); sys.exit()
                cl_tsz = np.interp(el, np.arange(len(cl_tsz_force)), cl_tsz_force)
            if 'tsz' in force_cl_dic:
                tsz_forced = True
                cl_tsz_force = force_cl_dic['tsz'][(freq1, freq2)]
                cl_tsz = np.interp(el, np.arange(len(cl_tsz_force)), cl_tsz_force)
            if 'cib' in force_cl_dic:
                cib_forced = True
                if (freq1, freq2) in force_cl_dic['cib']:
                    cl_dust_force = force_cl_dic['cib'][(freq1, freq2)]
                else:
                    cl_dust_force = force_cl_dic['cib'][(freq2, freq1)]
                #loglog(cl_dust_force); loglog(cl_dust); show(); sys.exit()
                cl_dust = np.interp(el, np.arange(len(cl_dust_force)), cl_dust_force)
            if 'tsz_cib' in force_cl_dic or 'cib_tsz' in force_cl_dic:
                tsz_cib_forced = True
                if 'tsz_cib' in force_cl_dic:
                    cl_tsz_cib_force = force_cl_dic['tsz_cib'][(freq1, freq2)]
                else:
                    cl_tsz_cib_force = force_cl_dic['cib_tsz'][(freq1, freq2)]
                #loglog(abs(cl_tsz_cib_force)); loglog(abs(cl_tsz_cib)); show(); sys.exit()
                cl_tsz_cib = np.interp(el, np.arange(len(cl_tsz_cib_force)), cl_tsz_cib_force)
            '''#20220507 - calculated tszxcib for mdpl2 and websky. so commenting these things out.       
            if 'y_cib' in force_cl_dic or 'cib_y' in force_cl_dic:
                #this is CIB at freq1 x y. So we need tsz_fac at freq2 to get cl_tszfreq2_cibfreq1
                cl_y_cib_force = force_cl_dic['y_cib'][freq1]
                #print(cl_y_cib_force, freq1); sys.exit()
                tsz_fac_freq2 = compton_y_to_delta_Tcmb(freq2*1e9)
                cl_tsz_cib_force = cl_y_cib_force * tsz_fac_freq2
                if (1):
                    cl_y_force = force_cl_dic['y']
                    tsz_fac_freq1 = compton_y_to_delta_Tcmb(freq1*1e9)
                    tsz_fac_freq2 = compton_y_to_delta_Tcmb(freq2*1e9)
                    cl_tsz_force_freq1_freq2 = cl_y_force * tsz_fac_freq1**2.
                    cl_tsz_force_freq2_freq2 = cl_y_force * tsz_fac_freq2**2.

                    cl_dust_force_freq1_freq1 = force_cl_dic['cib'][(freq1, freq1)]
                    cl_dust_force_freq2_freq2 = force_cl_dic['cib'][(freq2, freq2)]

                    cl_tsz_cib_force_v2 = fg.get_cl_tsz_cib_simple(freq1, freq2, cl_tsz_force_freq1_freq2, cl_tsz_force_freq2_freq2, cl_dust_force_freq1_freq1, cl_dust_force_freq2_freq2)

                if (1): #check mdpl2
                    import glob
                    mdpl2_fd = 'DRAFT/data/mdpl2/v0.3/tsz_and_cib/'
                    mdpl2_searchstr = '%s/cltszi*%s*%s*' %(mdpl2_fd, freq1, freq2)
                    mdpl2_fname = glob.glob(mdpl2_searchstr)[0]
                    mdpl2_tsz_cib_els, mdpl2_tsz_cib_cl = np.loadtxt(mdpl2_fname, unpack = True)
                    mdpl2_tsz_cib_dl_fac = mdpl2_tsz_cib_els * (mdpl2_tsz_cib_els+1)/2/np.pi

                if (1):

                    print(freq1, freq2)
                    clf()
                    cl_tsz_cib = cl_tsz_cib[:len(el)]
                    tmp_el, tmp_el_force = np.arange(len(cl_tsz_cib)), np.arange(len(cl_tsz_cib_force))
                    dl_fac, dl_fac_force = tmp_el * (tmp_el+1)/2/np.pi, tmp_el_force * (tmp_el_force+1)/2/np.pi
                    ax = subplot(111, yscale = 'log')
                    plot(tmp_el, dl_fac * cl_tsz_cib*2, color = 'black', ls = '-', label = r'G15 original'); plot(tmp_el, dl_fac * abs(cl_tsz_cib)*2, color = 'black', ls ='-.'); 
                    #plot(tmp_el, dl_fac * cl_dust, 'darkred', label = r'Websky CIB')
                    plot(tmp_el_force, dl_fac_force * cl_tsz_cib_force, color = 'red', ls = '-', label = r'Websky tSZ X CIB'); plot(tmp_el_force, dl_fac_force * abs(cl_tsz_cib_force), color = 'r', ls = '-.'); 
                    plot(tmp_el_force, dl_fac_force * cl_tsz_cib_force_v2*2, color = 'goldenrod', ls = '-', label = r'G15 parameterisation of websky tSZ X CIB'); plot(tmp_el_force, dl_fac_force * abs(cl_tsz_cib_force_v2)*2, color = 'goldenrod', ls = '-.'); 
                    plot(mdpl2_tsz_cib_els, mdpl2_tsz_cib_dl_fac * mdpl2_tsz_cib_cl, color = 'darkgreen', ls = '-', label = r'MDPL2 tSZ X CIB v0.3'); plot(mdpl2_tsz_cib_els, mdpl2_tsz_cib_dl_fac * abs(mdpl2_tsz_cib_cl), color = 'darkgreen', ls = '-.'); 
                    legend(loc = 1)
                    ylim(0.001, 1e3); 
                    xlim(100, 5000); title('%s x %s' %(freq1, freq2))
                    show(); #sys.exit()
                cl_tsz_cib = np.interp(el, np.arange(len(cl_tsz_cib_force)), cl_tsz_cib_force)

            #if force_cl_dic is not None and len(force_cl_dic)>0:#(cib_forced and tsz_forced) and 'y_cib' not in force_cl_dic and 'cib_y' not in force_cl_dic:
            if (cib_forced and tsz_forced) and 'y_cib' not in force_cl_dic and 'cib_y' not in force_cl_dic:
                #get tsz at freq1 and freq2
                cl_y_force = force_cl_dic['y']
                tsz_fac_freq1 = compton_y_to_delta_Tcmb(freq1*1e9)
                tsz_fac_freq2 = compton_y_to_delta_Tcmb(freq2*1e9)
                cl_tsz_force_freq1_freq2 = cl_y_force * tsz_fac_freq1**2.
                cl_tsz_force_freq2_freq2 = cl_y_force * tsz_fac_freq2**2.

                cl_dust_force_freq1_freq1 = force_cl_dic['cib'][(freq1, freq1)]
                cl_dust_force_freq2_freq2 = force_cl_dic['cib'][(freq2, freq2)]

                cl_tsz_cib_force = fg.get_cl_tsz_cib_simple(freq1, freq2, cl_tsz_force_freq1_freq2, cl_tsz_force_freq2_freq2, cl_dust_force_freq1_freq1, cl_dust_force_freq2_freq2)

                if (0):
                    clf()
                    ax = subplot(111, yscale = 'log')
                    plot(cl_tsz_cib, 'k-'); plot(abs(cl_tsz_cib), 'k-.'); 
                    plot(cl_tsz_cib_force, 'r-'); plot(abs(cl_tsz_cib_force), 'r-.'); 
                    show(); sys.exit()
                
                cl_tsz_cib = np.interp(el, np.arange(len(cl_tsz_cib_force)), cl_tsz_cib_force)
            '''
            if 'rad' in force_cl_dic or 'radio' in force_cl_dic:
                if 'rad' in force_cl_dic:
                    cl_radio_force = force_cl_dic['rad'][(freq1, freq2)]
                else:
                    cl_radio_force = force_cl_dic['radio'][(freq1, freq2)]
                cl_radio = np.interp(el, np.arange(len(cl_radio_force)), cl_radio_force)
            #print(cl_dust, cl_tsz, freq1, freq2); #sys.exit()
            #20220428 - force cl if force_cl_dic is supplied

            #if (1): #make a plot of CIB SPT x SPIRE interpolated + extended power spectra
            reqd_freq = 145 ##220 ##150 ##220 ##150 ##90
            if (0):#(freq1 == reqd_freq or freq2 == reqd_freq) and freq2>=freq1: 
                debug=True
                #if which_spec == 'TT' and (freq1==90): loglog(el, cl_dust, label = r'%s,%s' %(freq1,freq2)); 
                #if which_spec == 'TT' and (freq1==150): loglog(el, cl_dust, label = r'%s,%s' %(freq1,freq2)); 
                #if which_spec == 'TT' and (freq1==220): loglog(el, cl_dust, label = r'%s,%s' %(freq1,freq2)); 
                if (1):##which_spec == 'TT': 
                    freq_combs = []
                    color_dic = {}
                    shades_90 = [cm.Blues(int(d)) for d in np.linspace(100, 255, 5)]
                    shades_150 = [cm.Purples(int(d)) for d in np.linspace(100, 255, 4)]
                    shades_220 = [cm.Greens(int(d)) for d in np.linspace(100, 255, 3)]
                    shades_600 = ['red', 'maroon']
                    shades_857 = ['black']
                    shadearr = np.vstack( (shades_90, shades_150, shades_220 ))#, shades_600, shades_857) )
                    shadearr = shadearr.tolist()
                    shadearr.extend( shades_600 )
                    shadearr.extend( shades_857 )

                    for fcntr1, f1 in enumerate( freqarr ):
                        for fcntr2, f2 in enumerate( freqarr ):
                            if (f2, f1) in freq_combs: continue
                            freq_combs.append((f1,f2))
                    print(freq_combs, freq1, freq2)
                    freq_combs = np.asarray( freq_combs )
                    colorarr = [cm.jet(int(d)) for d in np.linspace(0, 255, len(freq_combs))]
                    colorarr = np.asarray( colorarr )
                    if freq1 == 90:
                        ls = ':'
                    elif freq1 == 150 or freq1 == 145:
                        ls = '-.'
                    elif freq1 == 220:
                        ls = '--'
                    else:
                        ls = '-'

                    ls = '-'

                    cind = np.where( (freq_combs[:,0] == freq1) & (freq_combs[:,1] == freq2) )[0][0]
                    #colorval = colorarr[cind]
                    if 345 in freqarr:
                        colorval = colorarr[cind]
                    else:
                        colorval = shadearr[cind]
                    #print(colorval)
                    try:
                        ax
                    except:
                        ax = subplot(111, yscale = 'log', xscale = 'log')
                    #print(freq1, freq2, cl_dust)
                    plot(el, cl_dust, color = colorval, ls = ls, label = r'%s,%s' %(freq1,freq2)); ylim(1e-7,1e6)
                    plot(el, cl_dust_ori, color = colorval, ls = ':')
                    #axvline(3000., ls = ':', lw = 0.25);axhline(1015., ls = ':', lw = 0.25)
                    #ax.yaxis.set_major_locator(MaxNLocator(nbins=10))
                    xlabel(r'Multipole $\ell$', fontsize = 14)
                    ylabel(r'C$_{\ell}$ $[\mu K^{2}]$', fontsize = 14)
                    title(r'%s spectra' %(tit))
                    ylim(1e-7, 1e5)
                    #show(); sys.exit()

            if cl_multiplier_dic is not None:
                if 'tsz' in cl_multiplier_dic:
                    cl_tsz = np.copy(cl_tsz) * cl_multiplier_dic['tsz']
                if 'radio' in cl_multiplier_dic:
                    cl_radio = np.copy(cl_radio) * cl_multiplier_dic['radio']
                if 'dust' in cl_multiplier_dic:
                    cl_dust = np.copy(cl_dust) * cl_multiplier_dic['dust']
                if 'tsz_cib' in cl_multiplier_dic:
                    cl_tsz_cib = np.copy(cl_tsz_cib) * cl_multiplier_dic['tsz_cib']
                if include_gal:
                    if 'gal_dust' in cl_multiplier_dic:
                        cl_gal_dust = np.copy(cl_gal_dust) * cl_multiplier_dic['gal_dust']
                    if 'gal_sync' in cl_multiplier_dic:
                        cl_gal_sync = np.copy(cl_gal_sync) * cl_multiplier_dic['gal_sync']

            if 'cmb' not in ignore_fg:
                if len(cl_cmb)<len(cl):
                    cl_cmb = np.interp(el, np.arange(len(cl_cmb)), cl_cmb)
                cl = cl + np.copy(cl_cmb[el])
            if 'ksz' not in ignore_fg:
                if len(cl_ksz)<len(cl):
                    cl_ksz = np.interp(el, np.arange(len(cl_ksz)), cl_ksz)
                #print('ksz')
                cl = cl + cl_ksz[el]
            if 'tsz' not in ignore_fg and 'y' not in ignore_fg:
                if len(cl_tsz)<len(cl):
                    cl_tsz = np.interp(el, np.arange(len(cl_tsz)), cl_tsz)
                #print('tsz')
                cl = cl + cl_tsz[el]
            if 'radio' not in ignore_fg:
                if len(cl_radio)<len(cl):
                    cl_radio = np.interp(el, np.arange(len(cl_radio)), cl_radio)
                #print('radio')
                cl = cl + cl_radio[el]
            if 'dust' not in ignore_fg:
                if len(cl_dust)<len(cl):
                    cl_dust = np.interp(el, np.arange(len(cl_dust)), cl_dust)
                #print('dust')
                cl = cl + cl_dust[el]

            #20220503 - add tszxcib if either tsz or cib is included.
            add_cl_tsz_cib = True
            #if ('dust' in ignore_fg and 'tsz' in ignore_fg) or 'tsz_cib' in ignore_fg or 'cib_tsz' in ignore_fg:
            if 'tsz_cib' in ignore_fg or 'cib_tsz' in ignore_fg or 'cib_y' in ignore_fg or 'y_cib' in ignore_fg:
                add_cl_tsz_cib = False 
            if add_cl_tsz_cib: #'dust' not in ignore_fg and 'tsz' not in ignore_fg and 'tsz_cib' not in ignore_fg:
                if len(cl_tsz_cib)<len(cl):
                    cl_tsz_cib = np.interp(el, np.arange(len(cl_tsz_cib)), cl_tsz_cib)
                cl = cl + cl_tsz_cib[el]

            if include_gal:# and not pol: #get galactic dust and sync

                if (0):#not pol:

                    loglog(cl_dg_po + cl_dg_clus, label = r'EG: %s, %s' %(freq1, freq2));
                    loglog(cl_gal_dust, label = r'Dust: %s, %s' %(freq1, freq2));title(r'%s' %(param_dict['which_gal_mask']))
                    loglog(cl_gal_sync, label = r'Sync: %s, %s' %(freq1, freq2));title(r'%s' %(param_dict['which_gal_mask']))
                    #cl = cl + cl_gal_dust + cl_gal_sync
                    #loglog(cl);
                    legend(loc = 1)
                    ylim(1e-5, 1e3);#show();sys.exit()
                    xlim(20,param_dict['lmax']);ylim(1e-8,1e6);
                    show();sys.exit()
                
                cl = cl + cl_gal_dust
                cl = cl + cl_gal_sync

            #make sure cl start from el=0 rather than el=10 which is the default for SPT G15 results
            lmin = min(el)
            cl = np.concatenate( (np.zeros(lmin), cl) )

            #change me: making cl to be simply dust
            #cl = cl_dust

            #noise auto power spectrum
            if nl_dic is not None:

                if (freq1,freq2) in nl_dic:
                    nl = nl_dic[(freq1,freq2)]
                else:
                    nl = nl_dic[freq1]
                    if freq1 != freq2: 
                        nl = np.copy(nl) * 0.

                if len(cl) > len(nl):
                    cl = cl[:len(nl)]
                elif len(cl) < len(nl):
                    nl = nl[:len(cl)]

                if (1): #20201121: remove very large numbers because of beam deconvolution                
                    ini_nl = np.median(nl[:100])
                    end_nl = np.median(nl[-100:])
                    if end_nl>ini_nl: #this implies beam deconvolution has made end nl pretty large
                        #having end_nl pretty large introduces covariance inversion issues                        
                        badinds = np.where(nl>=max_nl_value)[0]
                        nl[badinds] = max_nl_value
                        #print(ini_nl, end_nl)
                        #clf();loglog(nl); title('%s: %s,%s' %(which_spec, freq1,freq2)); show()

                el = np.arange(len(cl))

            else:
                nl = np.zeros(len(cl))
                #20191116 - fix this: there must be noise correlation in case of atmospheric noise
                print('\n\n\t\tfix me: there must be noise correlation in case of atmospheric noise')
                sys.exit()

            if cl_multiplier_dic is not None:
                if 'noise' in cl_multiplier_dic:
                    nl = np.copy(nl) * cl_multiplier_dic['noise']

            if 'noise' not in ignore_fg:
                if which_spec != 'TE':
                    cl = cl + np.copy(nl)

            if return_fg_spectra:
                if (1):#'cmb' not in ignore_fg:
                    if 'cmb' not in fg_cl_dic: fg_cl_dic['cmb'] = {}
                    fg_cl_dic['cmb'][(freq1, freq2)] = fg_cl_dic['cmb'][(freq2, freq1)] = cl_cmb
                if (1):#'ksz' not in ignore_fg:
                    if 'ksz' not in fg_cl_dic: fg_cl_dic['ksz'] = {}
                    fg_cl_dic['ksz'][(freq1, freq2)] = fg_cl_dic['ksz'][(freq2, freq1)] = cl_ksz
                if (1):#'tsz' not in ignore_fg:
                    if 'tsz' not in fg_cl_dic: fg_cl_dic['tsz'] = {}
                    fg_cl_dic['tsz'][(freq1, freq2)] = fg_cl_dic['tsz'][(freq2, freq1)] = cl_tsz
                if (1):#'radio' not in ignore_fg:
                    if 'radio' not in fg_cl_dic: fg_cl_dic['radio'] = {}
                    fg_cl_dic['radio'][(freq1, freq2)] = fg_cl_dic['radio'][(freq2, freq1)] = cl_radio
                if (1):#'dust' not in ignore_fg:
                    if 'cib' not in fg_cl_dic: fg_cl_dic['cib'] = {}
                    fg_cl_dic['cib'][(freq1, freq2)] = fg_cl_dic['cib'][(freq2, freq1)] = cl_dust
                if (1):#'dust' not in ignore_fg and 'tsz' not in ignore_fg and 'tsz_cib' not in ignore_fg:
                    if 'tsz_cib' not in fg_cl_dic: fg_cl_dic['tsz_cib'] = {}
                    ##fg_cl_dic['tsz_cib'][(freq1, freq2)] = fg_cl_dic['tsz_cib'][(freq2, freq1)] = cl_tsz_cib
                    fg_cl_dic['tsz_cib'][(freq1, freq2)] = cl_tsz_cib
                if (1):#'noise' not in ignore_fg and which_spec != 'TE':
                    if 'noise' not in fg_cl_dic: fg_cl_dic['noise'] = {}
                    fg_cl_dic['noise'][(freq1, freq2)] = fg_cl_dic['noise'][(freq2, freq1)] = nl
                if include_gal:
                    if 'galdust' not in fg_cl_dic: 
                        fg_cl_dic['galdust'] = {}
                        fg_cl_dic['galsync'] = {}
                    fg_cl_dic['galdust'][(freq1, freq2)] = fg_cl_dic['galdust'][(freq2, freq1)] = cl_gal_dust
                    fg_cl_dic['galsync'][(freq1, freq2)] = fg_cl_dic['galsync'][(freq2, freq1)] = cl_gal_sync

            cl[np.isnan(cl)] = 0.
            cl[np.isinf(cl)] = 0.

            if (0):##freq1 == 90 or freq2 == 90:
                print(cl[1000], cl_dust[1000], cl_tsz[1000], cl_radio[1000], cl_ksz[1000], cl_cmb[1000], cl_tsz_cib[1000], freq1, freq2)

            ##########################################################################################
            #20200516 - adjusting Nl when beam is too large (for 30/40 GHz bands)
            #removed this on 20201121 and setting nl to a large constant value
            adjust_for_large_beams = 0                   
            if adjust_for_large_beams:
                bl = bl_dic[freq1]
                if 'effective' in bl_dic:
                    bl_eff = bl_dic['effective']
                else:
                    bl_eff = bl_dic[145]
                bad_inds = np.where( (bl_eff/bl > beam_tol_for_ilc) )[0]
                print(freq1, freq2, bad_inds)
                cl[bad_inds] = 0.
            #20200516 - adjusting Nl when beam is too large (for 30/40 GHz bands)
            ##########################################################################################

            if (0):#which_spec == 'TT':
                loglog(cl); title('%s - %s' %(freq1, freq2)); show(); sys.exit()

            if (0):##which_spec == 'TE':
                #print(which_spec, freq1, freq2, cl)
                loglog(cl, label = r'%s - %s' %(freq1, freq2));
                legend(loc = 3, ncol = 2, fontsize = 6) 

            if (0):##freq1 != freq2: 
                print(freq1, freq2, cl.min(), cl.max(), nl.min(), nl.max()); sys.exit()

            cl_dic[(freq1, freq2)] = cl
            #print('\n\n\n')
    if debug: legend(loc = 1, ncol = 2, fontsize = 10); show(); sys.exit()
    if return_fg_spectra:
        return el, cl_dic, fg_cl_dic
    else:
        return el, cl_dic  

################################################################################################################
################################################################################################################
################################################################################################################

def get_acap(freqarr, final_comp = 'cmb', spec_index_rg = -0.76, freqcalib_fac = None):

    nc = len(freqarr)

    if freqcalib_fac is None: freqcalib_fac = np.ones(nc)

    if final_comp.lower() == 'cmb':
        freqscale_fac = np.ones(nc)

    elif final_comp.lower() == 'tsz' or final_comp.lower() == 'y':

        freqscale_fac = []
        for freq in sorted( freqarr ):
            freqscale_fac.append( compton_y_to_delta_Tcmb(freq * 1e9) )

        freqscale_fac = np.asarray( freqscale_fac )

    elif final_comp.lower() == 'cib' or final_comp.lower() == 'cibpo':
        freqscale_fac = []
        for freq in sorted( freqarr ):
            freqscale_fac.append( get_cib_freq_dep(freq * 1e9) )

        freqscale_fac = np.asarray( freqscale_fac )

    elif final_comp.lower() == 'cibclus':
        freqscale_fac = []
        for freq in sorted( freqarr ):
            freqscale_fac.append( get_cib_freq_dep(freq * 1e9, beta = 2.505) )

        freqscale_fac = np.asarray( freqscale_fac )

    elif final_comp.lower().find('misc_cib')>-1:

        #default values
        misc_tcib = 20.
        misc_beta = 1.505
        tcib_tmp =  re.findall('tcib\d*\.?\d+', final_comp.lower())
        if len(tcib_tmp)>0:
            tcib_tmp = tcib_tmp[0]
            misc_tcib = float(tcib_tmp.replace('tcib', ''))

        beta_tmp =  re.findall('beta\d*\.?\d+', final_comp.lower())
        if len(beta_tmp)>0:
            beta_tmp = beta_tmp[0]
            misc_beta = float(beta_tmp.replace('beta', ''))

        freqscale_fac = []
        for freq in sorted( freqarr ):
            freqscale_fac.append( get_cib_freq_dep(freq * 1e9, Tcib = misc_tcib, beta = misc_beta) )

        freqscale_fac = np.asarray( freqscale_fac )
        freqscale_fac /= np.max(freqscale_fac)

    elif final_comp.lower() == 'radio':
        freqscale_fac = []
        for freq in sorted( freqarr ):
            freqscale_fac.append( get_radio_freq_dep(freq, spec_index_rg = spec_index_rg) )

        freqscale_fac = np.asarray( freqscale_fac )

    acap = np.zeros(nc) + (freqscale_fac * freqcalib_fac) #assuming CMB is the same and calibrations factors are same for all channel

    acap = np.mat(acap).T #should be teb_len*nc x teb_len
    
    return acap

def create_clmat(freqarr, elcnt, cl_dic):
    """
    freqarr  = array of frequency channel
    elcnt = \el index
    cl_dic = cl_cmb + cl_FG auto and cross spectra of the frequency channels
    """
    nc = len(freqarr)
    clmat = np.zeros( (nc, nc) )

    for ncnt1, freq1 in enumerate(freqarr):
        for ncnt2, freq2 in enumerate(freqarr):
            clmat[ncnt2, ncnt1] = cl_dic[(freq1, freq2)][elcnt]

    return clmat

def get_clinv(freqarr, elcnt, cl_dic, which_spec):

    clmat_dic = {}
    if 'TT' in cl_dic:
        clmat_dic['TT'] = np.mat( create_clmat(freqarr, elcnt, cl_dic['TT']) )
    if 'EE' in cl_dic:
        clmat_dic['EE'] = np.mat( create_clmat(freqarr, elcnt, cl_dic['EE']) )
    if 'TE' in cl_dic:
        clmat_dic['TE'] = np.mat( create_clmat(freqarr, elcnt, cl_dic['TE']) )

    if (0):        
        if which_spec == 'TT' or which_spec == 'EE':
            clmat = clmat_dic[which_spec]
        else:
            nc = len(freqarr)
            clmat = np.zeros( ( nc * 2, nc * 2) )
            clmat[:nc, :nc] = clmat_dic['TT']
            clmat[nc:, nc:] = clmat_dic['EE']
            clmat[nc:, :nc] = clmat_dic['TE']
            clmat[:nc, nc:] = clmat_dic['TE']

            #imshow(clmat, vmin = 0., vmax = 100.); colorbar();sys.exit()
    
    if (1):
        clmat = clmat_dic[which_spec]

    if (0):##elcnt == 100:
        sbpldic = {'TT': 1, 'EE': 2, 'TE': 3}
        subplot(1,3,sbpldic[which_spec])
        imshow(clmat_dic[which_spec], vmin = None, vmax = None); colorbar(); title(which_spec);

    #clinv = sc.linalg.pinv2(clmat)
    clinv = np.linalg.pinv(clmat)
    #clinv = sc.linalg.inv(clmat)


    return clinv

def residual_power(param_dict, freqarr, el, cl_dic, which_spec, final_comp = 'cmb', spec_index_rg = -0.76, freqcalib_fac = None, lmin = 0, return_weights = 0, null_comp = None):

    try:
        lmin = param_dict['lmin']
    except:
        pass

    if (0):
        if which_spec == 'TE':
            te_len = 2
        else:
            te_len = 1

        freqarr_for_acap = np.tile(freqarr, te_len) #T and E, if need be
        acap = get_acap(freqarr_for_acap, final_comp = final_comp, spec_index_rg = spec_index_rg, freqcalib_fac = freqcalib_fac)

        nc = len(freqarr)
        cl_residual = np.zeros( (len(el)) )
        weightsarr = np.zeros( (te_len * nc, len( el ) ) )
    
    if (1):
        acap = get_acap(freqarr, final_comp = final_comp, spec_index_rg = spec_index_rg, freqcalib_fac = freqcalib_fac)

        nc = len(freqarr)
        cl_residual = np.zeros( (len(el)) )
        weightsarr = np.zeros( (nc, len( el ) ) )

    for elcnt, el in enumerate(el):
        if el <= lmin: continue ## or el>=lmax: continue
        clinv = get_clinv( freqarr, elcnt, cl_dic, which_spec )
        if (0):#el == lmin+1:
            print(which_spec, clinv.shape)

        nr = np.dot(clinv, acap)
        dr = np.dot( acap.T, np.dot(clinv, acap) )
        #drinv = sc.linalg.pinv2(dr)
        drinv = np.linalg.pinv(dr)
        weight = np.dot(nr, drinv)

        print(weight, el)
        sys.exit()

        #ILC residuals
        cl_residual[elcnt] = np.asarray(1./dr).squeeze()
        weightsarr[:, elcnt] = weight.squeeze()


    weightsarr = np.asarray( weightsarr )
    cl_residual = np.asarray( cl_residual )
    print(weightsarr.shape)

    #from IPython import embed; embed()

    cl_residual[np.isinf(cl_residual)] = 0.
    cl_residual[np.isnan(cl_residual)] = 0.
    
    if return_weights:
        return cl_residual, weightsarr
    else:
        return cl_residual


################################################################################################################
################################################################################################################
################################################################################################################
def get_effective_frequencies(experiment, band, component):
    """
    Function to get the SPT-SZ/SPT-SZ+SPTpol/SPT-3G/Planck effective freqeuncies for
    different foreground component.

    Parameters
    ----------
    experiment: str
        spt experiment for which we need the effective frequencies.
        can be one of: 'george', sptsz', 'reichardt', 'sptszpol', 'stp3g', '3g', 'planck'
        george/sptsz - Sec 5.7 of https://arxiv.org/pdf/1408.3161.pdf
        reichardt/sptszpol - Sec 5.1 of https://arxiv.org/pdf/2002.06197.pdf
        spt3g/3g - The effective band centres are calculated by Zhaodi
        (https://pole.uchicago.edu/spt3g/images/Effective_band_center_for_3G.pdf).
        Planck - calculated by Tom C. and included here by Srini.
            spt.uchicago.edu:/home/tcrawfor/spt_analysis/temp/nu_eff_planck_for_sr_batch_03may21.pro
        Galactic dust (GD): scratch/wlwu/script21_0413_solve_nueff_galdust.py
    band: str
        95GHz or 150GHz or 220GHz (SPT)
        100GHz or 143GHz or 217GHz or 353GHz or 545GHz or 857GHz (Planck)
        for which we need the effective band centre
    component : str
        The foreground component to use. Must be one of
        'tSZ', 'DG-Cl', 'DG-Po', 'DG', 'RG','GD'
    Returns
    -------
    band_centre: float for the desired foreground component.
    """
    band = '%sGHz' %(band)
    experiment = experiment.lower()
    component = component.lower()

    if component == 'cmb' or component == 'ksz':
        return band

    spt_exps = ['george', 'sptsz', 'reichardt', 'sptszpol', 'spt3g', '3g', 'spt4', 'spt4_c3', 'spt4_c4', 'spt3g+', 'spt3g_201920', 'spt3g_201920_plancksevenbands']
    cmb_exps = ['s4wide', 's4deepv3r025', 's4deep', 'cmbs4']
    assert experiment.lower() in spt_exps or experiment.lower() == 'planck' or experiment.lower() in cmb_exps
    assert component in ['tsz', 'dg-cl', 'dg-po', 'dg', 'rg', 'y']

    if band == '90GHz':
        band = '95GHz'

    if experiment in spt_exps:
        if band not in ['95GHz', '150GHz', '220GHz', '221GHz', '286GHz', '345GHz']:
            raise ValueError("band must be one of 90/95, 150, 220GHz for SPT")
        if component == 'dg-cl' or component == 'dg-po' or component == 'dg':
            component = 'dg'

    if experiment == 'planck':
        if band not in ['100GHz', '143GHz', '217GHz', '353GHz', '545GHz', '857GHz']:
            raise ValueError(
                "band must be one of 100, 143, 217, 353, 545, 857GHz for Planck"
            )
        if component == 'dg':
            component = 'dg-po'

    if experiment == 'george' or experiment == 'sptsz':
        eff_frequencies = {
            'dg': {'95GHz': 96.9, '150GHz': 153.4, '220GHz': 221.6},
            'rg': {'95GHz': 93.5, '150GHz': 149.5, '220GHz': 215.8},
            'tsz': {'95GHz': 96.6, '150GHz': 152.3, '220GHz': 220.1},
        }
    elif experiment == 'reichardt' or experiment == 'sptszpol':
        eff_frequencies = {
            'dg': {'95GHz': 96.9, '150GHz': 153.4, '220GHz': 221.6},
            'rg': {'95GHz': 93.5, '150GHz': 149.5, '220GHz': 215.8},
            'tsz': {'95GHz': 96.6, '150GHz': 152.3, '220GHz': 220.1},
        }
    elif experiment == 'spt3g' or experiment == '3g' or experiment == 'spt3g_201920':
        eff_frequencies = {
            'dg': {'95GHz': 95.96, '150GHz': 150.01, '220GHz': 222.76},
            'rg': {'95GHz': 93.52, '150GHz': 145.92, '220GHz': 213.34},
            'tsz': {'95GHz': 95.69, '150GHz': 148.85, '220GHz': 220.15},
            'gd':{'95GHz': 95.4791, '150GHz': 149.5289, '220GHz': 222.3065},
        }
    elif experiment.find('s4')>-1:
        eff_frequencies = {
            'tsz': {'93GHz': 95.90933019795854, '145GHz': 145.7404333477917, '155GHz': 155.78797681444541, '225GHz': 220.1482522978387, '278GHz': 269.2497424257618},
        }
    elif experiment.find('spt4')>-1 or experiment.find('spt3g+')>-1:
        eff_frequencies = {
            'tsz': {'95GHz': 95.69, '150GHz': 148.85, '220GHz': 220.15, '221GHz': 221.16065337121762, '286GHz': 285.25585052369, '345GHz': 344.61575620337703},
        }
    elif experiment == 'planck':
        eff_frequencies = {
            'dg-cl': {
                '100GHz': 103.471,
                '143GHz': 145.689,
                '217GHz': 224.882,
                '353GHz': 362.999,
                '545GHz': 552.200,
                '857GHz': 833.825,
            },
            'dg-po': {
                '100GHz': 103.194,
                '143GHz': 145.319,
                '217GHz': 224.399,
                '353GHz': 362.200,
                '545GHz': 550.585,
                '857GHz': 831.285,
            },
            'rg': {
                '100GHz': 100.895,
                '143GHz': 141.513,
                '217GHz': 216.254,
                '353GHz': 360.852,
                '545GHz': 543.832,
                '857GHz': 820.825,
            },
            'tsz': {
                '100GHz': 102.942,
                '143GHz': 144.423,
                '217GHz': 234.162,
                '353GHz': 355.093,
                '545GHz': 529.123,
                '857GHz': 784.019,
            },
        }
    eff_frequencies['y'] = eff_frequencies['tsz']

    return eff_frequencies[component][band]

def get_acap_new(freqarr, final_comp = 'cmb', spec_index_rg = -0.76, freqcalib_fac = None, teb_len = 1, experiment = None):

    #if experiment is not None and ( final_comp.lower() != 'cmb' and final_comp.lower() != 'ksz'):
    if experiment is not None and ( final_comp.lower() != 'cmb' and final_comp.lower() != 'ksz') and experiment.find('sptpolplusultradeepplus3gplusherschel') == -1 and experiment.find('spt3g_201920_planck') == -1 and experiment.find('ccat') == -1:
        if final_comp.lower() == 'tsz' or final_comp.lower() == 'y': #20220614 - only for tSZ for now.
            freqarr_mod = []
            for freq in sorted( freqarr ):
                freqarr_mod.append( get_effective_frequencies(experiment, freq, final_comp) )
            freqarr = np.copy(freqarr_mod)

    nc = len(freqarr)

    if freqcalib_fac is None: freqcalib_fac = np.ones(nc)

    if final_comp.lower() == 'cmb':
        freqscale_fac = np.ones(nc)

    elif final_comp.lower() == 'tsz' or final_comp.lower() == 'y':

        freqscale_fac = []
        for freq in sorted( freqarr ):
            freqscale_fac.append( compton_y_to_delta_Tcmb(freq * 1e9) )

        freqscale_fac = np.asarray( freqscale_fac )

    elif final_comp.lower() == 'tsz_cib' or final_comp.lower() == 'cib_tsz':
        pass
        '''
        tsz_cib_dic = {90: -0.0160040465551, 93: -0.0157431190535, 145: -0.0250242003427, 150: -0.0234920051701, 220: 0.00239951463962, 225: 0.00780531042154, 278: 0.157247256562, 345:1.0}
        freqscale_fac = []
        for freq in sorted( freqarr ):
            freqscale_fac.append( tsz_cib_dic[freq] )

        freqscale_fac = np.asarray( freqscale_fac )
        freqscale_fac /= np.max(freqscale_fac)
        '''

    elif final_comp.lower() == 'cib' or final_comp.lower() == 'cibpo':
        freqscale_fac = []
        for freq in sorted( freqarr ):
            freqscale_fac.append( get_cib_freq_dep(freq * 1e9) )

        freqscale_fac = np.asarray( freqscale_fac )
        freqscale_fac /= np.max(freqscale_fac)
        
    elif final_comp.lower() == 'cibclus':
        freqscale_fac = []
        for freq in sorted( freqarr ):
            freqscale_fac.append( get_cib_freq_dep(freq * 1e9, beta = 2.505) )

        freqscale_fac = np.asarray( freqscale_fac )
        freqscale_fac /= np.max(freqscale_fac)

    elif final_comp.lower().find('misc_cib')>-1:

        #default values
        misc_tcib = 20.
        misc_beta = 1.505
        tcib_tmp =  re.findall('tcib\d*\.?\d+', final_comp.lower())
        if len(tcib_tmp)>0:
            tcib_tmp = tcib_tmp[0]
            misc_tcib = float(tcib_tmp.replace('tcib', ''))

        beta_tmp =  re.findall('beta\d*\.?\d+', final_comp.lower())
        if len(beta_tmp)>0:
            beta_tmp = beta_tmp[0]
            misc_beta = float(beta_tmp.replace('beta', ''))

        #freqarr = [30, 44, 70, 100, 150, 217, 353, 545]
        freqscale_fac = []
        for freq in sorted( freqarr ):
            freqscale_fac.append( get_cib_freq_dep(freq * 1e9, Tcib = misc_tcib, beta = misc_beta) )

        freqscale_fac = np.asarray( freqscale_fac )
        freqscale_fac /= np.max(freqscale_fac)

        if (0):
            norm_freq = 150
            freqarr = np.asarray( freqarr )
            norm_freq_ind = np.argmin( abs(freqarr - norm_freq) )
            freqscale_fac /= freqscale_fac[norm_freq_ind]
            
            '''
            loglog(freqarr, freqscale_fac, 'go'); xlim(10., 1e3); ylim(0.05, 400.); show()
            print( freqscale_fac )
            sys.exit()
            '''

    elif final_comp.lower() == 'radio':
        freqscale_fac = []
        for freq in sorted( freqarr ):
            freqscale_fac.append( get_radio_freq_dep(freq, spec_index_rg = spec_index_rg) )

        freqscale_fac = np.asarray( freqscale_fac )

    acap = np.zeros(nc) + (freqscale_fac * freqcalib_fac) #assuming CMB is the same and calibrations factors are same for all channel

    if teb_len>1:
        acap_full = np.zeros( (teb_len, len(acap) * teb_len) )
        #acap_full[:len(acap), :] = acap

        #acap = np.mat(np.eye(len(acap)) * acap)
        acap_full[0,:len(acap)] = acap
        if final_comp.lower() == 'cmb':
            acap_full[1,len(acap):] = acap
        else: #polarisation weights are zero for other foregrounds
            acap_full[1,len(acap):] = 0.

        '''
        acap_full[0,len(acap):] = acap/10.
        acap_full[1,:len(acap)] = acap/10.
        '''

        acap_full = np.mat(acap_full).T #should be teb_len*nc x teb_len
        acap = acap_full
    else:
        acap = np.mat(acap).T #should be teb_len*nc x teb_len
    
    return acap

def get_cib_freq_dep(nu, Tcib = 20., Tcmb = 2.7255, h=6.62607004e-34, k_B=1.38064852e-23, beta = 1.505):
    if nu<1e4: nu *= 1e9

    bnu1 = fg.fn_BnuT(nu, temp = Tcib)
    dbdt = fg.fn_dB_dT(nu)
    value = (nu**beta) * bnu1 / dbdt

    return value

#def get_radio_freq_dep(nu, nu0 = 150., spec_index_rg = -0.9, null_highfreq_radio = 1):
def get_radio_freq_dep(nu, nu0 = 150., spec_index_rg = -0.76, null_highfreq_radio = 1):

    nr = fg.fn_dB_dT(nu0)
    dr = fg.fn_dB_dT(nu)
    epsilon_nu1_nu0 = nr/dr
    scaling = (nu/nu0)**spec_index_rg
    value = epsilon_nu1_nu0 * scaling

    if null_highfreq_radio and (nu>230):
        #print('\n\tthis extrapolation does not work for high freqeuncy radio. Making cl_radio = 0 for these bands.')
        value = 0.

    return value

def get_teb_spec_combination(cl_dic):
    pspec_arr = sorted( list( cl_dic.keys() ) )

    if pspec_arr == ['TT'] or pspec_arr == ['EE']: #only TT is supplied
        teb_len = 1
    elif pspec_arr == sorted(['TT', 'EE']) or pspec_arr == sorted(['TT', 'EE', 'TE']): #TT/EE/TE are supplied
        teb_len = 2
    elif pspec_arr == sorted(['TT', 'EE', 'BB']) or pspec_arr == sorted(['TT', 'EE', 'BB', 'TE', 'TB', 'EB']): #TT/EE/BB are supplied
        teb_len = 3
    else:
        teb_len = 1
        pspec_arr = None

    return teb_len, pspec_arr

def fn_corr_from_cov(covmat):
    diags = np.diag(covmat)
    corrmat = np.zeros_like(covmat)
    for i in range(covmat.shape[0]):
        for j in range(covmat.shape[0]):
            corrmat[i, j] = covmat[i, j] / np.sqrt( diags[i] *  diags[j] )
    #corrmat = covmat / np.outer(np.sqrt(diags), np.sqrt(diags))
    return corrmat

def create_clmat_new(freqarr, elcnt, cl_dic):
    """
    freqarr  = array of frequency channel
    elcnt = \el index
    cl_dic = cl_cmb + cl_FG auto and cross spectra of the frequency channels
    """
    nc = len(freqarr)
    teb_len, pspec_arr = get_teb_spec_combination(cl_dic)
    clmat = np.zeros( (teb_len * nc, teb_len * nc) )

    for pspecind, pspec in enumerate( pspec_arr ):
        curr_cl_dic = cl_dic[pspec]
        '''
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
        '''
        if teb_len == 1: #clmat for TT or EE or BB
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
        

    if (0): #get correlation matrix

        ##clmat = clmat + np.eye(len(clmat))+1.
        corr_mat = np.mat( fn_corr_from_cov( clmat ) )
        lower_tril = np.tril(corr_mat, k = -1)
        maxval = np.max(abs(lower_tril))
        if elcnt == 1000: ##maxval>=1.: 
            from IPython import embed; embed()
            #print('Nulling cov for this \ell = %s; max(corr) is %s' %(elcnt, maxval))
            #clmat *= 0.

    return clmat

def get_clinv_new(freqarr, elcnt, cl_dic, return_clmat = 0):
    clmat = np.mat( create_clmat_new(freqarr, elcnt, cl_dic) )
    #clmat = clmat + np.eye(len(clmat)) * np.min(np.diag(clmat)/1e3)
    #clinv = sc.linalg.pinv2(clmat)
    #clinv = sc.linalg.inv(clmat)
    clinv = np.linalg.pinv(clmat)

    if (0):##elcnt == 1032:
        from IPython import embed; embed()
        print(elcnt, clinv)
        #sys.exit()

    if return_clmat:
        return clinv, clmat
    else:
        return clinv

def residual_power_new(param_dict, freqarr, el, cl_dic, final_comp = 'cmb', spec_index_rg = -0.76, freqcalib_fac = None, lmin = 10, return_weights = 0, null_comp = None, experiment = None):

    #acap = get_acap(freqarr, final_comp = final_comp, freqcalib_fac = freqcalib_fac)
    teb_len, pspec_arr = get_teb_spec_combination(cl_dic) #20200527 - teb
    acap = get_acap_new(freqarr, final_comp = final_comp, spec_index_rg = spec_index_rg, freqcalib_fac = freqcalib_fac, teb_len = teb_len, experiment = experiment)

    if null_comp is not None:
        total_comp_to_null = 0
        if np.ndim(null_comp) == 0:
            bcap = get_acap_new(freqarr, final_comp = null_comp, spec_index_rg = spec_index_rg, freqcalib_fac = freqcalib_fac, teb_len = teb_len, experiment = experiment)
            total_comp_to_null += 1
        else:
            #from IPython import embed; embed()
            bcap = None
            for curr_null_comp in null_comp:
                curr_bcap = get_acap_new(freqarr, final_comp = curr_null_comp, spec_index_rg = spec_index_rg, freqcalib_fac = freqcalib_fac, teb_len = teb_len, experiment = experiment)
                if bcap is None:
                    bcap = np.copy( curr_bcap )
                else:
                    bcap = np.column_stack( (bcap, curr_bcap) )
                total_comp_to_null += 1
            bcap = np.mat(bcap)
            bcap_full = np.mat( np.copy( bcap ) )
            #from IPython import embed; embed()
            #bcap = get_acap_new(freqarr, final_comp = null_comp[0], spec_index_rg = spec_index_rg, freqcalib_fac = freqcalib_fac, teb_len = teb_len)

    nc = len(freqarr)
    #cl_residual = np.zeros( (len(el)) )
    #weightsarr = np.zeros( (nc, len( el ) ) )
    weightsarr = np.zeros( (teb_len * nc, teb_len, len( el ) ) )
    cl_residual = np.zeros( (3, len(el)) )

    cl_residual_tmp = []

    for elcnt, currel in enumerate(el):
        if currel <= lmin: continue ## or el>=lmax: continue
        #clinv = get_clinv_new( freqarr, elcnt, cl_dic )
        clinv, clmat = get_clinv_new( freqarr, elcnt, cl_dic, return_clmat = 1 )

        '''
        if (0):##el>200:
            from IPython import embed; embed()
            clf()
            tt_cov = clmat[:nc, :nc]
            ee_cov = clmat[nc:, nc:]
            te_cov = clmat[nc:, :nc]
            subplot(131); imshow(tt_cov); colorbar(); title(r'TT')
            subplot(132); imshow(ee_cov); colorbar(); title(r'EE')
            subplot(133); imshow(te_cov); colorbar(); title(r'TE')
            show()

            tt_var, ee_var, te_var = np.diag(tt_cov), np.diag(ee_cov), np.diag(te_cov)
            te_var_predction = 0.35 * np.sqrt( tt_var * ee_var )
            ax = subplot(111, yscale = 'log')
            plot(freqarr, tt_var, 'k'); plot(freqarr, ee_var, 'green'); plot(freqarr, te_var, 'red'); 
            plot(freqarr, te_var_predction, 'red', ls = '--', lw = 2.); 
            show()
            medval = np.median( np.asarray(clmat) )
            imshow(clmat, vmin = medval * 0.5, vmax = medval * 2.); colorbar(); title(r'$\ell$ = %s' %(el)); show()
        '''
        if null_comp is None:
            nr = np.dot(clinv, acap)
            dr = np.dot( acap.T, np.dot(clinv, acap) )

            #drinv = sc.linalg.pinv2(dr)
            drinv = np.linalg.pinv(dr)
            weight = np.dot(nr, drinv)
        else:

            '''
            G = np.column_stack( (acap, bcap) )
            ncap = np.mat([1., 0.]).T
            '''

            if (0):
                if currel>=2000:
                    bcap = bcap_full#[:,0]
                else:
                    bcap = bcap_full[:,1]

            G = np.column_stack( (acap, bcap) )
            G = np.mat(G)

            total_comps = G.shape[1]
            ###total_comps = int( G.shape[1]/teb_len ) #20210506 - check me
            ncap = np.zeros( total_comps )#total_comp_to_null + 1 )
            ncap[0] = 1.
            #ncap[1] = 0.5
            ncap = np.mat( ncap ).T

            nr = np.dot(clinv, G)
            dr = np.dot( G.T, np.dot(clinv, G) )

            #drinv = np.dot( sc.linalg.pinv2(dr), ncap)
            drinv = np.dot( np.linalg.pinv(dr), ncap)
            weight = np.dot(nr, drinv)

            if (0):#currel>1500:
                #from IPython import embed; embed()
                acap_sum = np.sum( np.asarray(acap) * np.asarray(weight) )
                bcap_1_sum = np.sum( np.asarray(bcap[:, 0]) * np.asarray(weight) )
                print(currel, acap_sum, bcap_1_sum)
                '''
                if total_comp_to_null == 2:
                    bcap_2_sum = np.sum( np.asarray(bcap[:, 1]) * np.asarray(weight) )
                    print(weight, G, ncap, acap_sum, bcap_1_sum, bcap_2_sum)
                elif total_comp_to_null == 3:
                    bcap_2_sum = np.sum( np.asarray(bcap[:, 1]) * np.asarray(weight) )
                    bcap_3_sum = np.sum( np.asarray(bcap[:, 2]) * np.asarray(weight) )
                    print(weight, G, ncap, acap_sum, bcap_1_sum, bcap_2_sum, bcap_3_sum)
                else:
                    print(weight, G, ncap, acap_sum, bcap_1_sum)
                sys.exit()
                '''

        #ILC residuals
        #cl_residual[elcnt] = np.asarray(1./dr).squeeze()
        ###from IPython import embed; embed(); sys.exit()
        if teb_len>1:
            cl_residual_tt, cl_residual_ee, cl_residual_te = drinv[0,0], drinv[1,1], drinv[0,1]
            cl_residual[:, elcnt] = cl_residual_tt, cl_residual_ee, cl_residual_te
        else:
            cl_residual[0, elcnt] = drinv[0]


        #weightsarr[:, elcnt] = np.asarray(nr/dr).squeeze()
        weightsarr[:, :, elcnt] = weight

        cl_residual_tmp.append( drinv )

    '''
    if (1):
        cl_residual_tmp = np.asarray( cl_residual_tmp )
        tt, ee, te = cl_residual_tmp[:,0,0], cl_residual_tmp[:,1,1], cl_residual_tmp[:,0,1]
        ax  =subplot(111, yscale = 'log')
        plot(tt, 'k'); plot(ee, 'green'); plot(te, 'red'); ylim(1e-15, 1e2); show(); 

        from IPython import embed; embed()
        clf()
        weightsarr = np.asarray( weightsarr )
        tt_w, ee_w = np.asarray( [weightsarr[:nc, 0], weightsarr[nc:, 1]] )
        te_tt_w, te_ee_w = np.asarray( [weightsarr[nc:, 0], weightsarr[:nc, 1]] )
        w_for_plot = [tt_w, ee_w]
        colordic = {27:'indigo', 39:'royalblue', 93: 'lightgreen', 145: 'darkgreen', 225: 'goldenrod', 278: 'darkred'}
        for iter in range(2):
            ax  =subplot(3,1,iter+1); 
            for i in range(nc):
                nu = freqarr[i]
                plot(w_for_plot[iter][i], color = colordic[nu], label = r'%s' %nu);
            axhline(0., lw = 0.1)
            axhline(1., lw = 0.1)
            plot(np.sum(w_for_plot[iter], axis=0), 'k--')
            if iter == 0:
                legend(loc = 4, ncol = 5, fontsize = 4)

        ax = subplot(3,1,3)
        te_w_for_plot = [te_tt_w, te_ee_w]
        for iter in range(2):
            if iter == 0:
                ls = '-'
            else:
                ls = '-.'
            for i in range(nc):
                nu = freqarr[i]
                #plot(te_w_for_plot[iter][i], color = colordic[nu], label = r'%s' %nu);
            axhline(0., lw = 0.1)
            axhline(1., lw = 0.1)
            plot(np.sum(te_w_for_plot[iter], axis=0), color = 'k', ls = ls )
        show()

        dummyel = 200
        tmp = weightsarr[:,:,dummyel]
        print(tmp)
        print(tt_w[:,dummyel], np.sum(tt_w[:,dummyel]))
        print(ee_w[:,dummyel], np.sum(ee_w[:,dummyel]))
        print(te_tt_w[:,dummyel], np.sum(te_tt_w[:,dummyel]))
        print(te_ee_w[:,dummyel], np.sum(te_ee_w[:,dummyel]))

    sys.exit()
    '''

    ###from IPython import embed; embed()

    weightsarr = np.asarray( weightsarr )
    cl_residual = np.asarray( cl_residual )

    #from IPython import embed; embed()

    cl_residual[np.isinf(cl_residual)] = 0.
    cl_residual[np.isnan(cl_residual)] = 0.
    
    if return_weights:
        return cl_residual, weightsarr
    else:
        return cl_residual

################################################################################################################

def coth(x):
    return (np.exp(x) + np.exp(-x)) / (np.exp(x) - np.exp(-x))

################################################################################################################

#def compton_y_to_delta_Tcmb(nu, Tcmb = 2.73, h=6.62607004e-34, k_B=1.38064852e-23):
def compton_y_to_delta_Tcmb(nu, Tcmb = 2.73, h=6.626e-34, k_B=1.38e-23):

    if nu<1e3: nu *= 1e9

    x = h * nu / k_B / Tcmb
    g_nu = x * coth(x/2.) - 4.

    return Tcmb * np.mean(g_nu)

################################################################################################################

def get_ilc_map(final_comp, el, map_dic, bl_dic, nside, lmax, cl_dic = None, spec_index_rg = -0.76, nl_dic = None, lmin = 10, freqcalib_fac = None, ignore_fg = [], full_sky = 0, estimate_covariance = 0, mapparams = None, apod_mask = None, weightsarr = None):

    """
    inputs:
    final_comp: 'CMB' or 'tsz' or 'y'

    map_dic: dicionary containing
    map_dic[95] = map_95 where map_95 is either a healpix map or a SPT3G map object or a simple flat-sky 2D numpy array
    map_dic[150] = map_150 
    map_dic[220] = map_220

    bl_dic = {} #dictionary containing Bl of different frequencies - either a 1d array for healpix maps or 2D array for SPT3G/flat-sky maps
    bl_dic[95], bl_dic[150], bl_dic[220]

    nside = healpix nside set to None for flat-sky maps

    lmax = lmax used for ILC (covariance)

    cl_dic = same as beam dic and contains Cl of different frequencies - either a 1d array for healpix maps or 2D array for SPT3G/flat-sky maps
    nl_dic = same as beam dic and contains Nl of different frequencies - either a 1d array for healpix maps or 2D array for SPT3G/flat-sky maps

    lmin = minimum \el to use

    freqcalib_fac = array containing calibration factors for different frequencies. for example: freqcalib_fac = [1., 1., 1.]

    """
    spt3g_maps = 0
    freqarr, maparr = [], []
    for keyname in sorted(map_dic):
        freqarr.append( keyname )
        curr_map = map_dic[keyname]
        '''
        if isinstance(curr_map, core.G3Frame):
            spt3g_maps = 1
            #get the apodisation mask before converting the map inot an array
            if str(apod_mask) == 'from_weight':
                apod_mask = apodmask.make_border_apodization(curr_map, radius_arcmin=90.0)
            curr_map = np.array( curr_map['T'] ) ## / core.G3Units.uK
            curr_map[np.isnan(curr_map)] = 0.
            curr_map[np.isinf(curr_map)] = 0.
        '''
        if apod_mask is not None:
            curr_map = curr_map * apod_mask

        maparr.append( curr_map )

    #get covariance
    if cl_dic is None:
        '''
        if spt3g_maps:
            el, cl_dic = get_spt3g_covariance_dic(map_dic, lmin, lmax)
        else:
            el, cl_dic = get_analytic_covariance(freqarr, nl_dic = nl_dic, ignore_fg = ignore_fg)
        '''
        if estimate_covariance:
            el, cl_dic = get_map_covariance(map_dic, lmin, lmax) 
        else:
            el, cl_dic = get_analytic_covariance(freqarr, nl_dic = nl_dic, ignore_fg = ignore_fg)


    #get weights
    if weightsarr is None:
        weightsarr = get_multipole_weightsarr(final_comp, freqarr, el, cl_dic, lmin, freqcalib_fac, spec_index_rg = spec_index_rg)#, ignore_fg)
    weightsarr_1d = np.copy(weightsarr)

    #convert weights to 2D if flat-sky
    if not full_sky:
        assert mapparams is not None
        weightsarr_2D = []
        for currW in weightsarr:
            el = np.arange(len(currW))
            currW_2D = flatsky.cl_to_cl2d(el, currW, mapparams) 
            weightsarr_2D.append(currW_2D)
        weightsarr = np.asarray( weightsarr_2D )

    #rebeaming
    rebeamarr = misc.rebeam( bl_dic )

    #modify weights to include rebeam
    weightsarr = weightsarr * rebeamarr
    '''
    plot(weightsarr[0], 'k-'); plot(weightsarr[1], 'r-'); plot(weightsarr[2], 'g-')
    weightsarr = weightsarr * rebeamarr
    plot(weightsarr[0], 'k--'); plot(weightsarr[1], 'r-'); plot(weightsarr[2], 'g-')
    show(); sys.exit()
    sys.exit()
    '''

    #get ilc map now
    ilc_map = apply_ilc_weightsarr(maparr, weightsarr, nside, lmax, full_sky = full_sky)

    if apod_mask is not None:
        ilc_map = ilc_map * apod_mask

    return ilc_map, weightsarr_1d

################################################################################################################

def apply_ilc_weightsarr(maparr, weightsarr, nside, lmax, full_sky = 0, verbose = 0):

    '''
    clf()
    freqs = [95, 150]#, 220]
    colordic = {95: 'darkblue', 150: 'green', 220: 'darkred'}
    for frqcntr, freq in enumerate( freqs ):
        plot(weightsarr[frqcntr], color = colordic[freq], label = r'%s' %(freq))
    plot(np.sum(weightsarr, axis = 0), 'k', label = r'Sum')
    legend(loc = 1);show(); sys.exit()
    '''

    #get the ilc combined map now
    weighted_almarr = []
    for mm in range(len(maparr)):
        if full_sky:
            map_alm = H.map2alm(maparr[mm], lmax = lmax)
            curr_weight = weightsarr[mm][:lmax]
            
            map_alm_weighted = H.almxfl(map_alm, curr_weight)

            '''
            map_alm_weighted = np.zeros_like(map_alm)
            for el in range(len(curr_weight)):
                alm_inds = H.Alm.getidx(lmax, el, np.arange(el + 1))
                map_alm_weighted[alm_inds] =  curr_weight[el] * map_alm[alm_inds]
            '''
            #plot(map_alm_weighted.real, 'r');show();sys.exit()
            ###map_weighted = H.alm2map(map_alm_weighted, nside = nside, verbose = verbose, lmax = lmax)
            #clf();plot(curr_weight); show()
            #H.mollview(maparr[mm], sub = (1,2,1)); H.mollview(map_weighted, sub = (1,2,2)); show()
        else:
            curr_map = maparr[mm]
            weightsarr[mm][np.isnan(weightsarr[mm])]=0.
            weightsarr[mm][np.isinf(weightsarr[mm])]=0.
            map_alm_weighted = np.fft.fft2(curr_map) * weightsarr[mm]
        weighted_almarr.append(map_alm_weighted)

    ilc_alm = np.sum(weighted_almarr, axis = 0)

    if full_sky:
        ilc_map = H.alm2map(ilc_alm, nside = nside, verbose = verbose, lmax = lmax)
    else:
        ilc_map = np.fft.ifft2( ilc_alm ).real

    return ilc_map

################################################################################################################

def get_multipole_weightsarr(final_comp, freqarr, el, cl_dic, lmin, freqcalib_fac, spec_index_rg = -0.76):#, ignore_fg):

    #acap = get_acap(freqarr, final_comp = final_comp, freqcalib_fac = freqcalib_fac)
    teb_len, pspec_arr = get_teb_spec_combination(cl_dic) #20200527 - teb
    acap = get_acap(freqarr, final_comp = final_comp, spec_index_rg = spec_index_rg, freqcalib_fac = freqcalib_fac, teb_len = teb_len)

    nc = len(freqarr)

    #get weightsarr
    if np.ndim(el) == 1:
        #weightsarr = np.zeros( (nc, len( el ) ) )
        weightsarr = np.zeros( (teb_len * nc, len( el ) ) )
        for elcnt, curr_el in enumerate( el ):
            if curr_el <= lmin: continue ## or el>=lmax: continue
            clinv = get_clinv( freqarr, elcnt, cl_dic )
            
            nr = np.dot(clinv, acap)
            dr = np.dot( acap.T, np.dot(clinv, acap) )

            weight_mat = np.asarray(nr/dr).squeeze()

            weightsarr[:, elcnt] = weight_mat        
    else:
        #2d- PSDs
        from IPython import embed; embed()
        sys.exit()
        lx, ly = el
        for elcnt1, curr_lx in enumerate( lx ):
            for elcnt2, curr_ly in enumerate( ly ):
                curr_el = np.sqrt(curr_lx**2. + curr_ly**2.)
                if curr_el <= lmin: continue ## or el>=lmax: continue
                clinv = get_clinv( freqarr, elcnt, cl_dic )
            
            nr = np.dot(clinv, acap)
            dr = np.dot( acap.T, np.dot(clinv, acap) )

            weight_mat = np.asarray(nr/dr).squeeze()

            weightsarr[:, elcnt] = weight_mat

    return weightsarr

################################################################################################################
def get_map_covariance(map_dic, lmax, bl_dic, apod_mask = None):#, lmin = 2):

    print('Estimating covariance from maps now')

    freqarr = sorted(map_dic.keys())
    cl_dic = {}
    for cntr1, freq1 in enumerate( freqarr ):
        for cntr2, freq2 in enumerate( freqarr ):
            if (freq2, freq1) in cl_dic:
                cl_dic[(freq1, freq2)] = cl_dic[(freq2, freq1)]
                continue
            print((freq2, freq1))
            map1, map2 = map_dic[freq1], map_dic[freq2]

            cl = H.anafast(map1, map2, lmax = lmax - 1)# lmax + lmin - 1)
            blsquare = bl_dic[freq1] * bl_dic[freq2]
            cl /= blsquare
            #cl = cl[lmin: lmax]

            #account for the mask
            if apod_mask is not None:
                fsky = np.mean(apod_mask)**2.
                cl /= fsky

            cl[np.isnan(cl)] = 0.
            cl[np.isinf(cl)] = 0.
            cl_dic[(freq1, freq2)] = cl 
        
            #el = np.arange(lmin, lmax + lmin)
            el = np.arange(lmax)

    return el, cl_dic

################################################################################################################
'''
def get_spt3g_covariance_dic(map_dic, lmin, lmax, apod_mask = 'from_weight', return_2d = 1):

    freqarr = sorted(map_dic.keys())
    cl_dic = {}
    for cntr1, freq1 in enumerate( freqarr ):
        for cntr2, freq2 in enumerate( freqarr ):
            if (freq2, freq1) in cl_dic:
                cl_dic[(freq1, freq2)] = cl_dic[(freq2, freq1)]
                continue
            map1, map2 = map_dic[freq1], map_dic[freq2]

            if return_2d:
                mapps = map_analysis.calculate_powerspectra(map1, map2, apod_mask = apod_mask, return_2d = 1)['TT']

                mapps[np.isnan(mapps)] = 0.
                mapps[np.isinf(mapps)] = 0.
                cl_dic[(freq1, freq2)] = mapps 

                from IPython import embed; embed()

                el = mapps.get_lxly()

            else:
                cl = map_analysis.calculate_powerspectra(map1, map2, lmin = lmin, lmax = lmax, apod_mask = apod_mask)['TT']
                cl = np.concatenate( (np.zeros(lmin), cl) )
                cl[np.isnan(cl)] = 0.
                cl[np.isinf(cl)] = 0.
                cl_dic[(freq1, freq2)] = cl 
            
                el = np.arange(lmin, lmax)
                el = np.concatenate( (np.zeros(lmin), el) )

    return el, cl_dic
'''

################################################################################################################
