import numpy as np, sys, os, scipy as sc, healpy as H, foregrounds as fg, misc
from pylab import *
################################################################################################################
def get_analytic_covariance(param_dict, freqarr, nl_dic = None, bl_dic = None, ignore_fg = [], which_spec = 'TT', pol_frac_per_cent_dust = 0.02, pol_frac_per_cent_radio = 0.03, pol_frac_per_cent_tsz = 0., pol_frac_per_cent_ksz = 0., include_gal = 0, beam_tol_for_ilc = 1000., cib_corr_coeffs = None):

    #ignore_fg = foreground terms that must be ignored
    possible_ignore_fg = ['cmb', 'tsz', 'ksz', 'radio', 'dust']
    if len(ignore_fg)>0:
        if not all( [ currfg in possible_ignore_fg for currfg in ignore_fg] ):
            print( '\n\t Alert: Elements of ignore_fg should be one of the following: %s\n\n' %(np.array2string(np.asarray(possible_ignore_fg))) )
            sys.exit()

    el, cl_cmb = fg.get_foreground_power_spt('CMB', freq1 = param_dict['freq0'], freq2 = param_dict['freq0'])
    el, cl_ksz = fg.get_foreground_power_spt('kSZ', freq1 = param_dict['freq0'], freq2 = param_dict['freq0'])
    if which_spec == 'EE':
        cl_ksz = cl_ksz * pol_frac_per_cent_ksz**2.
    if which_spec == 'TE':
        cl_ksz = cl_ksz * 0.

    cl_dic = {}
    cl_ori = np.zeros(len(el))
    for freq1 in freqarr:
        for freq2 in freqarr:

            if (freq2, freq1) in cl_dic:
                cl_dic[(freq1, freq2)] = cl_dic[(freq2, freq1)]
                continue

            #get dust
            el,  cl_dg_po, cl_dg_clus = fg.get_cl_dust(freq1, freq2, freq0 = param_dict['freq0'], fg_model = param_dict['fg_model'], spec_index_dg_po = param_dict['spec_index_dg_po'], spec_index_dg_clus = param_dict['spec_index_dg_clus'], Tcib = param_dict['Tcib'])
            if which_spec == 'EE':
                cl_dg_po = cl_dg_po * pol_frac_per_cent_dust**2.
                cl_dg_clus = cl_dg_clus * pol_frac_per_cent_dust**2.
            elif which_spec == 'TE':
                cl_dg_po = cl_dg_po * 0.
                cl_dg_clus = cl_dg_clus * 0.

            #include CIB decorrelation if available
            if cib_corr_coeffs is not None:
                if freq1 == freq2:
                    corr_coeff = 1.
                else:
                    if (freq1, freq2) in cib_corr_coeffs:
                        corr_coeff = cib_corr_coeffs[(freq1, freq2)]
                    elif (freq2, freq1) in cib_corr_coeffs:
                        corr_coeff = cib_corr_coeffs[(freq2, freq1)]
                cl_dg_po *= corr_coeff
                cl_dg_clus *= corr_coeff

            #get tsz
            el, cl_tsz = fg.get_cl_tsz(freq1, freq2, freq0 = param_dict['freq0'], fg_model = param_dict['fg_model'])
            if which_spec == 'EE':
                cl_tsz = cl_tsz * pol_frac_per_cent_tsz**2.
            elif which_spec == 'TE':
                cl_tsz = cl_tsz * 0.

            #get radio
            el, cl_radio = fg.get_cl_radio(freq1, freq2, freq0 = param_dict['freq0'], fg_model = param_dict['fg_model'], spec_index_rg = param_dict['spec_index_rg'])
            if which_spec == 'EE':
                cl_radio = cl_radio * pol_frac_per_cent_radio**2.
            elif which_spec == 'TE':
                cl_radio = cl_radio * 0.

            #get tSZ x CIB
            el, cl_tsz_cib = fg.get_cl_tsz_cib(freq1, freq2, freq0 = param_dict['freq0'], fg_model = param_dict['fg_model'], spec_index_dg_po = param_dict['spec_index_dg_po'], spec_index_dg_clus = param_dict['spec_index_dg_clus'], Tcib = param_dict['Tcib'])
            if which_spec == 'EE' or which_spec == 'TE':
                cl_tsz_cib = cl_tsz_cib * 0.
            if (0):
                freq_combs = [ (90, 90), (90,150), (90, 220), (150, 150), (150, 220), (220, 220)]
                colorarr = ['navy', 'blue', 'royalblue', 'green', 'lime', 'darkred']
                ax = subplot(111, yscale = 'log');
                for cntr, f1f2 in enumerate( freq_combs ):
                    freq1, freq2 = f1f2
                    el, cl_tsz_cib = fg.get_cl_tsz_cib(freq1, freq2, freq0 = param_dict['freq0'], fg_model = param_dict['fg_model'], spec_index_dg_po = param_dict['spec_index_dg_po'], spec_index_dg_clus = param_dict['spec_index_dg_clus'], Tcib = param_dict['Tcib'])
                    el, cl_tsz_cib_v2 = fg.get_foreground_power_spt('tSZ-CIB', freq1 = freq1, freq2 = freq2) 
                    plot(el, cl_tsz_cib, color = colorarr[cntr], lw = 1.); plot(el, -cl_tsz_cib_v2, color = colorarr[cntr], lw = 2., ls = '--'); 
                xlim(100, 5000)
                show()
                sys.exit()

            cl = np.copy( cl_ori )
            if 'cmb' not in ignore_fg:
                cl = cl + np.copy(cl_cmb)
            if 'ksz' not in ignore_fg:
                #print('ksz')
                cl = cl + cl_ksz             
            if 'tsz' not in ignore_fg:
                #print('tsz')
                cl = cl + cl_tsz
            if 'radio' not in ignore_fg:
                #print('radio')
                cl = cl + cl_radio
            if 'dust' not in ignore_fg:
                #print('dust')
                cl = cl + cl_dg_po + cl_dg_clus
            if 'dust' not in ignore_fg and 'tsz' not in ignore_fg:
                cl = cl + cl_tsz_cib


            if include_gal:# and not pol: #get galactic dust and sync

                el, cl_gal_dust = fg.get_cl_galactic(param_dict, 'dust', freq1, freq2, which_spec, bl_dic = bl_dic, el = el)
                el, cl_gal_sync = fg.get_cl_galactic(param_dict, 'sync', freq1, freq2, which_spec, bl_dic = bl_dic, el = el)

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

                el = np.arange(len(cl))

            else:
                nl = np.zeros(len(cl))
                #20191116 - fix this: there must be noise correlation in case of atmospheric noise
                print('\n\n\t\tfix me: there must be noise correlation in case of atmospheric noise')
                sys.exit()

            if which_spec != 'TE':
                cl = cl + nl                

            cl[np.isnan(cl)] = 0.
            cl[np.isinf(cl)] = 0.

            ##########################################################################################
            #20200516 - adjusting Nl when beam is too large (for 30/40 GHz bands)
            adjust_for_large_beams = 1                    
            if adjust_for_large_beams:
                bl = bl_dic[freq1]
                if 'effective' in bl_dic:
                    bl_eff = bl_dic['effective']
                else:
                    bl_eff = bl_dic[145]
                bad_inds = np.where( (bl_eff/bl > beam_tol_for_ilc) )[0]
                #print(nuval1, nuval2, bad_inds)
                cl[bad_inds] = 0.
            #20200516 - adjusting Nl when beam is too large (for 30/40 GHz bands)
            ##########################################################################################

            if (0):#which_spec == 'TT':
                loglog(cl); title('%s - %s' %(freq1, freq2)); show(); sys.exit()

            if (0):##which_spec == 'TE':
                #print(which_spec, freq1, freq2, cl)
                loglog(cl, label = r'%s - %s' %(freq1, freq2));
                legend(loc = 3, ncol = 2, fontsize = 6) 

            cl_dic[(freq1, freq2)] = cl

    return el, cl_dic  

################################################################################################################
################################################################################################################
################################################################################################################

def get_acap(freqarr, final_comp = 'cmb', freqcalib_fac = None):

    nc = len(freqarr)

    if freqcalib_fac is None: freqcalib_fac = np.ones(nc)

    if final_comp.lower() == 'cmb':
        freqscale_fac = np.ones(nc)

    elif final_comp.lower() == 'tsz' or final_comp.lower() == 'y':

        freqscale_fac = []
        for freq in sorted( freqarr ):
            freqscale_fac.append( compton_y_to_delta_Tcmb(freq * 1e9) )

        freqscale_fac = np.asarray( freqscale_fac )

        if final_comp.lower() == 'tsz': #tsz at 150 GHz
            freqarr = np.asarray( freqarr )
            if 150 in freqarr:
                freqind = np.where(freqarr == 150)[0]
            elif 145 in freqarr:
                freqind = np.where(freqarr == 145)[0]
            freqscale_fac = freqscale_fac / freqscale_fac[freqind]


    elif final_comp.lower() == 'cib':
        freqscale_fac = []
        for freq in sorted( freqarr ):
            freqscale_fac.append( get_cib_freq_dep(freq * 1e9) )

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

    clinv = sc.linalg.pinv2(clmat)
    #clinv = sc.linalg.inv(clmat)


    return clinv

def residual_power(param_dict, freqarr, el, cl_dic, which_spec, final_comp = 'cmb', freqcalib_fac = None, lmin = 0, return_weights = 0, null_comp = None):

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
        acap = get_acap(freqarr_for_acap, final_comp = final_comp, freqcalib_fac = freqcalib_fac)

        nc = len(freqarr)
        cl_residual = np.zeros( (len(el)) )
        weightsarr = np.zeros( (te_len * nc, len( el ) ) )
    
    if (1):
        acap = get_acap(freqarr, final_comp = final_comp, freqcalib_fac = freqcalib_fac)

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
        drinv = sc.linalg.pinv2(dr)
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

def get_acap_new(freqarr, final_comp = 'cmb', freqcalib_fac = None, teb_len = 1):

    nc = len(freqarr)

    if freqcalib_fac is None: freqcalib_fac = np.ones(nc)

    if final_comp.lower() == 'cmb':
        freqscale_fac = np.ones(nc)

    elif final_comp.lower() == 'tsz' or final_comp.lower() == 'y':

        freqscale_fac = []
        for freq in sorted( freqarr ):
            freqscale_fac.append( compton_y_to_delta_Tcmb(freq * 1e9) )

        freqscale_fac = np.asarray( freqscale_fac )

        if final_comp.lower() == 'tsz': #tsz at 150 GHz
            freqarr = np.asarray( freqarr )
            if 150 in freqarr:
                freqind = np.where(freqarr == 150)[0]
            elif 145 in freqarr:
                freqind = np.where(freqarr == 145)[0]
            freqscale_fac = freqscale_fac / freqscale_fac[freqind]

    elif final_comp.lower() == 'cib':
        freqscale_fac = []
        for freq in sorted( freqarr ):
            val = get_cib_freq_dep(freq * 1e9)
            freqscale_fac.append( val )
        freqscale_fac = np.asarray( freqscale_fac )
        freqscale_fac /= max(freqscale_fac)

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

def get_cib_freq_dep(nu, Tcib = 20., Tcmb = 2.7255, h=6.62607004e-34, k_B=1.38064852e-23, spec_index_dg = 1.505):
    if nu<1e3: nu *= 1e9

    bnu1 = fg.fn_BnuT(nu, temp = Tcib)
    dbdt = fg.fn_dB_dT(nu)
    value = (nu**spec_index_dg) * bnu1 / dbdt

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

def get_clinv_new(freqarr, elcnt, cl_dic, return_clmat = 0):
    clmat = np.mat( create_clmat_new(freqarr, elcnt, cl_dic) )
    #clmat = clmat + np.eye(len(clmat)) * np.min(np.diag(clmat)/1e3)
    clinv = sc.linalg.pinv2(clmat)
    #clinv = sc.linalg.inv(clmat)

    if return_clmat:
        return clinv, clmat
    else:
        return clinv

def residual_power_new(param_dict, freqarr, el, cl_dic, final_comp = 'cmb', freqcalib_fac = None, lmin = 10, return_weights = 0, null_comp = None):

    #acap = get_acap(freqarr, final_comp = final_comp, freqcalib_fac = freqcalib_fac)
    teb_len, pspec_arr = get_teb_spec_combination(cl_dic) #20200527 - teb
    acap = get_acap_new(freqarr, final_comp = final_comp, freqcalib_fac = freqcalib_fac, teb_len = teb_len)

    if null_comp is not None:
        total_comp_to_null = 0
        if np.ndim(null_comp) == 0:
            bcap = get_acap_new(freqarr, final_comp = null_comp, freqcalib_fac = freqcalib_fac, teb_len = teb_len)
            total_comp_to_null += 1
        else:
            #from IPython import embed; embed()
            bcap = None
            for curr_null_comp in null_comp:
                curr_bcap = get_acap_new(freqarr, final_comp = curr_null_comp, freqcalib_fac = freqcalib_fac, teb_len = teb_len)
                if bcap is None:
                    bcap = np.copy( curr_bcap )
                else:
                    bcap = np.column_stack( (bcap, curr_bcap) )
                total_comp_to_null += 1
            bcap = np.mat(bcap)
            #from IPython import embed; embed()
            #bcap = get_acap_new(freqarr, final_comp = null_comp[0], freqcalib_fac = freqcalib_fac, teb_len = teb_len)
            #print(bcap)
            #sys.exit()


    nc = len(freqarr)
    #cl_residual = np.zeros( (len(el)) )
    #weightsarr = np.zeros( (nc, len( el ) ) )
    weightsarr = np.zeros( (teb_len * nc, teb_len, len( el ) ) )
    cl_residual = np.zeros( (3, len(el)) )

    cl_residual_tmp = []

    for elcnt, el in enumerate(el):
        if el <= lmin: continue ## or el>=lmax: continue
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

            drinv = sc.linalg.pinv2(dr)
            weight = np.dot(nr, drinv)
        else:

            '''
            G = np.column_stack( (acap, bcap) )
            ncap = np.mat([1., 0.]).T
            '''
            G = np.column_stack( (acap, bcap) )
            G = np.mat(G)
            ncap = np.zeros( total_comp_to_null + 1 )
            ncap[0] = 1.
            ncap = np.mat( ncap ).T

            nr = np.dot(clinv, G)
            dr = np.dot( G.T, np.dot(clinv, G) )

            drinv = np.dot( sc.linalg.pinv2(dr), ncap)
            weight = np.dot(nr, drinv)

            if (0):
                acap_sum = np.sum( np.asarray(acap) * np.asarray(weight) )
                bcap_1_sum = np.sum( np.asarray(bcap[:, 0]) * np.asarray(weight) )
                if total_comp_to_null == 2:
                    bcap_2_sum = np.sum( np.asarray(bcap[:, 1]) * np.asarray(weight) )
                    print(weight, G, ncap, acap_sum, bcap_1_sum, bcap_2_sum)
                else:
                    print(weight, G, ncap, acap_sum, bcap_1_sum)
                sys.exit()

        #ILC residuals
        #cl_residual[elcnt] = np.asarray(1./dr).squeeze()
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

    from IPython import embed; embed()

    sys.exit()
    '''

    weightsarr = np.asarray( weightsarr )
    cl_residual = np.asarray( cl_residual )

    if final_comp.lower() == 'tsz': #tsz at 150 GHz
        freqarr = np.asarray( freqarr )
        if 150 in freqarr:
            freqind = np.where(freqarr == 150)[0]
        elif 145 in freqarr:
            freqind = np.where(freqarr == 145)[0]
        ysz_Tsz_conv_fac = compton_y_to_delta_Tcmb(freqarr[freqind] * 1e9)
        cl_residual = cl_residual / (ysz_Tsz_conv_fac**2.)

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

def compton_y_to_delta_Tcmb(nu, Tcmb = 2.73, h=6.62607004e-34, k_B=1.38064852e-23):

    if nu<1e3: nu *= 1e9

    x = h * nu / k_B / Tcmb
    g_nu = x * coth(x/2.) - 4.

    return Tcmb * np.mean(g_nu)

################################################################################################################

def get_ilc_map(final_comp, el, map_dic, bl_dic, nside, lmax, cl_dic = None, nl_dic = None, lmin = 10, freqcalib_fac = None, ignore_fg = [], full_sky = 1, estimate_covariance = 0, mapparams = None, apod_mask = None):

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
    weightsarr = get_multipole_weightsarr(final_comp, freqarr, el, cl_dic, lmin, freqcalib_fac)#, ignore_fg)
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

def get_multipole_weightsarr(final_comp, freqarr, el, cl_dic, lmin, freqcalib_fac):#, ignore_fg):

    #acap = get_acap(freqarr, final_comp = final_comp, freqcalib_fac = freqcalib_fac)
    teb_len, pspec_arr = get_teb_spec_combination(cl_dic) #20200527 - teb
    acap = get_acap(freqarr, final_comp = final_comp, freqcalib_fac = freqcalib_fac, teb_len = teb_len)

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
