from pylab import *
from matplotlib import rc;rc('text', usetex=True);rc('font', weight='bold');matplotlib.rcParams['text.latex.preamble'] = [r'\boldmath']
import os
rcParams['font.family'] = 'serif'
rc('text.latex',preamble=r'\usepackage{/Users/sraghunathan/.configs/apjfonts}')


if (0): #make plots of S4 ILC curves
    s4_ilc_folder = 'results/20200601/s4like_mask//TT-EE-TE/baseline/'

    colorarr = ['black', 'darkred']
    include_gal_arr = [0, 1]
    for include_gal in include_gal_arr:
        if not include_gal:
            fname = '%s/S4_ilc_galaxy0_27-39-93-145-225-278_TT-EE-TE_lf2-mf12-hf5_AZ.npy' %(s4_ilc_folder)
        else:
            galmask = 3
            fname = '%s/S4_ilc_galaxy1_27-39-93-145-225-278_TT-EE-TE_lf2-mf12-hf5_galmask%s_AZ.npy' %(s4_ilc_folder, galmask)

        dic = np.load(fname, allow_pickle = 1).item()
        el_nl, cl_residual = dic['el'], dic['cl_residual']
        Nl_TT, Nl_EE = cl_residual['TT'], cl_residual['EE']

        clf()
        fsval = 14

        ax = subplot(111, yscale = 'log', xscale = 'log')
        plot(Nl_TT, ls = '-', color = 'black', label = r'ILC: TT')
        plot(Nl_EE, ls = '-', color = 'orangered', label = r'ILC: EE')

        if (1):
            #camb CMB
            camb_file = '%s/%s' %(dic['param_dict']['data_folder'], dic['param_dict']['Dlfile_len'])
            Tcmb= dic['param_dict']['T_cmb']
            el_camb = np.loadtxt(camb_file, usecols = [0])
            dl_camb = np.loadtxt(camb_file, usecols = [1,2,3,4])
            cl_camb = ( Tcmb**2. * dl_camb * 2 * np.pi ) / ( el_camb[:,None] * (el_camb[:,None] + 1) )
            cl_camb *= 1e12
            Dls_fac_camb = el_camb * (el_camb + 1) / 2/ np.pi
            plot(el_camb, cl_camb[:,0], color = 'gray', lw = 0.1)
            plot(el_camb, cl_camb[:,1], color = 'navy', lw = 0.1)
            #plot(el_camb, cl_camb[:,1], color = 'navy', lw = 0.5)

        xlim(10, 7000); #ylim(0.9, 1.1)
        xlabel(r'Multipole $\ell$', fontsize = fsval)
        galstr = ''
        if include_gal:
            galstr = ': Galaxy = %s: Galmask = %s' %(include_gal, galmask)
        title(r'S4 LAT ILC residuals (DSR specs)%s' %(galstr))
        legend(loc = 3, fontsize = fsval - 3, fancybox = 1)
        ylabel(r'$C_{\ell}\ [\mu K^{2}]$', fontsize = fsval)

        plname = 'reports/s4_lat_freq_dep/s4_LAT_ilc_curves_galaxy%s.png' %(include_gal)
        savefig(plname, dpi = 150)
    sys.exit()
    #show(); sys.exit()

if (0): #make raito plots of baseline vs other configurations
    clf()
    fig = figure(figsize=(8, 5))
    subplots_adjust(wspace = 0.05)
    fsval = 14
    s4_ilc_folder = 'results/20200601/s4like_mask/TT-EE-TE/'

    f1 = '%s/baseline//S4_ilc_galaxy1_27-39-93-145-225-278_TT-EE-TE_lf2-mf12-hf5_galmask1_AZ.npy' %(s4_ilc_folder)
    f2 = '%s/tubes_mod//S4_ilc_galaxy1_27-39-93-145-225-278_TT-EE-TE_lf1-mf12-hf6_galmask1_AZ.npy' %(s4_ilc_folder)
    f3 = '%s/tubes_mod//S4_ilc_galaxy1_27-39-93-145-225-278_TT-EE-TE_lf1-mf13-hf5_galmask1_AZ.npy' %(s4_ilc_folder)
    f4 = '%s/tubes_mod//S4_ilc_galaxy1_27-39-93-145-225-278_TT-EE-TE_lf1-mf14-hf4_galmask1_AZ.npy' %(s4_ilc_folder)
    f5 = '%s/tubes_mod//S4_ilc_galaxy1_27-39-93-145-225-278_TT-EE-TE_lf2-mf13-hf4_galmask1_AZ.npy' %(s4_ilc_folder)
    f6 = '%s/tubes_mod//S4_ilc_galaxy1_27-39-93-145-225-278_TT-EE-TE_lf3-mf12-hf4_galmask1_AZ.npy' %(s4_ilc_folder)

    farr= [f1, f2, f3, f4, f5, f6]

    colorarr = ['black', 'navy', 'darkgreen', 'goldenrod', 'orangered', 'darkred']
    for fcntr, fname in enumerate( farr ):
        dic = np.load(fname, allow_pickle = 1).item()
        el_nl, cl_residual = dic['el'], dic['cl_residual']
        Nl_TT, Nl_EE, Nl_TE = cl_residual['TT'], cl_residual['EE'], cl_residual['TE']

        if fcntr == 0:
            Nl_TT_baseline, Nl_EE_baseline, Nl_TE_baseline = Nl_TT, Nl_EE, Nl_TE

        labval = fname.split('/')[-1].split('_')[5].upper()
        if fcntr == 0:
            labval = '%s (Baseline)' %(labval)
        ax = subplot(121)
        plot(Nl_TT/Nl_TT_baseline, ls = '-', color = colorarr[fcntr], label = r'%s' %(labval))
        xlim(100, 7000);ylim(0.8, 1.3)
        xlabel(r'Multipole $\ell$', fontsize = fsval); ylabel(r'$N_{\ell}/N_{\ell}^{\rm baseline}$', fontsize = fsval)
        legend(loc = 2, fontsize = fsval - 3, fancybox = 1)
        title(r'TT')
        ax = subplot(122)
        plot(Nl_EE/Nl_EE_baseline, ls = '--', color = colorarr[fcntr])
        xlim(100, 7000);ylim(0.8, 1.3)
        xlabel(r'Multipole $\ell$', fontsize = fsval)
        #plot(Nldic['TT']/Nl_TT_baseline)                    
        title(r'EE')
        setp(ax.get_yticklabels(which = 'both'), visible=False)

    plname = 'reports/s4_lat_freq_dep/s4_LAT_baseline_otherconfigs_ilc_curves.png'
    savefig(plname, dpi = 200)
    sys.exit()
    #show(); sys.exit()


if (0): #make ratio plot of ILC for TTEE vs TTEETE
    clf()
    #fig = figure(figsize=(8, 5))
    #subplots_adjust(wspace = 0.05)
    fsval = 14
    s4_ilc_folder = 'results/20200610/s4like_mask/'

    colorarr = ['black', 'darkred']
    include_gal_arr = [0, 1]
    for include_gal in include_gal_arr:
        if not include_gal:
            f1 = '%s/S4_ilc_galaxy0_93-145-225-278_TT-EE-TE.npy' %(s4_ilc_folder)
            f2 = '%s/S4_ilc_galaxy0_93-145-225-278_TT-EE.npy' %(s4_ilc_folder)
        else:
            f1 = '%s/S4_ilc_galaxy1_93-145-225-278_TT-EE-TE_galmask3_AZ.npy' %(s4_ilc_folder)
            f2 = '%s/S4_ilc_galaxy1_93-145-225-278_TT-EE_galmask3_AZ.npy' %(s4_ilc_folder)

        farr= [f1, f2]

        for fcntr, fname in enumerate( farr ):
            dic = np.load(fname, allow_pickle = 1).item()
            el_nl, cl_residual = dic['el'], dic['cl_residual']
            Nl_TT, Nl_EE = cl_residual['TT'], cl_residual['EE']

            print(Nl_EE, fcntr)

            if fcntr == 0:
                Nl_TT_baseline, Nl_EE_baseline = Nl_TT, Nl_EE
                continue


            '''
            ax = subplot(121)
            plot(Nl_TT/Nl_TT_baseline, ls = '-', color = colorarr[fcntr])#, label = r'%s' %(labval))
            xlim(100, 7000);ylim(0.8, 1.3)
            xlabel(r'Multipole $\ell$', fontsize = fsval); ylabel(r'$N_{\ell}/N_{\ell}^{\rm baseline}$', fontsize = fsval)
            legend(loc = 2, fontsize = fsval - 3, fancybox = 1)
            title(r'TT')
            ax = subplot(122)
            plot(Nl_EE/Nl_EE_baseline, ls = '--', color = colorarr[fcntr])
            xlim(100, 7000);ylim(0.8, 1.3)
            xlabel(r'Multipole $\ell$', fontsize = fsval)
            #plot(Nldic['TT']/Nl_TT_baseline)                    
            title(r'EE')
            setp(ax.get_yticklabels(which = 'both'), visible=False)
            '''

            labval = r'Galaxy = %s' %(include_gal)
            ax = subplot(111)
            plot(Nl_EE/Nl_EE_baseline, ls = '--', color = colorarr[include_gal], label = r'%s' %(labval))
            xlim(100, 7000);ylim(0.9, 1.1)
            xlabel(r'Multipole $\ell$', fontsize = fsval)
            #plot(Nldic['TT']/Nl_TT_baseline)                    
            title(r'Ratio $N_{\ell}^{\rm EE}$: Separate (TT and EE), Joint (TT and EE)')
            legend(loc = 2, fontsize = fsval - 3, fancybox = 1)
            ylabel(r'$N_{\ell}^{\rm sep}/N_{\ell}^{\rm joint}$', fontsize = fsval)

    plname = 'reports/s4_lat_freq_dep/s4_LAT_sepTTEE_jointTTEE_baseline_4bands_ilc_curves.png'
    savefig(plname, dpi = 200)
    sys.exit()
    #show(); sys.exit()

if (0): #make ILC plot got Gil similar to Fig. 3 of https://arxiv.org/pdf/2006.06594.pdf

    def inv_var_noise(nl1, nl2):

        padzeros_len = len(nl1) - len(nl2)
        nl2_pad = np.zeros(padzeros_len) + 1e10
        nl2 = np.concatenate( (nl2, nl2_pad) )
        nl_ref = 1./ ( (1./nl1) + (1./nl2) )

        return nl_ref

    sys.path.append('/Users/sraghunathan/Research/SPTPol/analysis/git/DRAFT/modules/')
    import misc
    clf()
    fsval = 14
    s4_ilc_folder = 'results/20200610/s4like_mask/'
    

    colorarr = ['orangered', 'darkred']
    include_gal = 0

    if not include_gal:
        fname_wide = '%s/s4wide_ilc_galaxy0_27-39-93-145-225-278_TT-EE-TE.npy' %(s4_ilc_folder)
    else:
        fname_wide = '%s/S4_ilc_galaxy1_27-39-93-145-225-278_TT-EE-TE_galmask3_AZ.npy' %(s4_ilc_folder)

    fname_deep = '%s/s4deep_ilc_galaxy0_27-39-93-145-225-278_TT-EE-TE.npy' %(s4_ilc_folder)
    fname_arr = [fname_wide, fname_deep]

    #camb CMB
    dic = np.load(fname_wide, allow_pickle = 1).item()
    camb_file = '%s/%s' %(dic['param_dict']['data_folder'], dic['param_dict']['Dlfile_len'])
    Tcmb= dic['param_dict']['T_cmb']
    el_camb = np.loadtxt(camb_file, usecols = [0])
    dl_camb = np.loadtxt(camb_file, usecols = [1,2,3,4])
    cl_camb = ( Tcmb**2. * dl_camb * 2 * np.pi ) / ( el_camb[:,None] * (el_camb[:,None] + 1) )
    cl_camb *= 1e12
    Dls_fac_camb = el_camb * (el_camb + 1) / 2/ np.pi

    #kSZ
    ksz_file = '%s/dl_ksz_CSFplusPATCHY_13sep2011_norm1_fake12000.txt' %(dic['param_dict']['data_folder'])
    el_ksz, dl_ksz = np.loadtxt(ksz_file, unpack = 1)
    dl_ksz *= 3.

    #planck SMICA
    smica_file = '/Users/sraghunathan/Research/SPTPol/analysis/git/fisher_cosmo/data/Planck/DR2/Nl_smica_halfmission_HDHM.npy'
    smica_nl = np.load(smica_file, allow_pickle=1)

    ax = subplot(111, yscale = 'log')
    plot(el_camb,  cl_camb[:,0] * Dls_fac_camb, ls = '-', color = 'k', label = r'Lensed CMB')
    plot(el_ksz,  dl_ksz, ls = '-', color = 'goldenrod', label = r'kSZ')


    for fcntr, fname in enumerate( fname_arr ):
        if fcntr == 0:
            colorval = 'orangered'
            labval = 'ILC: S4-Wide'
        else:
            colorval = 'navy'
            labval = 'ILC: S4-Ultra deep'
        #ILC noise
        dic = np.load(fname, allow_pickle = 1).item()
        el_nl, cl_residual = dic['el'], dic['cl_residual']
        try:
            Nl_TT, Nl_EE = cl_residual['TT'], cl_residual['EE']
        except:
            Nl_TT, Nl_EE = cl_residual['T'], cl_residual['P']
        Dls_fac = el_nl * (el_nl + 1) / 2/ np.pi

        #inv variance noise
        Nl_TT_with_planck = inv_var_noise(Nl_TT, smica_nl)

        lwval = 1.
        #plot(el_nl, Nl_TT * Dls_fac, ls = '--', color = colorval, label = labval, lw = lwval)
        #plot(el_nl, Nl_TT_with_planck * Dls_fac, ls = '-', color = colorval, label = r'+{\it Planck}', lw = lwval)
        plot(el_nl, Nl_TT_with_planck * Dls_fac, ls = '-', color = colorval, label = r'%s + {\it Planck}' %(labval), lw = lwval)

        #inv variance noise
        #nl_145_with_planck = inv_var_noise(nl_145, smica_nl)
        Nl_TT_with_planck = inv_var_noise(Nl_TT, smica_nl)

        if (0):
            #150 GHz white noise level
            beamval, noiseval = dic['beam_noise_dic']['T'][145]
            elknee, alphaknee = elknee_dic = dic['elknee_dic']['T'][145]
            nl_145 = misc.get_nl(noiseval, el_nl, beamval, elknee = elknee, alphaknee = alphaknee)
            #nl_145 = misc.get_nl(noiseval, el_nl, beamval, elknee = -1, alphaknee = 0.)

            #inv variance noise
            #nl_145_with_planck = inv_var_noise(nl_145, smica_nl)
            Nl_TT_with_planck = inv_var_noise(Nl_TT, smica_nl)

            plot(el_nl, nl_145 * Dls_fac, ls = '--', color = 'navy', label = r'Noise: S4 145 GHz')
            plot(el_nl, nl_145_with_planck * Dls_fac, ls = '-', color = 'navy', label = r'+ {\it Planck}')


    xlim(100, 8000);ylim(1e-1, 1e4)
    xlabel(r'Multipole $\ell$', fontsize = fsval)
    legend(loc = 1, fontsize = fsval - 3, fancybox = 1)
    ylabel(r'$D_{\ell}$ [$\mu K^{2}$]', fontsize = fsval)

    plname = 'reports/s4_wide_deep_6bands_ilc_curves_for_Gil.png'
    savefig(plname, dpi = 200)
    #show()
    sys.exit()
    #show(); sys.exit()


if (1): #make ILC plot for SPT

    sys.path.append('/Users/sraghunathan/Research/SPTPol/analysis/git/DRAFT/modules/')
    import misc
    clf()
    subplots_adjust(hspace = 0.1)
    fsval = 14
    spt_ilc_folder = 'results/spt/'
    

    colorarr = ['orangered', 'darkred']
    include_gal = 0

    sptpolultradeep = '%s/sptpolultradeep_ilc_90-150-220_TT-EE.npy' %(spt_ilc_folder)
    #sptpolultradeepplus3g = '%s/sptpolultradeepplus3g_ilc_90-150-220_TT-EE.npy' %(spt_ilc_folder)
    sptpolplusultradeep = '%s/sptpolplusultradeep_ilc_90-150-220_TT-EE.npy' %(spt_ilc_folder)
    sptpolplusultradeepplus3g = '%s/sptpolplusultradeepplus3g_ilc_90-150-220_TT-EE.npy' %(spt_ilc_folder)
    sptpolplusultradeepplus3gfull = '%s/sptpolplusultradeepplus3gfull_ilc_90-150-220_TT-EE.npy' %(spt_ilc_folder)

    fname_arr = [sptpolultradeep, sptpolplusultradeep, sptpolplusultradeepplus3g, sptpolplusultradeepplus3gfull]
    lab_arr = [r'Ultradeep', r'Ultradeep+SPTpol', r'Ultradeep+SPTpol+3G', r'Ultradeep+SPTpol+3G-full']

    #camb CMB
    dic = np.load(sptpolultradeep, allow_pickle = 1).item()
    camb_file = '%s/%s' %(dic['param_dict']['data_folder'], dic['param_dict']['Dlfile_len'])
    Tcmb= dic['param_dict']['T_cmb']
    el_camb = np.loadtxt(camb_file, usecols = [0])
    dl_camb = np.loadtxt(camb_file, usecols = [1,2,3,4])
    cl_camb = ( Tcmb**2. * dl_camb * 2 * np.pi ) / ( el_camb[:,None] * (el_camb[:,None] + 1) )
    cl_camb *= 1e12
    Dls_fac_camb = el_camb * (el_camb + 1) / 2/ np.pi

    #kSZ
    ksz_file = '%s/dl_ksz_CSFplusPATCHY_13sep2011_norm1_fake12000.txt' %(dic['param_dict']['data_folder'])
    el_ksz, dl_ksz = np.loadtxt(ksz_file, unpack = 1)
    dl_ksz *= 3.

    tr, tc = 10, 1
    rspan = 6
    rspan2 = tr-rspan

    ax = subplot2grid((tr,tc), (0,0), rowspan = rspan, yscale = 'log')
    plot(el_camb,  cl_camb[:,0] * Dls_fac_camb, ls = '-', color = 'gray')#, label = r'Lensed CMB')
    plot(el_ksz,  dl_ksz, ls = '-', color = 'royalblue', label = r'kSZ')


    Nl_TT_arr, colorarr = [], []
    for fcntr, fname in enumerate( fname_arr ):
        if fcntr == 0:
            colorval = 'darkgreen'
        elif fcntr == 1:
            colorval = 'orangered'
        elif fcntr == 2:
            colorval = 'black'
        else:
            colorval = 'm'
        labval = lab_arr[fcntr]
        #ILC noise
        dic = np.load(fname, allow_pickle = 1).item()
        el_nl, cl_residual = dic['el'], dic['cl_residual']
        try:
            Nl_TT, Nl_EE = cl_residual['TT'], cl_residual['EE']
        except:
            Nl_TT, Nl_EE = cl_residual['T'], cl_residual['P']
        Dls_fac = el_nl * (el_nl + 1) / 2/ np.pi
        Nl_TT_arr.append(Nl_TT)
        colorarr.append(colorval)

        lwval = 1.
        plot(el_nl, Nl_TT * Dls_fac, ls = '-', color = colorval, label = labval, lw = lwval)
        #plot(el_nl, Nl_TT_with_planck * Dls_fac, ls = '-', color = colorval, label = r'+{\it Planck}', lw = lwval)
        #plot(el_nl, Nl_TT_with_planck * Dls_fac, ls = '-', color = colorval, label = r'%s + {\it Planck}' %(labval), lw = lwval)

    if (0):
        #Ultradeep only
        dic = np.load(sptpolultradeep, allow_pickle = 1).item()
        nl_fg_150 = dic['cl_dic']['TT'][(150, 150)]
        plot(el_nl, nl_fg_150 * Dls_fac, ls = '-.', color = colorval, label = labval, lw = lwval)

    xmin, xmax = 2000, 10000
    xlim(xmin, xmax);ylim(1., 1e4)
    #xlabel(r'Multipole $\ell$', fontsize = fsval)
    setp(ax.get_xticklabels(which = 'both'), visible=False)
    legend(loc = 2, fontsize = fsval - 4, fancybox = 1)#, ncol = 3)
    ylabel(r'$D_{\ell}$ [$\mu K^{2}$]', fontsize = fsval)
    title(r'SPTpol 100d ultradeep noise levels')

    ax = subplot2grid((tr,tc), (rspan,0), rowspan = rspan2)#, yscale = 'log')

    if (1): #Ultradeep only
        dic = np.load(sptpolultradeep, allow_pickle = 1).item()
        nl_fg_150 = dic['cl_dic']['TT'][(150, 150)]
        plot(el_nl, nl_fg_150 * Dls_fac, ls = '-.', color = colorval, label = labval, lw = lwval)


    Nl_TT_arr = np.asarray(Nl_TT_arr)
    Nl_TT_arr = Nl_TT_arr / nl_fg_150
    for ncntr in range(len(Nl_TT_arr)):
        plot(el_nl, Nl_TT_arr[ncntr], ls = '-', color = colorarr[ncntr], lw = lwval)
    xlim(xmin, xmax)#;ylim(1e-1, 1e3)
    ylim(0.5, 1.05)
    axhline(1.,lw=0.35)
    ylabel(r'ILC/C$_{\ell}^{150}$', fontsize = fsval)
    #plname = 'reports/s4_wide_deep_6bands_ilc_curves_for_Gil.png'
    #savefig(plname, dpi = 200)
    show()
    sys.exit()
    #show(); sys.exit()

if (0):
    clf(); 
    fsval = 8
    lwval = 0.75
    subplots_adjust(wspace=0.1, hspace = 0.1)
    xmin, xmax = 20, 7000
    which_spec_arr = ['TT', 'EE', 'TE']
    for cntr, which_spec in enumerate( which_spec_arr ):
        ax = subplot(1,3,cntr+1, xscale = 'log', yscale = 'log')
        plot(el, cl_residual[which_spec], 'black', lw = 2., label = r'Residual')
        xlim(xmin, xmax);
        ylim(1e-8,1e6);
        xlabel(r'Multipole $\ell$')
        if cntr == 0: 
            ylabel(r'$C_{\ell}$')
            ##legend(loc = 1, fontsize = 6, ncol = 2, handlelength = 2., handletextpad = 0.1)
        else:
            setp(ax.get_yticklabels(which = 'both'), visible=False)
        for label in ax.get_xticklabels(): label.set_fontsize(fsval)
        for label in ax.get_yticklabels(): label.set_fontsize(fsval)
        title(r'%s' %(which_spec))
    show()
