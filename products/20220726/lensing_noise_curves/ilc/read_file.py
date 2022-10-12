import pickle, gzip, numpy as np, glob, sys, os
from pylab import *

plot_all_lensing_estimators = True


expname_arr = ['spt3g', 'sobaseline', 'sogoal', 's4wide', 's4deepv3r025']
expname_arr = ['spt3g', 's4wide']
expname_dic = {'sobaseline': 'SO-Baseline', 'sogoal': 'SO-Goal', 'spt3g': 'SPT-3G', 's4wide': 'S4-Wide', 's4deepv3r025': 'S4-Ultradeep'}
color_dic = {'sobaseline': 'royalblue', 'sogoal': 'darkorange', 'spt3g': 'goldenrod', 's4wide': 'darkgreen', 's4deepv3r025': 'darkred'}

reqd_lmaxtt = 5000

#read lensing noise curves
pl_dic = {}
expname_arr_final = []
for expcntr, expname in enumerate( expname_arr ):
    fname_arr = glob.glob('%s_lmin*_lmax*.npy' %(expname))
    pl_dic[expname] = {}
    for fname in fname_arr:
        fname_split = fname.split('_')
        lmax = fname_split[2]
        if fname.find('lmaxtt')>-1:
            lmaxtt = fname_split[3]
        else:
            lmaxtt = lmax
        lmax, lmaxtt = int(lmax.replace('lmax', '').replace('.npy', '')), int(lmaxtt.replace('lmax', '').replace('tt', '').replace('.npy', ''))
        if lmaxtt != reqd_lmaxtt: continue
        print(fname_split, lmax, lmaxtt)
        dic = np.load(fname, allow_pickle = 1, encoding='latin1').item()
        #print(dic.keys()); sys.exit()
        els, cl_kk, nl_tt, nl_eb, nl_mv, nl_mvpol = dic['els'], dic['cl_kk'], dic['Nl_TT'], dic['Nl_EB'], dic['Nl_MV'], dic['Nl_MVpol']
        nl_te, nl_tb, nl_ee = dic['Nl_ET'], dic['Nl_TB'], dic['Nl_EE']
        if fname.find('galmask') == -1:
            pl_dic[expname]['no_galaxy'] = els, cl_kk, nl_tt, nl_te, nl_tb, nl_ee, nl_eb, nl_mv, nl_mvpol
        else:
            pl_dic[expname]['with_galaxy'] = els, cl_kk, nl_tt, nl_te, nl_tb, nl_ee, nl_eb, nl_mv, nl_mvpol
        expname_arr_final.append( expname )

print('\nmake plots now\n')
if not plot_all_lensing_estimators:

    #make plot
    fsval = 14
    #xmin, xmax = min(els)+10, max(els)
    xmin, xmax = 0., 5000.
    ymin, ymax = 1e-9, 1e-5

    for which_est in ['MV', 'TT', 'MVpol']:
        clf()
        ax = subplot(111, yscale = 'log')
        for expcntr, expname in enumerate( expname_arr_final ):
            print(expname)
            els, cl_kk, nl_tt, nl_te, nl_tb, nl_ee, nl_eb, nl_mv, nl_mvpol = pl_dic[expname]['no_galaxy']
            labval = expname_dic[expname]
            if expname.find('s4')>-1:
                labval = r'%s ($\ell_{\rm max} = 5000; \ell_{\rm max}^{\rm TT} = %s$)' %(labval, reqd_lmaxtt)
            else:
                labval = r'%s ($\ell_{\rm max} = \ell_{\rm max}^{\rm TT} = %s$)' %(labval, reqd_lmaxtt)
            if which_est == 'MV':
                nl_to_plot = nl_mv
            elif which_est == 'TT':
                nl_to_plot = nl_tt
            elif which_est == 'MVpol':
                nl_to_plot = nl_mvpol

            ##print(els.shape, nl_to_plot.shape)
            nl_to_plot = nl_to_plot.real
            nl_to_plot[np.isnan(nl_to_plot)] = 0.

            plot(els, nl_to_plot, ls = '-', lw = 2., color = color_dic[expname], label = r'%s' %(labval))
            #plot(els, nl_tt, ls = '-.', lw = 1., color = colorarr[expcntr])
            #plot(els, nl_mvpol, ls = ':', lw = 1., color = colorarr[expcntr])

        #plot([], [], '-', label = r'MV', color = 'black')
        #plot([], [], '-.', label = r'TT', color = 'black')
        #plot([], [], ':', label = r'MV-Pol', color = 'black')
        plot(els, cl_kk, lw = 2., color = 'gray', alpha = 0.5)
        title(r'Lensing noise curves (%s)' %(which_est), fontsize = fsval)
        legend(loc = 2, fontsize = fsval-5)
        ylim(ymin, ymax); xlim(xmin, xmax)
        #setp(ax.get_xticklabels(which = 'both'), visible=False)
        ylabel(r'$L(L+1)]^{2}$C$_{L}^{\phi \phi}/2\pi$', fontsize = fsval)
        xlabel(r'Multipole $L$', fontsize = fsval)

        grid(True, which='major', axis = 'x', lw = 0.5, alpha = 0.3)
        grid(True, which='major', axis = 'y', lw = 0.5, alpha = 0.3)
        #grid(True, which='minor', axis = 'x', lw = 0.5, alpha = 0.2)
        grid(True, which='minor', axis = 'y', lw = 0.5, alpha = 0.2)

        plname = 'lensing_noise_curves_%s.png' %(which_est)
        savefig(plname, dpi = 200.)
        #show(); 
else:
    #make plot
    fsval = 14
    #xmin, xmax = min(els)+10, max(els)
    ymin, ymax = 1e-9, 1e-5
    color_dic_est = {'TT': 'darkgreen', 'TE': 'navy', 'TB': 'royalblue', 'EE': 'purple', 'EB': 'goldenrod', 'MVpol': 'orangered', 'MV': 'black'}
    for expcntr, expname in enumerate( expname_arr_final ):    

        print(expname)
        els, cl_kk, nl_tt, nl_te, nl_tb, nl_ee, nl_eb, nl_mv, nl_mvpol = pl_dic[expname]['no_galaxy']
        xmin, xmax = 0., max(els)
        print(els)

        titval = expname_dic[expname]
        if expname.find('s4')>-1 and xmax != reqd_lmaxtt:
            titval = r'%s ($\ell_{\rm max} = 5000; \ell_{\rm max}^{\rm TT} = %s$)' %(titval, reqd_lmaxtt)
        else:
            titval = r'%s ($\ell_{\rm max} = \ell_{\rm max}^{\rm TT} = %s$)' %(titval, reqd_lmaxtt)

        lensing_est_dict = {'TT': nl_tt, 'TE': nl_te, 'TB': nl_tb, 'EE': nl_ee, 'EB': nl_eb, 'MVpol': nl_mvpol, 'MV': nl_mv}
        clf()
        ax = subplot(111, yscale = 'log')
        for which_est in lensing_est_dict:
            ##print(els.shape, nl_to_plot.shape)
            nl_to_plot = lensing_est_dict[which_est]
            nl_to_plot = nl_to_plot.real
            nl_to_plot[np.isnan(nl_to_plot)] = 0.

            if which_est == 'MV' or which_est == 'MVpol':
                lwval = 2.
            else:
                lwval = 1.

            plot(els, nl_to_plot, ls = '-', lw = lwval, color = color_dic_est[which_est], label = r'%s' %(which_est))

        plot(els, cl_kk, lw = 2., color = 'gray', alpha = 0.5)
        title(r'Lensing noise curves: %s' %(titval), fontsize = fsval)
        legend(loc = 4, ncol = 4, fontsize = fsval-5)
        xmax=5000
        ylim(ymin, ymax); xlim(xmin, xmax)
        #setp(ax.get_xticklabels(which = 'both'), visible=False)
        ylabel(r'$L(L+1)]^{2}$C$_{L}^{\phi \phi}/2\pi$', fontsize = fsval)
        xlabel(r'Multipole $L$', fontsize = fsval)

        grid(True, which='major', axis = 'x', lw = 0.25, alpha = 0.3)
        grid(True, which='major', axis = 'y', lw = 0.25, alpha = 0.3)
        #grid(True, which='minor', axis = 'x', lw = 0.25, alpha = 0.2)
        grid(True, which='minor', axis = 'y', lw = 0.25, alpha = 0.2)

        plname = 'plots_all_estimators/lensing_noise_curves_%s_allestimators_lmax%s.png' %(expname, reqd_lmaxtt)
        savefig(plname, dpi = 200.)
        #show(); 
        close()
sys.exit()    