from pylab import *
from matplotlib import rc;rc('text', usetex=True);rc('font', weight='bold'); #matplotlib.rcParams['text.latex.preamble'] = r'\boldmath'
import os
#rc('text.latex',preamble=r'\usepackage{/Users/sraghunathan/.configs/apjfonts}')

import argparse, sys, numpy as np, scipy as sc, warnings, os, glob

sys.path.append('/Users/sraghunathan/Research/SPTPol/analysis/git/DRAFT/modules/')
import foregrounds as fg, misc, flatsky, misc, exp_specs

#import matplotlib.cbook
warnings.filterwarnings('ignore',category=RuntimeWarning)
#warnings.filterwarnings('ignore', category=DeprecationWarning)
#warnings.filterwarnings('ignore', category=matplotlib.cbook.mplDeprecation)


expdic = {\
's4wideclean': ['s4like_mask_v2/TT-EE/baseline/', 'galmask2', 'S4-Wide', 'S4-Clean', 'black', '-'], \
's4widedirty': ['s4like_mask_v2/TT-EE/baseline/', 'galmask5', 'S4-Wide', 'S4-Galactic Plane', 'black', '--'], \
's4deep': ['s4delensing_mask/TT-EE/baseline/', '', 'S4-Ultra deep', 'S4-Ultra deep', 'black', '-'], \
}
reqd_yearval = 7.
res_cl_dic = {}
beam_noise_dic, elknee_dic = {}, {}
for expcntr, expname in enumerate( expdic ):
    #fd = expdic[expname][0]
    #searchstr = '%s/*%dyears*.npy' %(fd, reqd_yearval)
    fd, searchstr_tmp, tit, explabel, colorval, lsval = expdic[expname]
    searchstr = '%s/*%s*%dyears*.npy' %(fd, searchstr_tmp, reqd_yearval)
    fname = sorted( glob.glob(searchstr) )[0]

    print(fname)
    resdic = np.load(fname, allow_pickle = True).item()

    if expcntr == 0:
        param_dict = resdic['param_dict']
        camb_file = '%s/%s' %(param_dict['data_folder'], param_dict['Dlfile_len'])
        el_camb = np.loadtxt(camb_file, usecols = [0])
        dl_camb = np.loadtxt(camb_file, usecols = [1,2,3,4])
        Tcmb = 2.73
        cl_camb = ( Tcmb**2. * dl_camb * 2 * np.pi ) / ( el_camb[:,None] * (el_camb[:,None] + 1) )
        cl_camb *= 1e12
        cl_TT, cl_EE, cl_BB, cl_TE = cl_camb.T

    el = resdic['el']
    cl_residual = resdic['cl_residual']
    which_spec_arr = cl_residual.keys()
    res_cl_dic[expname] = cl_residual
    beam_noise_dic[expname] = resdic['beam_noise_dic']
    elknee_dic[expname] = resdic['elknee_dic']

clf()
figure(figsize=(10., 6.3))
subplots_adjust(wspace = 0.02, hspace = 0.05)
fsval = 13
sbpl = 1
delta_el = 100
el_mod = np.arange(min(el), max(el), delta_el)
dl_or_cl = 'cl'
if dl_or_cl == 'dl':
    dl_fac = el * (el+1)/2/np.pi
    dl_fac_camb = el_camb * (el_camb+1)/2/np.pi
    dl_fac_mod = el_mod * (el_mod+1)/2/np.pi
    ymin, ymax = 0.01, 1e4
    ylabval = 'D_{\ell}'
else:
    dl_fac = 1.
    dl_fac_camb = 1.
    dl_fac_mod = 1.
    ymin, ymax = 1e-8,1e-1
    ylabval = 'C_{\ell}'
xmin, xmax = 250, 5000
for cntr, which_spec in enumerate( which_spec_arr ):
    for expcntr, expname in enumerate( expdic ):
        if expname.find('s4wide')>-1:
            if which_spec == 'TT':
                sbpl = 1
            elif which_spec == 'EE':
                sbpl = 3
        elif expname.find('s4deep')>-1:
            if which_spec == 'TT':
                sbpl = 2
            elif which_spec == 'EE':
                sbpl = 4
        ax = subplot(2, len(which_spec_arr), sbpl, yscale = 'log')#, xscale = 'log')
        #explabel, colorval = r'%s' %(expdic[expname][1]), expdic[expname][2]
        fd, searchstr_tmp, tit, explabel, colorval, lsval = expdic[expname]
        labval = r'%s' %(which_spec)
        if which_spec == 'TT':
            plot(el_camb, cl_TT * dl_fac_camb, color = 'gray', alpha = 0.5)
            TPval = 'T'
        elif which_spec == 'EE':
            plot(el_camb, cl_EE * dl_fac_camb, color = 'gray', alpha = 0.5)
            TPval = 'P'

        if (1): #plot individual band noise
            nuarr = sorted(beam_noise_dic[expname][TPval].keys())
            #cmap = cm.gist_rainbow #Spectral #cm.rainbow #cm.jet
            cmap = cm.jet
            colorarr = [cmap(int(d)) for d in np.linspace(0, 255, len(nuarr))]
            #colorarr = ['navy', 'blue', 'green', 'goldenrod', 'red', 'darkred']
            for nucntr, nu in enumerate( nuarr ):
                beamval, noiseval = beam_noise_dic[expname][TPval][nu]
                elknee, alphaknee = elknee_dic[expname][TPval][nu]
                nl = misc.get_nl(noiseval, el, beamval, elknee = elknee, alphaknee = alphaknee)
                #plot(el, nl * dl_fac, label = r'%s GHz' %(nu), color = colorarr[nucntr], ls = '-', lw = 0.5)
                if sbpl == 3 and expcntr == 0:
                    labval = r'$N_{\ell}^{%s}$' %(nu)
                else:
                    labval = None
                plot(el, nl * dl_fac, label = labval, color = colorarr[nucntr], ls = '-', lw = 0.5)

        cl_residual = res_cl_dic[expname][which_spec]
        if (1):
            cl_residual_mod = np.interp(el_mod, el, cl_residual)        

        if sbpl == 1:
            labval = r'%s' %(explabel)
        else:
            labval = None
        plot(el_mod, cl_residual_mod * dl_fac_mod, color = colorval, ls = lsval, lw = 2., label = labval)
        xlim(xmin, xmax); 
        ylim(ymin, ymax);

        if cntr == 0:
            title(r'%s' %(tit), fontsize = fsval+2)
            setp(ax.get_xticklabels(which = 'both'), visible=False)
        else:
            xlabel(r'Multipole $\ell$', fontsize = fsval)

        if (sbpl-1)%2 == 0: 
            ylabel(r'$%s^{%s}$ [$\mu$K$^{2}$]' %(ylabval, which_spec), fontsize = fsval)
        else:
            setp(ax.get_yticklabels(which = 'both'), visible=False)

        #if cntr == 0 and expcntr == 0:
        if sbpl==1 or sbpl==3:
            legend(loc = 3, fontsize = fsval - 2, ncol = 3, handletextpad = 0.15, columnspacing = 0.5, framealpha = 0., handlelength = 2.)


        for label in ax.get_xticklabels(): label.set_fontsize(fsval-2)
        for label in ax.get_yticklabels(): label.set_fontsize(fsval-2)
        grid(True, which='major', axis = 'x', lw = 0.25, alpha = 0.2)
        grid(True, which='both', axis = 'y', lw = 0.25, alpha = 0.2)
        ax.xaxis.set_major_locator(MaxNLocator(nbins=5))
        ax.xaxis.set_minor_locator(MaxNLocator(nbins=20))
        ax.tick_params(axis = 'x', direction='in', length=4., width=1, which = 'major')
        ax.tick_params(axis = 'x', direction='in', length=2., width=1, which = 'minor')

        #sbpl += 1
        #show(); sys.exit()


plname = 'residual_ilc_curves_s4wide_ultradeep_202102designtoolinputgalsims_SPTextragal.pdf'
savefig(plname, dpi = 200.)
show(); sys.exit()