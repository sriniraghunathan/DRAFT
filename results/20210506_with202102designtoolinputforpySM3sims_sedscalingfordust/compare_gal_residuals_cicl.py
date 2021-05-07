from pylab import *
from matplotlib import rc;rc('text', usetex=True);rc('font', weight='bold'); #matplotlib.rcParams['text.latex.preamble'] = r'\boldmath'
import os
#rc('text.latex',preamble=r'\usepackage{/Users/sraghunathan/.configs/apjfonts}')

import argparse, sys, numpy as np, scipy as sc, warnings, os, glob, re

sys.path.append('/Users/sraghunathan/Research/SPTPol/analysis/git/DRAFT/modules/')
import foregrounds as fg, misc, flatsky, misc, exp_specs

#import matplotlib.cbook
warnings.filterwarnings('ignore',category=RuntimeWarning)
#warnings.filterwarnings('ignore', category=DeprecationWarning)
#warnings.filterwarnings('ignore', category=matplotlib.cbook.mplDeprecation)

cilc_fd_arr = glob.glob('nulled*/s4like_mask_v2/TT-EE/baseline/')
ilc_fd_arr = glob.glob('s4like_mask_v2/TT-EE/baseline/')
fd_arr = np.concatenate( (ilc_fd_arr, cilc_fd_arr))

reqd_yearval = 7.

coloar_arr =['black', 'royalblue', 'darkgreen', 'goldenrod', 'red']

pl_dic = {}
for fdcntr, fd in enumerate( fd_arr ):
    searchstr = '%s/*%dyears*.npy' %(fd, reqd_yearval)
    fname = sorted( glob.glob(searchstr) )[0]
    resdic = np.load(fname, allow_pickle = True).item()

    if fdcntr == 0:
        param_dict = resdic['param_dict']
        camb_file = '%s/%s' %(param_dict['data_folder'], param_dict['Dlfile_len'])
        el_camb = np.loadtxt(camb_file, usecols = [0])
        dl_camb = np.loadtxt(camb_file, usecols = [1,2,3,4])
        Tcmb = 2.73
        cl_camb = ( Tcmb**2. * dl_camb * 2 * np.pi ) / ( el_camb[:,None] * (el_camb[:,None] + 1) )
        cl_camb *= 1e12
        cl_TT, cl_EE, cl_BB, cl_TE = cl_camb.T

    el = resdic['el']
    cl_residual_dic = resdic['cl_residual']
    fg_res_dic = resdic['fg_res_dic']
    which_spec_arr = cl_residual_dic.keys()

    if fd.find('nulled') == -1:
        fdlab = 'MV ILC'
    else:
        fdlab = fd.split('/')[0].replace('nulled_misc_cib_', '')
        tcib_tmp = float( re.findall('tcib\d*\.?\d+', fdlab.lower())[0].replace('tcib', '') )
        beta_tmp = float( re.findall('beta\d*\.?\d+', fdlab.lower())[0].replace('beta', '') )
        fdlab = r'cILC: T$_{\rm d}$ = %g, $\beta_{\rm d}$=%g' %(tcib_tmp, beta_tmp)

    pl_dic[fd] = {}
    for which_spec in which_spec_arr:
        pl_dic[fd][which_spec] = [el, cl_residual_dic[which_spec], fg_res_dic[which_spec]['galdust'], fdlab, coloar_arr[fdcntr]]


clf()
figure(figsize=(8., 4.5))
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
    ymin, ymax = 1e-8,1e-2
    ylabval = 'C_{\ell}'
xmin, xmax = 250, 5000
for whichspeccntr, which_spec in enumerate( which_spec_arr ):
    ax = subplot(1, len(which_spec_arr), whichspeccntr+1, yscale = 'log')#, xscale = 'log')        
    for fdcntr, fd in enumerate( fd_arr ):
        labval = r'%s' %(which_spec)
        if fdcntr==0:
            if which_spec == 'TT':
                plot(el_camb, cl_TT * dl_fac_camb, color = 'gray', alpha = 0.5)
                TPval = 'T'
            elif which_spec == 'EE':
                plot(el_camb, cl_EE * dl_fac_camb, color = 'gray', alpha = 0.5)
                TPval = 'P'
            plot([], [], 'k-', label = r'Total residual')
            plot([], [], 'k-.', label = r'Gal. dust')

        el, cl_residual, gal_dust_res, fdlab, colorval = pl_dic[fd][which_spec]
        if (1):
            cl_residual_mod = np.interp(el_mod, el, cl_residual)
            gal_dust_res_mod = np.interp(el_mod, el, gal_dust_res)

        plot(el_mod, cl_residual_mod * dl_fac_mod, color = colorval, ls = '-', lw = 1., label = fdlab)
        plot(el_mod, gal_dust_res_mod * dl_fac_mod, color = colorval, ls = '-.', lw = 1.)
        xlim(xmin, xmax); 
        ylim(ymin, ymax);

        title(r'%s' %(which_spec), fontsize = fsval+2)
        xlabel(r'Multipole $\ell$', fontsize = fsval)

        if whichspeccntr == 0: 
            ylabel(r'$%s$ [$\mu$K$^{2}$]' %(ylabval), fontsize = fsval)
        else:
            setp(ax.get_yticklabels(which = 'both'), visible=False)

        #if cntr == 0 and expcntr == 0:
        if whichspeccntr==1:
            legend(loc = 1, fontsize = fsval - 3, ncol = 1, handletextpad = 0.2, columnspacing = 0.5, framealpha = 0., handlelength = 2.)

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


plname = 'residual_ilc_cilc_curves_plus_galdustresiduals.pdf'
#savefig(plname, dpi = 200.)
show(); sys.exit()