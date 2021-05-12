from pylab import *
from matplotlib import rc;rc('text', usetex=True);rc('font', weight='bold');
import os
import argparse, sys, numpy as np, scipy as sc, warnings, os, glob

#import matplotlib.cbook
warnings.filterwarnings('ignore',category=RuntimeWarning)
#warnings.filterwarnings('ignore', category=DeprecationWarning)
#warnings.filterwarnings('ignore', category=matplotlib.cbook.mplDeprecation)

foldername = 's4like_mask_v2/TT-EE/baseline/' #change this accordingly
searchstr = '%s/*.npy' %(foldername)
flist = sorted( glob.glob(searchstr) )

xmin, xmax = 100, 5000
colorarr = [cm.jet(int(d)) for d in np.linspace(0, 255, len(flist))]
res_cl_dic = {}
for fcntr, fname in enumerate( flist ):
    if fname.find('years')>-1:
        yearval = float(fname.split('_')[-1].replace('for','').replace('years','').replace('.npy',''))
    else:
        yearval = 10.

    resdic = np.load(fname, allow_pickle = True).item()

    if fcntr == 0:
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
    res_cl_dic[yearval] = cl_residual

clf()
figure(figsize=(10., 5.))
fsval = 14
for cntr, which_spec in enumerate( which_spec_arr ):
    ax = subplot(1, len(which_spec_arr), cntr+1, yscale = 'log')#, xscale = 'log')
    for yearcntr, yearval in enumerate( sorted(res_cl_dic) ):
        labval = r'%g years' %(yearval)
        if yearval == 1.0:
            if which_spec == 'TT':
                plot(el_camb, cl_TT, color = 'gray')
            elif which_spec == 'EE':
                plot(el_camb, cl_EE, color = 'gray')
            labval = r'%g year' %(yearval)

        plot(el, res_cl_dic[yearval][which_spec], color = colorarr[yearcntr], lw = 1., label = labval)
    xlim(xmin, xmax);
    ylim(1e-7,1e-1);
    xlabel(r'Multipole $\ell$', fontsize = fsval)
    title(r'%s' %(which_spec), fontsize = fsval)
    if cntr == 0: 
        ylabel(r'$C_{\ell}$ [$\mu$K$^{2}$]', fontsize = fsval)
    if cntr == 1:
        legend(loc = 1, fontsize = fsval - 3, ncol = 2, handlelength = 2., handletextpad = 0.1)

    for label in ax.get_xticklabels(): label.set_fontsize(fsval-2)
    for label in ax.get_yticklabels(): label.set_fontsize(fsval-2)
plname = 'residual_ilc_curves.png'
savefig(plname, dpi = 200.)
show(); sys.exit()