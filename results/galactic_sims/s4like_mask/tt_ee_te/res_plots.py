import numpy as np, glob, sys, re

from pylab import *
from matplotlib import rc;rc('text', usetex=True);rc('font', weight='bold');matplotlib.rcParams['text.latex.preamble'] = [r'\boldmath']
import os
rc('text.latex',preamble=r'\usepackage{/Users/sraghunathan/.configs/apjfonts}')
rcParams['font.family'] = 'serif'
rcParams['figure.dpi'] = 150
rcParams["figure.facecolor"] = 'white'

basefolder = ''
basefolder = 'with_gal/'
baselinesearchstr = 'lf2-mf12-hf5'

if (0):
    searchstr = 'hf2'
    opfd = 'plots/consolidated/'#hf2/'

if (0):
    searchstr = 'hf3'
    opfd = 'plots/consolidated/'#hf3/'

if (0):
    searchstr = 'hf1'
    opfd = 'plots/consolidated/'#hf3/'

if (0):
    searchstr = 'hf4'
    opfd = 'plots/consolidated/'#hf3/'

if (1):
    searchstr = 'lf1'
    opfd = 'plots/consolidated/'#lf1/'

if (0):
    searchstr = 'lf0'
    opfd = 'plots/consolidated/'#lf0/'

opfd = '%s%s' %(basefolder, opfd)
coloarr = ['black', 'navy', 'darkgreen', 'goldenrod', 'darkred']
also_te = 1
if not also_te:
    basefname = glob.glob('%s*TT-EE_*%s*' %(basefolder,baselinesearchstr))
    fnames = sorted( glob.glob('%s*TT-EE_*%s*' %(basefolder,searchstr)) )
else:
    basefname = glob.glob('%s*TT-EE-TE_*%s*' %(basefolder,baselinesearchstr))
    fnames = sorted( glob.glob('%s*TT-EE-TE_*%s*' %(basefolder,searchstr)) )

flist = np.concatenate( (basefname, fnames))
fsval = 8
ratio_plots = 1
clf()
fig = figure(figsize=(5.,2.8))
subplots_adjust(wspace = 0.01)
for fcntr, f in enumerate( flist ):
    #if f.find('TT-EE-TE')>-1: continue
    print(f)
    dic = np.load(f, allow_pickle = 1).item()['cl_residual']
    if also_te:
        cl_tt, cl_ee, cl_te = dic['TT'],dic['EE'], dic['TE']
        cl_arr = [cl_tt, cl_ee, cl_te]
    else:
        cl_tt, cl_ee = dic['TT'],dic['EE']
        cl_arr = [cl_tt, cl_ee]
    colorval = coloarr[fcntr]

    labval = re.findall(r'lf[0-9]*-mf[0-9]*-hf[0-9]*', f)[0]
    if f.find(baselinesearchstr)>-1:
        colorval = 'black'
        lsval = 'solid'
        labval = '%s (Baseline)' %(labval)
        cl_arr_baseline = np.copy(cl_arr)
    else:
        lsval = 'dashed'
    lwval = 1.
    for clcntr, (cl, base_cl) in enumerate( zip(cl_arr, cl_arr_baseline) ):
        if ratio_plots:
            ax = subplot(1, len(cl_arr), clcntr+1)
            plot(cl/base_cl, color = colorval, ls = lsval, lw = lwval, label = r'%s' %(labval))
        else:
            ax = subplot(1, len(cl_arr), clcntr+1, xscale = 'log', yscale = 'log')
            plot(cl, color = colorval, ls = lsval, lw = lwval, label = r'%s' %(labval))

        if clcntr == 0:
            tit = 'TT'
        elif clcntr == 1:
            tit = 'EE'
        elif clcntr == 2:
            tit = 'TE'
        title(r'%s' %(tit), fontsize = fsval)

        xlim(20, 6500);
        #ylim(1e-8,1e6);
        if not ratio_plots:
            ylim(1e-7,1.);
        else:
            ylim(0.5, 1.5);
        xlabel(r'Multipole $\ell$')
        if clcntr == 0: 
            if not ratio_plots:
                ylabel(r'$C_{\ell}\ [\mu K^{2}]$')
            else:
                ylabel(r'$C_{\ell}/C_{\ell}^{\rm baseline}$')
        else:
            if also_te:
                if clcntr == 2:
                    legend(loc = 1, fontsize = 8, ncol = 1, handlelength = 2., handletextpad = 0.1)
            else:
                legend(loc = 1, fontsize = 8, ncol = 1, handlelength = 2., handletextpad = 0.1)
            setp(ax.get_yticklabels(which = 'both'), visible=False)
        for label in ax.get_xticklabels(): label.set_fontsize(fsval)
        for label in ax.get_yticklabels(): label.set_fontsize(fsval)
if also_te:
    plname = '%s/tteete_combined_%s.png' %(opfd, searchstr)
else:
    plname = '%s/ttee_combined_%s.png' %(opfd, searchstr)
if ratio_plots:
    plname = plname.replace('combined','ratio')
savefig(plname)
#show()
sys.exit()