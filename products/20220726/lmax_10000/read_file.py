import numpy as np, pickle, gzip, glob
from pylab import *

flist = glob.glob('*.npy')
which_spec_arr = ['TT', 'EE']
expname_arr = ['spt3g', 'sobaseline', 'sogoal', 's4wide', 's4deepv3r025']
expname_dic = {'sobaseline': 'SO-Baseline', 'sogoal': 'SO-Goal', 'spt3g': 'SPT-3G', 's4wide': 'S4-Wide', 's4deepv3r025': 'S4-Ultradeep'}
color_dic = {'sobaseline': 'royalblue', 'sogoal': 'darkorange', 'spt3g': 'goldenrod', 's4wide': 'darkgreen', 's4deepv3r025': 'darkred'}

pl_dic = {}
for fname in flist:
    res_dic = np.load(fname, allow_pickle = 1, encoding = 'latin1').item()
    expname = fname.split('_')[0]
    if expname not in expname_arr: continue
    print(fname, expname, res_dic.keys())
    cl_residual_dic = res_dic['cl_residual']
    els = res_dic['el']
    for which_spec in which_spec_arr:
        if which_spec not in pl_dic:
            pl_dic[which_spec] = {}

        pl_dic[which_spec][expname] = cl_residual_dic[which_spec]


#make plot now
clf()
figure(figsize=(10., 5.2))
subplots_adjust(wspace=0.1)
fsval = 14.
xmin, xmax = 0, 10000. 
ymin, ymax = 1e-8, .1
for which_spec_cntr, which_spec in enumerate( which_spec_arr ):
    ax = subplot(1, 2, which_spec_cntr+1, yscale = 'log')
    for expname in expname_arr:
        curr_ilc_cl = pl_dic[which_spec][expname]
        plot(els, curr_ilc_cl, color = color_dic[expname], label = r'%s' %(expname_dic[expname]), lw = 2.)

    xlim(xmin, xmax);
    ylim(ymin, ymax);
    xlabel(r'Multipole $\ell$', fontsize = fsval)
    if which_spec_cntr == 0:
        ylabel(r'$C_{\ell}$ [$\mu$K$^{2}$]', fontsize = fsval)
        legend(loc = 1, fontsize = fsval, ncol = 1, handlelength = 2., handletextpad = 0.1)
    else:
        #pass
        setp(ax.get_yticklabels(which = 'both'), visible=False)
    for label in ax.get_xticklabels(): label.set_fontsize(fsval)
    for label in ax.get_yticklabels(): label.set_fontsize(fsval)

    title(r'%s' %(which_spec), fontsize = fsval)
    grid(True, which='major', axis = 'x', lw = 0.5, alpha = 0.3)
    grid(True, which='major', axis = 'y', lw = 0.5, alpha = 0.3)
    #grid(True, which='minor', axis = 'x', lw = 0.5, alpha = 0.2)
    grid(True, which='minor', axis = 'y', lw = 0.5, alpha = 0.2)
suptitle(r'Standard ILC residuals', fontsize = fsval + 4)
plname = 'ilc_residuals.png'
savefig(plname, dpi = 200.)
#show()