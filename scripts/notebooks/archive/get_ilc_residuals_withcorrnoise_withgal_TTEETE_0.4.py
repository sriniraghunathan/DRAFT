#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '2')

#%pylab notebook
get_ipython().run_line_magic('matplotlib', 'inline')
from pylab import *
from matplotlib import rc;rc('text', usetex=True);rc('font', weight='bold');matplotlib.rcParams['text.latex.preamble'] = [r'\boldmath']
import os
rc('text.latex',preamble=r'\usepackage{/Users/sraghunathan/.configs/apjfonts}')


# In[2]:


rcParams['figure.dpi'] = 150
rcParams['font.family'] = 'serif'
rcParams["figure.facecolor"] = 'white'


# In[3]:


import argparse, sys, numpy as np, scipy as sc, warnings, os
sys.path.append('/Users/sraghunathan/Research/SPTPol/analysis/git/DRAFT/modules/')
import flatsky, misc, exp_specs
import ilc, foregrounds as fg

#import matplotlib.cbook
warnings.filterwarnings('ignore',category=RuntimeWarning)
#warnings.filterwarnings('ignore', category=DeprecationWarning)
#warnings.filterwarnings('ignore', category=matplotlib.cbook.mplDeprecation)


# In[4]:


#some constants
h=6.62607004e-34 #Planck constant in m2 kg / s
k_B=1.38064852e-23 #Boltzmann constant in m2 kg s-2 / K-1
Tcmb = 2.73 #Kelvin


# In[5]:


#params
paramfile = 'params.ini'

# read and store param dict
param_dict = misc.fn_get_param_dict(paramfile)
el = np.arange(param_dict['lmax'])
include_gal = param_dict['include_gal'] ##1
s4like_mask = param_dict['s4like_mask']
s4like_mask_v2 = param_dict['s4like_mask_v2']
if not include_gal:
    if s4like_mask:
        param_dict['which_gal_mask'] = 3
    elif s4like_mask_v2:
        param_dict['which_gal_mask'] = 2
which_gal_mask = param_dict['which_gal_mask']
try:
    remove_atm = param_dict['remove_atm']
except:
    remove_atm = 0    


# In[6]:


#S4 specs
expname = 's4wide'
#expname = 's4deep'
specs_dic, corr_noise_bands, rho, corr_noise = exp_specs.get_exp_specs(expname, remove_atm = remove_atm)
freqarr = sorted( specs_dic.keys() )
nc = len( freqarr )
freqcalib_fac = None
final_comp = 'cmb'
TParr = ['T', 'P']
#which_spec_arr = ['TT', 'EE', 'TE']
which_spec_arr = ['TT', 'EE']
#which_spec_arr = ['TE']
##include_gal = 1


# In[7]:


#beam and noise arr
beamarr = []
noisearr_T, elkneearr_T, alphakneearr_T = [], [], []
noisearr_P, elkneearr_P, alphakneearr_P = [], [], []
for freq in freqarr:
    beam_arcmins, white_noise_T, elknee_T, alphaknee_T, whitenoise_P, elknee_P, alphaknee_P = specs_dic[freq]
    beamarr.append(beam_arcmins)
    noisearr_T.append(white_noise_T)
    noisearr_P.append(whitenoise_P)
    elkneearr_T.append(elknee_T)
    elkneearr_P.append(elknee_P)
    alphakneearr_T.append(alphaknee_T)
    alphakneearr_P.append(alphaknee_P)    

print(elkneearr_T)


# In[8]:


#collect beam and noise into a dic; elknee and alpha into a dic
beam_noise_dic = {}
elknee_dic = {}
for TP in TParr:
    beam_noise_dic[TP] = {}
    elknee_dic[TP] = {} 
    if TP == 'T':
        freqarr, beamarr, noisearr, elkneearr, alphakneearr = freqarr, beamarr, noisearr_T, elkneearr_T, alphakneearr_T
    elif TP == 'P':
        freqarr, beamarr, noisearr, elkneearr, alphakneearr = freqarr, beamarr, noisearr_P, elkneearr_P, alphakneearr_P

    for (freq, beam, noise, elknee, alphaknee) in zip(freqarr, beamarr, noisearr, elkneearr, alphakneearr):
        beam_noise_dic[TP][freq] = [beam, noise]
        elknee_dic[TP][freq] = [elknee, alphaknee]


# In[9]:


#get beam deconvolved noise nls
nl_dic = {}
for TP in TParr:
    nl_dic[TP]={}
    for freq1 in freqarr:
        beamval1, noiseval1 = beam_noise_dic[TP][freq1]
        elknee1, alphaknee1 = elknee_dic[TP][freq1]
        for freq2 in freqarr:        
            beamval2, noiseval2 = beam_noise_dic[TP][freq2]
            elknee2, alphaknee2 = elknee_dic[TP][freq2]
            
            if freq1 == freq2:
                nl = misc.get_nl(noiseval1, el, beamval1, elknee = elknee1, alphaknee = alphaknee1)
            else:
                if freq2 in corr_noise_bands[freq1]:
                    nl = misc.get_nl(noiseval1, el, beamval1, elknee = elknee1, alphaknee = alphaknee1,                                      beamval2 = beamval2, noiseval2 = noiseval2, elknee2 = elknee2, alphaknee2 = alphaknee2, rho = rho)
                else:
                    nl = np.zeros( len(el) )
            nl[el<=param_dict['lmin']] = 0.
            ##nl[nl == 0.] = np.min(nl[nl!=0.])/1e3
            nl_dic[TP][(freq1, freq2)] = nl
print(nl_dic['T'].keys())


# In[10]:


#get beams
bl_dic = misc.get_beam_dic(freqarr, beam_noise_dic['T'], param_dict['lmax'])
print(bl_dic.keys())
if (0):
    for freq in freqarr:
        plot(bl_dic[freq], label = freq)
    legend(loc = 1)


# In[11]:


#get the CMB, noise, and foreground covriance
try:
    ignore_fg = param_dict['ignore_fg']
except:
    ignore_fg = []

ignore_fg.append(final_comp.lower()) #the required component need not go into the covariance matrix.
print(ignore_fg)

#freqarr = [145]
#param_dict['which_gal_mask'] = 0
cl_dic = {}
for which_spec in which_spec_arr:
    if which_spec == 'TT':
        el, cl_dic[which_spec] = ilc.get_analytic_covariance(param_dict, freqarr,                 nl_dic = nl_dic['T'], ignore_fg = ignore_fg, include_gal = include_gal, bl_dic = bl_dic)
    else:
        el, cl_dic[which_spec] = ilc.get_analytic_covariance                    (param_dict, freqarr, nl_dic = nl_dic['P'], ignore_fg = ignore_fg, which_spec = which_spec,                     pol_frac_per_cent_dust = param_dict['pol_frac_per_cent_dust'],                     pol_frac_per_cent_radio = param_dict['pol_frac_per_cent_radio'],                     pol_frac_per_cent_tsz = param_dict['pol_frac_per_cent_tsz'],                     pol_frac_per_cent_ksz = param_dict['pol_frac_per_cent_ksz'],                     include_gal = include_gal, bl_dic = bl_dic)
print(cl_dic.keys(), cl_dic.keys())


# In[21]:


if (0):
    el,  cl_dg_po, cl_dg_clus = fg.get_cl_dust(145, 145, freq0 = param_dict['freq0'], fg_model = param_dict['fg_model'], spec_index_dg_po = param_dict['spec_index_dg_po'], spec_index_dg_clus = param_dict['spec_index_dg_clus'], Tcib = param_dict['Tcib'])
    el, cl_radio = fg.get_cl_radio(145, 145, freq0 = param_dict['freq0'], fg_model = param_dict['fg_model'], spec_index_rg = param_dict['spec_index_rg'])
    ax = subplot(111, yscale='log')
    plot(el, cl_dg_po, color = 'orangered', ls = '-', label = r'Dusty galaxies: Poisson')
    plot(el, cl_dg_clus, color = 'orangered', ls = '-.', label = r'Dusty galaxies: Clustered')
    plot(el, cl_radio, color = 'lime', ls = '-', label = r'Radio galaxies: Poisson')
    plot(el, cl_dg_po + cl_radio, color = 'black', ls = '-', label = r'Poisson: Dusty + radio galaxies')
    axhline(5e-6, ls = '--', color = 'k', label = r'S4 prediction')
    legend(loc = 1, fancybox = 1, fontsize = 8)
    xlabel(r'Multipole $\ell$')
    ylabel(r'C$_{\ell}$ [$\mu K^{2}$]')
    title(r'Point source power @ 145  GHZ')
    ylim(1e-7, 1e-4)
    xlim(100, 5000)
if (1):
    el, cl_tsz_cib = fg.get_foreground_power_spt('tSZ-CIB', freq1 = param_dict['freq0'], freq2 = param_dict['freq0'])
    el, cl_dg_po, cl_dg_clus = fg.get_cl_dust(145, 145, freq0 = param_dict['freq0'], fg_model = param_dict['fg_model'], spec_index_dg_po = param_dict['spec_index_dg_po'], spec_index_dg_clus = param_dict['spec_index_dg_clus'], Tcib = param_dict['Tcib'])
    el, cl_radio = fg.get_cl_radio(145, 145, freq0 = param_dict['freq0'], fg_model = param_dict['fg_model'], spec_index_rg = param_dict['spec_index_rg'])
    el, cl_tsz = fg.get_cl_tsz(145, 145, freq0 = param_dict['freq0'], fg_model = param_dict['fg_model'])
    dls_fac = el * (el+1)/2/np.pi
    cl_point_source = cl_radio + cl_dg_po
    tr, tc = 10, 1
    rspan = 5
    rspan2 = tr - rspan
    clf()
    subplots_adjust(hspace=0.1)
    ax = subplot2grid((tr,tc), (0, 0), rowspan = rspan, yscale='log')
    plot(el, cl_point_source * dls_fac, color = 'green', ls = '-', label = r'Poisson point sources')
    plot(el, cl_tsz * dls_fac, color = 'navy', ls = '-', label = r'tSZ')
    plot(el, cl_dg_clus * dls_fac, color = 'orangered', ls = '-', label = r'Dusty galaxies: Clustered')
    plot([], [], 'k-', label=r'DRAFT')
    if (1): #WAFTT
        waftt_cl_dg_clus = np.loadtxt('data/WAFTT/Cls_cibc_145x145.txt')
        waftt_cl_tsz = np.loadtxt('data/WAFTT/Cls_tsz_145x145.txt')
        waftt_el = np.arange( len(waftt_cl_tsz) )
        waftt_cl_point_source = np.tile(5e-6, len(waftt_el))
        waftt_dls_fac = waftt_el * (waftt_el+1)/2/np.pi
        plot(waftt_el, waftt_cl_point_source * waftt_dls_fac, color = 'green', ls = '-.')#, label = r'WAFTT: Poisson')
        plot(waftt_el, waftt_cl_tsz * waftt_dls_fac, color = 'navy', ls = '-.')#, label = r'WAFTT: tSZ')
        plot(waftt_el, waftt_cl_dg_clus * waftt_dls_fac, color = 'orangered', ls = '-.')#, label = r'WAFTT: Dusty galaxies: Clustered')
        plot([], [], 'k-.', label=r'WAFTT')
    #plot(el, -cl_tsz_cib * dls_fac, 'g-')
    #plot(el, -cl_tsz_cib_90 * dls_fac, 'm-')
    legend(loc = 2, fancybox = 1, fontsize = 6)
    ylabel(r'D$_{\ell}$ [$\mu K^{2}$]')
    title(r'Point sources / tSZ power @ 145 GHZ')
    ylim(1e-1, 1e3)
    xlim(100, 5000)
    setp(ax.get_xticklabels(which = 'both'), visible=False)

    ax = subplot2grid((tr,tc), (rspan+1, 0), rowspan = rspan2)
    el_inds = np.arange(5000)
    plot(el_inds, cl_tsz[el_inds]/waftt_cl_tsz[el_inds], color = 'navy', lw = 0.5, alpha = 0.5)
    plot(el_inds, cl_point_source[el_inds]/waftt_cl_point_source[el_inds], color = 'green', lw = 0.5, alpha = 0.5)
    plot(el_inds, cl_dg_clus[el_inds]/waftt_cl_dg_clus[el_inds], color = 'orangered', lw = 0.5, alpha = 0.5)

    cl_total_draft  = cl_tsz + cl_point_source + cl_dg_clus
    cl_total_waftt  = waftt_cl_tsz + waftt_cl_point_source + waftt_cl_dg_clus
    labval = r'CIB + radio + tSZ'
    plot(el_inds, cl_total_draft[el_inds]/cl_total_waftt[el_inds], color = 'black', lw = 1., label = labval)

    legend(loc = 2, fancybox = 1, fontsize = 6)
    xlabel(r'Multipole $\ell$')
    ylabel(r'DRAFT/WAFTT', fontsize = 8)
    ylim(0.5, 1.5)
    xlim(100, 5000)
    axhline(1., lw = 0.5)


# In[50]:


#colordic = {27:'indigo', 39:'royalblue', 93: 'lightgreen', 145: 'darkgreen', 225: 'goldenrod', 278: 'darkred'}
colordic = {}
dummy_freqarr = [20, 27, 39, 93, 145, 225, 278]
colorarr = [cm.jet(int(d)) for d in np.linspace(0, 255, len(dummy_freqarr) )]
for fcntr, f in enumerate( dummy_freqarr ):
    colordic[f] = colorarr[fcntr]
colordic[27] = 'navy'
colordic[39] = 'royalblue'
colordic[93] = 'lightseagreen'
colordic[145] = 'darkgreen'


# In[51]:


'''
#get the residual power now
weights_dic, cl_residual = {}, {}
for which_spec in which_spec_arr:
    print(which_spec)
    cl_residual[which_spec], weights_dic[which_spec] = ilc.residual_power(param_dict, freqarr, el, cl_dic, which_spec, final_comp = final_comp, freqcalib_fac = freqcalib_fac, return_weights = 1)
'''


# In[52]:


'''
#plot weights now
clf()
fig = figure(figsize=(6, 6))
for cntr, which_spec in enumerate( which_spec_arr ): 
    #if which_spec == 'TE': continue
    ax = subplot ( len(which_spec_arr), 1, cntr+1)
    if which_spec == 'TE':
        tot_teiter = 1 ##2
    else:
        tot_teiter = 1 
    
    for teiter in range(tot_teiter):
        shift = teiter * len(freqarr)
        if teiter == 0:
            lsval = '-'
            lwval = 1.
        else:
            lsval = ':'
            lwval = 0.5
        for frqcntr, freq in enumerate( freqarr ):            
            weights_arr = weights_dic[which_spec][frqcntr+shift]
            plot(weights_arr, color = colordic[freq], ls = lsval, lw = lwval, label = r'%s' %(freq))
        plot(np.sum(weights_dic[which_spec], axis = 0), 'k--', lw = 0.5, label = r'Sum')    
    axhline(lw=0.3);
    if cntr == len(which_spec):
        xlabel(r'Multipole $\ell$');
    else:
        setp(ax.get_xticklabels(which = 'both'), visible=False)
    if cntr == 0:
        legend(loc = 4, ncol = 5, fontsize = 6)
        
    ylabel(r'Weight $W_{\ell}$')
    ylim(-2., 2.);
show()#;sys.exit()
'''


# In[53]:


#get the residual power now
cl_residual_arr, weights_arr = ilc.residual_power_new(param_dict, freqarr, el, cl_dic, final_comp = final_comp, freqcalib_fac = freqcalib_fac, return_weights = 1)


# In[54]:


cl_residual, weights_dic = {}, {}
for which_spec in which_spec_arr:
    if which_spec == 'TT':
       cl_res = cl_residual_arr[0]
       weights = weights_arr[:nc, 0]
    elif which_spec == 'EE':
       cl_res = cl_residual_arr[1]
       weights = weights_arr[nc:, 1]
    elif which_spec == 'TE':
       cl_res = cl_residual_arr[2]
       weights = np.asarray( [weights_arr[nc:, 0], weights_arr[:nc, 1]] )
    cl_residual[which_spec], weights_dic[which_spec] = cl_res, weights
print(cl_residual.keys())


# In[55]:


#plot and results file name
freqarr_str = '-'.join( np.asarray( freqarr ).astype(str) )
which_spec_arr_str = '-'.join( np.asarray( which_spec_arr ).astype(str) )
#opfname = 'results/galactic_sims/S4_ilc_20204020_galaxy%s_%s.npy' %(include_gal, freqarr_str)
#opfname = 'results/galactic_sims/S4_ilc_zonca_sims_20204028_galaxy%s_%s_%s.npy' %(include_gal, freqarr_str, which_spec_arr_str)
#parent_folder = 'results/20200610'
parent_folder = 'results/20200701'
opfname = '%s/%s_ilc_galaxy%s_%s_%s.npy' %(parent_folder, expname, include_gal, freqarr_str, which_spec_arr_str)

if not corr_noise:
    opfname = opfname.replace('.npy', '_nocorrnoise.npy')
    
if s4like_mask:
    opfname = opfname.replace(parent_folder, '%s/s4like_mask/' %(parent_folder))

if s4like_mask_v2:
    opfname = opfname.replace(parent_folder, '%s/s4like_mask_v2/' %(parent_folder))

if include_gal:
    opfname = opfname.replace('.npy', '_galmask%s.npy' %(which_gal_mask))

if remove_atm:
    opfname = opfname.replace('.npy', '_noatmnoise.npy')

if include_gal:    
    cl_gal_folder = param_dict['cl_gal_folder']
    if cl_gal_folder.find('CUmilta')>-1:
        opfname = opfname.replace('.npy', '_CU.npy')
    else:
        opfname = opfname.replace('.npy', '_AZ.npy')

try:
    param_dict['cl_gal_dic_sync_fname_forced']
    opfname = opfname.replace('.npy', '_forcingsynctoCU.npy')
except:
    pass

#plname = opfname.replace('.npy', '.png').replace('S4_ilc', 'plot_S4_ilc')
if s4like_mask:
    plname = opfname.replace('.npy', '.png').replace('%s/s4like_mask/' %(parent_folder), '%s/s4like_mask/plots/' %(parent_folder))
elif s4like_mask_v2:
    plname = opfname.replace('.npy', '.png').replace('%s/s4like_mask_v2/' %(parent_folder), '%s/s4like_mask_v2/plots/' %(parent_folder))
else:
    plname = opfname.replace('.npy', '.png').replace(parent_folder, '%s/plots/' %(parent_folder))
plfolder = '/'.join(plname.split('/')[:-1])
os.system('mkdir -p %s' %(plfolder))
os.system('mkdir -p %s' %(plfolder))
print(opfname)
print(plname)


# In[56]:


freq0, lmax = param_dict['freq0'], param_dict['lmax']
if include_gal:
    foregrounds_to_plot = ['kSZ', 'tSZ', 'DG-Po', 'DG-Cl', 'RG', 'dust', 'sync']
    pol_foregrounds_to_plot = ['dust', 'sync']
else:
    foregrounds_to_plot = ['kSZ', 'tSZ', 'DG-Po', 'DG-Cl', 'RG']
    pol_foregrounds_to_plot = []

#CAMB output for plotting
#camb_file = param_dict['Dlfile_len']
camb_file = '%s/%s' %(param_dict['data_folder'], param_dict['Dlfile_len'])
el_camb = np.loadtxt(camb_file, usecols = [0])
dl_camb = np.loadtxt(camb_file, usecols = [1,2,3,4])

cl_camb = ( Tcmb**2. * dl_camb * 2 * np.pi ) / ( el_camb[:,None] * (el_camb[:,None] + 1) )
cl_camb *= 1e12
cl_TT, cl_EE, cl_BB, cl_TE = cl_camb.T

clf(); 
fig = figure(figsize = (6,3))
fsval = 8
lwval = 0.75
plot_weights = 0
xmin, xmax = 100, 5000
if plot_weights:
    tr, tc = 6, len(which_spec_arr)
    subplots_adjust(wspace=0.1, hspace = 0.1)
    #first plot weights
    rspan, cspan = 2, 1
    curr_row = 0
    for cntr, which_spec in enumerate( which_spec_arr ):
        ax = subplot2grid((tr,tc), (curr_row, cntr), rowspan = rspan, colspan = cspan)#, xscale = 'log')#, yscale = 'log')
        if which_spec == 'TE':
            continue
            tot_teiter = 2
            for teiter in range(tot_teiter):
                if teiter == 0:
                    lsval = '-'
                else:
                    lsval = ':'
                shift = teiter * len(freqarr)
                for frqcntr, freq in enumerate( freqarr ):
                    #plot(weights_dic[which_spec][teiter][frqcntr], color = colordic[freq], label = r'%s' %(freq), lw = lwval, ls = lsval)
                    plot(weights_dic[which_spec][frqcntr+shift], color = colordic[freq], label = r'%s' %(freq), lw = lwval, ls = lsval)
                plot(np.sum(weights_dic[which_spec][teiter], axis = 0), 'k--', label = r'Sum', lw = lwval, ls = lsval)        
        else:
            for frqcntr, freq in enumerate( freqarr ):
                plot(weights_dic[which_spec][frqcntr], color = colordic[freq], label = r'%s' %(freq), lw = lwval)
            plot(np.sum(weights_dic[which_spec], axis = 0), 'k--', label = r'Sum', lw = lwval)
        axhline(lw=0.3);
        #xlabel(r'Multipole $\ell$');
        setp(ax.get_xticklabels(which = 'both'), visible=False)
        if cntr == 0:
            ylabel(r'Weight $W_{\ell}$')
            legend(loc = 1, fontsize = 5, ncol = 4, handlelength = 2., handletextpad = 0.1)
        else:
            setp(ax.get_yticklabels(which = 'both'), visible=False)
        ylim(-3., 3.);
        xlim(xmin, xmax);
        for label in ax.get_xticklabels(): label.set_fontsize(fsval)
        for label in ax.get_yticklabels(): label.set_fontsize(fsval)        

        title(r'%s' %(which_spec))#, fontsize = 10)

    curr_row = rspan
    rspan = tr - rspan
for cntr, which_spec in enumerate( which_spec_arr ):
    if plot_weights:
        #ax = subplot(1,2,cntr+1, xscale = 'log', yscale = 'log')
        ax = subplot2grid((tr,tc), (curr_row, cntr), rowspan = rspan, colspan = cspan, xscale = 'log', yscale = 'log')
    else:
        ax = subplot(1, len(which_spec_arr), cntr+1, xscale = 'log', yscale = 'log')
    plot(el, cl_residual[which_spec], 'black', lw = 2., label = r'Residual')
    if which_spec == 'TT':
        plot(el_camb, cl_TT, 'gray', lw = 1., label = r'TT')
        '''
        cl_fg = np.zeros(len(el))
        for curr_fg in foregrounds_to_plot:
            if curr_fg == 'dust':
                el, cl_curr_fg = fg.get_cl_galactic(param_dict, 'dust', 145, 145, 'TT', bl_dic = bl_dic, el = el)
            elif curr_fg == 'sync':
                el, cl_curr_fg = fg.get_cl_galactic(param_dict, 'sync', 145, 145, 'TT', bl_dic = bl_dic, el = el)
            else:
                el_, cl_curr_fg = fg.get_foreground_power_spt(curr_fg, freq1 = freq0, lmax = lmax)
            #plot(el, cl_curr_fg, lw = 0.5, ls = '--', label = r'150: %s' %(curr_fg), alpha = 0.4)
            cl_fg += cl_curr_fg
        plot(el, cl_fg, lw = 5., ls = '--', label = r'150: All foregrounds', alpha = 1.)
        '''
    elif which_spec == 'EE':
        plot(el_camb, cl_EE, 'gray', lw = 0.5)#, label = r'EE')
        '''
        for curr_fg in pol_foregrounds_to_plot:
            if curr_fg == 'dust':
                el, cl_curr_fg = fg.get_cl_galactic(param_dict, 'dust', 145, 145, 'EE', bl_dic = bl_dic, el = el)
            elif curr_fg == 'sync':
                el, cl_curr_fg = fg.get_cl_galactic(param_dict, 'sync', 145, 145, 'EE', bl_dic = bl_dic, el = el)
            plot(el, cl_curr_fg, lw = 0.5, ls = '--', label = r'150: %s' %(curr_fg), alpha = 0.4)
        '''
    elif which_spec == 'TE':
        plot(el_camb, cl_TE, 'gray', ls = '-', lw = 0.5)#, label = r'TE')        
        plot(el_camb, abs( cl_TE ), 'gray', ls = '--', lw = 0.5) 
        '''
        for curr_fg in pol_foregrounds_to_plot:
            if curr_fg == 'dust':
                el, cl_curr_fg = fg.get_cl_galactic(param_dict, 'dust', 145, 145, 'EE', bl_dic = bl_dic, el = el)
            elif curr_fg == 'sync':
                el, cl_curr_fg = fg.get_cl_galactic(param_dict, 'sync', 145, 145, 'EE', bl_dic = bl_dic, el = el)
            plot(el, cl_curr_fg, lw = 0.5, ls = '--', label = r'150: %s' %(curr_fg), alpha = 0.4)
        '''
    #for freq in freqarr:
    #    plot(el, cl_dic[which_spec][(freq,freq)], color = colordic[freq], lw = 0.5, ls = '-', label = r'%s' %(freq), alpha = 1.)        
    mv_comb_arr = []
    for freq in freqarr:
        if which_spec == 'TT':
            nl = nl_dic['T'][(freq, freq)]
        elif which_spec == 'EE':
            nl = nl_dic['P'][(freq, freq)]
        elif which_spec == 'TE':
            nl = nl_dic['T'][(freq, freq)] * 0.
        plot(el, nl, color = colordic[freq], lw = lwval, ls = '--', label = r'Noise: %s' %(freq))#, alpha = 0.5)
        mv_comb_arr.append(1./nl)
    mv_comb = 1./(np.sum(mv_comb_arr, axis = 0))
    plot(el, mv_comb, color = 'k', lw = 0.5, ls = '--', label = r'Noise: MV')
    #legend(loc=3, fancybox=1, ncol = 4, fontsize = 6);
    
    xlim(xmin, xmax);
    ylim(1e-7,1e-1);
    xlabel(r'Multipole $\ell$')
    if cntr == 0: 
        ylabel(r'$C_{\ell}$')
        legend(loc = 3, fontsize = 6, ncol = 2, handlelength = 2., handletextpad = 0.1)
    else:
        pass#setp(ax.get_yticklabels(which = 'both'), visible=False)
    for label in ax.get_xticklabels(): label.set_fontsize(fsval)
    for label in ax.get_yticklabels(): label.set_fontsize(fsval)
        
#tit = 'Galaxy = %s; Mask = %s; Bands = %s' %(include_gal, param_dict['which_gal_mask'], str(freqarr))
if remove_atm:
    tit = 'Galaxy = %s; Mask = %s; Bands = %s; no 1/f' %(include_gal, param_dict['which_gal_mask'], str(freqarr))
else:
    if include_gal:
        tit = 'Galaxy = %s; Mask = %s; Bands = %s' %(include_gal, param_dict['which_gal_mask'], str(freqarr))    
    else:
        tit = 'Galaxy = %s; Mask = N/A; Bands = %s' %(include_gal, str(freqarr))    
if (0):#not corr_noise:
    tit = '%s; No corr. noise' %(tit)
suptitle(r'%s' %tit, x = 0.53, y = 1.)
######savefig(plname)
show()


# In[57]:


#first plot weights
fig = figure(figsize = (6,3))
rspan, cspan = 2, 1
curr_row = 0
for cntr, which_spec in enumerate( which_spec_arr ):
    ax = subplot(1,2,cntr+1)#, xscale = 'log')#, yscale = 'log')
    if which_spec == 'TE':
        continue
        tot_teiter = 2
        for teiter in range(tot_teiter):
            if teiter == 0:
                lsval = '-'
            else:
                lsval = ':'
            shift = teiter * len(freqarr)
            for frqcntr, freq in enumerate( freqarr ):
                #plot(weights_dic[which_spec][teiter][frqcntr], color = colordic[freq], label = r'%s' %(freq), lw = lwval, ls = lsval)
                plot(weights_dic[which_spec][frqcntr+shift], color = colordic[freq], label = r'%s' %(freq), lw = lwval, ls = lsval)
            plot(np.sum(weights_dic[which_spec][teiter], axis = 0), 'k--', label = r'Sum', lw = lwval, ls = lsval)        
    else:
        for frqcntr, freq in enumerate( freqarr ):
            plot(weights_dic[which_spec][frqcntr], color = colordic[freq], label = r'%s' %(freq), lw = lwval)
        plot(np.sum(weights_dic[which_spec], axis = 0), 'k--', label = r'Sum', lw = lwval)
    axhline(lw=0.3);
    xlabel(r'Multipole $\ell$');
    #setp(ax.get_xticklabels(which = 'both'), visible=False)
    if cntr == 0:
        ylabel(r'Weight $W_{\ell}$')
        legend(loc = 3, fontsize = 5, ncol = 4, handlelength = 2., handletextpad = 0.1)
    else:
        setp(ax.get_yticklabels(which = 'both'), visible=False)
    ylim(-1., 2.);
    xlim(xmin, xmax);
    for label in ax.get_xticklabels(): label.set_fontsize(fsval)
    for label in ax.get_yticklabels(): label.set_fontsize(fsval)        

    title(r'%s' %(which_spec))#, fontsize = 10)


# In[ ]:





# In[60]:


cl_gal_dic_dust_fname = param_dict['cl_gal_dic_dust_fname']
try:
    cl_gal_folder = param_dict['cl_gal_folder']
    cl_gal_dic_dust_fname = '%s/%s' %(cl_gal_folder, cl_gal_dic_dust_fname)
except:
    pass
print(cl_gal_dic_dust_fname)
galdustsims_cl = np.load(cl_gal_dic_dust_fname, allow_pickle=1, encoding = 'latin1').item()
if not include_gal:
    fsky_val = 0.68
else:
    fsky_val = galdustsims_cl['fsky_arr'][param_dict['which_gal_mask']]
opdic = {}
opdic['el'] = el
opdic['cl_residual'] = cl_residual
opdic['freqcalib_fac'] = freqcalib_fac
opdic['param_dict'] = param_dict
opdic['fsky_val'] = fsky_val
opdic['which_gal_mask'] = which_gal_mask
opdic['weights'] = weights_dic
#opdic['nl_dic'] = nl_dic
opdic['beam_noise_dic'] = beam_noise_dic
opdic['elknee_dic'] = elknee_dic
np.save(opfname, opdic)
print(opfname)


# In[97]:


farr = ['results/20200701/s4like_mask_v2//s4wide_ilc_galaxy0_27-39-93-145-225-278_TT-EE.npy',    'results/20200701/s4like_mask_v2//s4wide_ilc_galaxy1_27-39-93-145-225-278_TT-EE_galmask2_AZ.npy',    'results/20200701/s4like_mask_v2//s4wide_ilc_galaxy1_27-39-93-145-225-278_TT-EE_galmask5_AZ.npy',]


# In[116]:


miny, maxy = 0., 1.2
ilc_results_dic = {}
for fcntr, f in enumerate( farr ):
    ilc_dic = np.load(f, allow_pickle=1).item()['cl_residual']
    ilc_results_dic[fcntr] = [ ilc_dic['TT'], ilc_dic['EE'] ]

pl_comb_arr = [(0,1), (0,2), (1,2)]
lab_arr = [('gal0', 'galmask2'), ('gal0', 'galmask5'), ('galmask2', 'galmask5')]
colorarr = ['orangered', 'green', 'k']
alphaarr = [0.2, 0.2, 1.]
for cntr, pl_comb in enumerate( pl_comb_arr ):
    p1,p2 = pl_comb
    ax  =subplot(121)
    labval = r'%s/%s' %(lab_arr[cntr][0], lab_arr[cntr][1])
    alphaval = alphaarr[cntr]
    plot(ilc_results_dic[p1][0]/ilc_results_dic[p2][0], color = colorarr[cntr], label = r'%s' %labval, alpha = alphaval)
    ylabel(r'Ratio')
    axhline(1., lw=0.5);xlim(30, 5000);ylim(miny, maxy)
    legend(loc = 1, fancybox = 1, fontsize = 8)
    ax  =subplot(122)
    plot(ilc_results_dic[p1][1]/ilc_results_dic[p2][1], color = colorarr[cntr], alpha = alphaval)
    xlabel(r'Multipole $\ell$')
    axhline(1., lw=0.5);xlim(30, 5000);ylim(miny, maxy)


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:


f1 = 'results/20200610/s4like_mask//s4wide_ilc_galaxy1_27-39-93-145-225-278_TT-EE_galmask3_AZ.npy'
f2 = 'results/20200701/DSR/s4like_mask//s4wide_ilc_galaxy1_27-39-93-145-225-278_TT-EE_galmask3_AZ.npy'
f3 = 'results/20200701/s4like_mask//s4wide_ilc_galaxy1_27-39-93-145-225-278_TT-EE_galmask3_AZ.npy'
ilc_dic1 = np.load(f1, allow_pickle=1).item()['cl_residual']
ilc_dic2 = np.load(f2, allow_pickle=1).item()['cl_residual']
ilc_dic3 = np.load(f3, allow_pickle=1).item()['cl_residual']
tt1, ee1 = ilc_dic1['TT'], ilc_dic1['EE']
tt2, ee2 = ilc_dic2['TT'], ilc_dic2['EE']
tt3, ee3 = ilc_dic3['TT'], ilc_dic3['EE']
tt_dic={1:tt1, 2: tt2, 3:tt3}
ee_dic={1:ee1, 2: ee2, 3:ee3}
pl_comb_arr = [(3,2)] #,(1,2)]
#pl_comb_arr = [(3,1)] #,(1,2)]
"""
tt_arr=[ ee1, ee2, ee3]
#ilc_ee_old = np.load('dummy.npy', allow_pickle=1)
#plot(ilc_ee/ilc_ee_old, lw = 2.)
subplot(121);
plot(tt3/tt2, color = 'k')
#plot(tt2/tt1, color = 'purple')
title(r'TT')
ylim(0.8, 1.2)
subplot(122);
plot(ee3/ee2, color = 'k')
#plot(ee2/ee1, color = 'purple')
title(r'EE')
xlim(0, 5000)
ylim(0.9, 1.1)
"""
for pl_comb in pl_comb_arr:
    p1,p2 = pl_comb
    plot(tt_dic[p1]/tt_dic[p2], color = 'k', label = r'TT')
    plot(ee_dic[p1]/ee_dic[p2], color = 'orangered', label = r'EE')
xlabel(r'Multipole $\ell$')
ylabel(r'Ratio [current/Old]')
axhline(1., lw=0.5)
xlim(30, 5000)
ylim(0.8, 1.2)
legend(loc = 1)


# In[ ]:





# In[ ]:





# In[ ]:




