#!/usr/bin/env python
# coding: utf-8

# In[1]:


from pylab import *
rcParams["figure.facecolor"] = 'white'


# In[2]:


rcParams['figure.dpi'] = 150


# In[3]:


import sys, numpy as np, scipy as sc, warnings, os#, healpy as H
sys.path.append('/Users/sraghunathan/Research/SPTPol/analysis/git/DRAFT/modules/')
import foregrounds as fg, misc

#import matplotlib.cbook
warnings.filterwarnings('ignore',category=RuntimeWarning)
#warnings.filterwarnings('ignore', category=DeprecationWarning)
#warnings.filterwarnings('ignore', category=matplotlib.cbook.mplDeprecation)


# In[19]:


#params
TParr = ['T', 'P']
#which_spec = 'TT'
which_spec = 'EE'
#which_spec = 'BB'
#which_spec = 'TE'
#which_spec = 'EB'

param_dict = {}


#spt3g masks
TParr = ['T', 'P']
nside, lmax = 2048, 6000
data_folder = '/Volumes/data_PHD_WD_babbloo/s4/cmbs4/map_based_simulations/202002_foregrounds_extragalactic_cmb_tophat/4096/'
###nside, lmax = 2048, 2000 ##6000
###data_folder = '/Volumes/data_PHD_WD_babbloo/s4/cmbs4/map_based_simulations/202102_design_tool_input/4096/'
param_dict['cl_gal_dic_dust_fname'] = '%s/dust/0000/spt3g/cls_galactic_sims_dust_nside%s_lmax%s.npy' %(data_folder, nside, lmax)
param_dict['cl_gal_dic_sync_fname'] = '%s/synchrotron/0000/spt3g/cls_galactic_sims_sync_nside%s_lmax%s.npy' %(data_folder, nside, lmax)
param_dict['cl_gal_dic_freefree_fname'] = '%s/freefree/0000/spt3g/cls_galactic_sims_freefree_nside%s_lmax%s.npy' %(data_folder, nside, lmax)
extra_str = 'pySM3/spt3g/'

param_dict['lmax'] = lmax
param_dict['Dlfile_len'] = '/Users/sraghunathan/Research/SPTPol/analysis/git/DRAFT/data/output_planck_r_0.0_2015_cosmo_lensedCls.dat'


# In[20]:


#CAMB output for plotting
camb_file = param_dict['Dlfile_len']
Tcmb = 2.725
el_camb = np.loadtxt(camb_file, usecols = [0])
dl_camb = np.loadtxt(camb_file, usecols = [1,2,3,4])
dl_camb = dl_camb * Tcmb**2. * 1e12

cl_camb = ( dl_camb * 2 * np.pi ) / ( el_camb[:,None] * (el_camb[:,None] + 1) )
cl_TT, cl_EE, cl_BB, cl_TE = cl_camb.T
dl_TT, dl_EE, dl_BB, dl_TE = dl_camb.T


# In[24]:

freqarr = [93, 145]#, 225]
which_gal_mask = 0
use_sed_scaling = False
bl_dic = None

'''
clf()
fig = figure(figsize=(8., 4.))
subplots_adjust(wspace=0.05)
sbpl = 1
plotted = []
for freq1 in freqarr:
    for freq2 in freqarr:

        if (freq1, freq2) in plotted or (freq2, freq1) in plotted: continue
        if freq1 != freq2: continue
            
        print(freq1, freq2)

        el, cl_dust = fg.get_cl_galactic(param_dict, 'dust', freq1, freq2, which_spec = which_spec, which_gal_mask = which_gal_mask, use_sed_scaling = use_sed_scaling, bl_dic = bl_dic)
        el, cl_sync = fg.get_cl_galactic(param_dict, 'sync', freq1, freq2, which_spec = which_spec, which_gal_mask = which_gal_mask, use_sed_scaling = use_sed_scaling, bl_dic = bl_dic)
        el, cl_freefree = fg.get_cl_galactic(param_dict, 'freefree', freq1, freq2, which_spec = which_spec, which_gal_mask = which_gal_mask, use_sed_scaling = use_sed_scaling, bl_dic = bl_dic)

        ax = subplot(1,2,sbpl, yscale = 'log')
        dl_fac = (el* (el+1))/2./np.pi
        if which_spec == 'TT':
            plot(el_camb, dl_TT, color = 'gray', label = r'CMB')
        elif which_spec == 'EE':
            plot(el_camb, dl_EE, color = 'gray', label = r'CMB')
        plot(el, dl_fac * cl_dust, color = 'tab:orange', label = r'Dust')
        plot(el, dl_fac * cl_sync, color = 'tab:green', label = r'Sync')
        plot(el, dl_fac * cl_freefree, color = 'tab:red', label = r'FF')
        title(r'%s GHz' %(freq1), fontsize = 10)
        xlim(-50., 800.)
        if which_spec == 'TT':
            ylim(1e-4, 1e4)
        elif which_spec == 'EE':
            ylim(1e-4, 1e2)
        if sbpl == 1:
            legend(loc = 1, fontsize = 10)
            ylabel(r'D$_{\ell}^{\rm %s}$ [$\mu$K$^{2}$]' %(which_spec), fontsize = 12)
        else:
            setp(ax.get_yticklabels(which = 'both'), visible=False)
        sbpl += 1
        xlabel(r'Multipole $\ell$', fontsize = 12)
        plotted.append((freq1, freq2))
suptitle(r'Sim set = 202102 design tool input (Beam = Not deconvolved)', fontsize = 12, y = 1.)
###savefig('pysm3_202102_design_tool_input_spt3gmask_%s.png' %(which_spec), dpi = 200.)
show(); 
sys.exit()
'''



freqarr = [93, 145] ##278]
#alpha_dic = {93: 0.3, 145: 0.6, 225: 0.8, 278: 1.}
alpha_dic = {93: .5, 145: 1., 225: 1., 278: 1.}
lw_dic = {93: 1., 145: 1., 225: 1.5, 278: 2.}
ls_dic = {93: '-.', 145: '-', 225: '--', 278: ':'}

tot_mask_iter = 4
color_dic = {0:'navy', 1: 'green', 2: 'goldenrod', 3: 'darkred'}

tot_mask_iter = 5
color_dic = {0:'navy', 1: 'green', 2: 'goldenrod', 3: 'orangered', 4: 'darkred'}

tot_mask_iter = 4
color_dic = {0:'navy', 1: 'green', 2: 'goldenrod', 3: 'darkred'}
reqd_masks = []
if extra_str.find('s4like_mask_v2')>-1 or extra_str.find('s4delensing')>-1:
    tot_mask_iter = 6
    if extra_str.find('s4delensing')>-1:
        reqd_masks = [0]
    else:
        reqd_masks = [0, 1, 2, 3, 4, 5]
    #color_dic = {0:'navy', 1: 'darkblue', 2: 'royalblue', 3: 'orangered', 4: 'darkred', 5: 'maroon'}
    cmap = cm.jet
    color_dic = {0:cmap(0), 1: cmap(15), 2: cmap(30), 3: cmap(230), 4: cmap(240), 5: cmap(255)}
    freqarr = [145]
    which_spec = 'EE'#'TT' #EE' #TT'
if extra_str.find('spt3g')>-1:
    tot_mask_iter = 4
    reqd_masks = [0, 1, 2, 3]
    #color_dic = {0:'navy', 1: 'darkblue', 2: 'royalblue', 3: 'orangered', 4: 'darkred', 5: 'maroon'}
    cmap = cm.jet
    color_dic = {0: 'navy', 1: 'darkred', 2: 'darkgreen', 3: 'goldenrod'}
    mask_str_dic = {0: 'Winter', 1: 'Summer: el1c-el2c', 2: 'Summer: el1b-el2b', 3: 'Summer: el1-el5'}
    freqarr = [145]
    which_spec = 'TT' #'EE' #'TT'


#xscale_val = 'log' ##None#'log' #None
xscale_val = None
for mask_iter in range(tot_mask_iter):
    if len(reqd_masks)>0:
        if mask_iter not in reqd_masks: continue
    plot_done = []
    for freq1 in freqarr:
        for freq2 in freqarr:
            

            if (freq1, freq2) in plot_done: continue
                
            if freq1 != freq2: continue
                
            try:
                el, cl_dust = fg.get_cl_galactic(param_dict, 'dust', freq1, freq2, which_spec = which_spec, which_gal_mask = mask_iter, use_sed_scaling = False)
            except:
                cl_dust = None


            try:
                el, cl_sync = fg.get_cl_galactic(param_dict, 'sync', freq1, freq2, which_spec = which_spec, which_gal_mask = mask_iter)
            except:
                cl_sync = None
                
            if freq1 == 93:
                delta_ell = 50
            elif freq1 == 145:
                delta_ell = 100                         
            el_to_plot = np.arange(min(el), max(el)+1, delta_ell)

            plot_done.append((freq1, freq2))
            lab = '%s: (%s,%s)' %(mask_iter, freq1, freq2)

            lwval = lw_dic[freq1]
            colorval = color_dic[mask_iter]
            alphaval = alpha_dic[freq1]
            lsval = ls_dic[freq1]
            #lwval = 1.
            lab = None
            if cl_dust is not None:
                cl_dust_to_plot = np.interp(el_to_plot, el, cl_dust)
                ax = subplot(1,2,1, yscale = 'log', xscale=xscale_val)
                plot(el, cl_dust, label = lab, lw = 0.5, linestyle = lsval, color = colorval, alpha = alphaval)
                #plot(el_to_plot, cl_dust_to_plot, label = lab, lw = lwval, linestyle = lsval, color = colorval, alpha = alphaval)
            if cl_sync is not None:
                cl_sync_to_plot = np.interp(el_to_plot, el, cl_sync)                            
                ax = subplot(1,2,2, yscale = 'log', xscale=xscale_val)
                plot(el, cl_sync, label = lab, lw = 0.5, linestyle = lsval, color = colorval, alpha = alphaval)
                #plot(el_to_plot, cl_sync_to_plot, label = lab, lw = lwval, linestyle = lsval, color = colorval, alpha = alphaval)

for cntr, which_comp in enumerate( ['dust', 'sync'] ):
    
    ax = subplot(1,2,cntr+1, yscale = 'log', xscale=xscale_val)
    if which_spec == 'TT':
        plot(el_camb, cl_TT, 'gray', lw = 2., label = r'%s' %(which_spec), alpha = 0.3, zorder = -1000.)
    elif which_spec == 'EE':
        plot(el_camb, cl_EE, 'gray', lw = 2., label = r'%s' %(which_spec), alpha = 0.3, zorder = -1000.)
    elif which_spec == 'TE':
        plot(el_camb, cl_TE, 'gray', lw = 2., label = r'%s' %(which_spec), alpha = 0.3, zorder = -1000.)
        plot(el_camb, abs( cl_TE ), 'gray', ls = '--', lw = 2., alpha = 0.3, zorder = -1000.)        
    elif which_spec == 'BB':
        plot(el_camb, cl_BB, 'gray', lw = 2., label = r'%s' %(which_spec), alpha = 0.3, zorder = -1000.)
    if cntr == 1: ##0:

        galdustsims_cl = np.load(param_dict['cl_gal_dic_dust_fname'], allow_pickle=1, encoding = 'latin1').item()
        fsky_arr = galdustsims_cl['fsky_arr']
        if len(freqarr)>1:
            for freq1 in absfreqarr:
                plot([], [], color = 'k', alpha = alpha_dic[freq1], ls = ls_dic[freq1], label = r'%s GHz' %(freq1))
        for mask_iter in range(tot_mask_iter):
            if len(reqd_masks)>0:
                if mask_iter not in reqd_masks: continue
            fsky_val = fsky_arr[mask_iter]
            labval = r'Mask %s: f$_{\rm sky} = %.2f$' %(mask_iter, fsky_val)
            if mask_str_dic is not None:
                labval = r'Mask %s (%s): f$_{\rm sky} = %.2f$' %(mask_iter, mask_str_dic[mask_iter], fsky_val)
            plot([], [], color = color_dic[mask_iter], label = labval)
            
        #legend(loc = 3, ncol = 2, fontsize = 8)
        legend(loc = 1, ncol = 1, fontsize = 8)
    xlabel(r'Multipole $\ell$', fontsize = 14)
    if cntr == 0:
        ylabel(r'$C_{\ell}\ [\mu {\rm K}^{2}]$', fontsize = 14)
    else:
        setp(ax.get_yticklabels(which = 'both'), visible=False)
    extra_str_tmp = extra_str.split('/')[0].replace('_', '\_')
    if len(freqarr)>1:
        title(r'%s: %s (%s)' %(which_spec, which_comp.capitalize(), extra_str_tmp), fontsize = 14)
    else:
        title(r'%s: %s (%s): %s GHz' %(which_spec, which_comp.capitalize(), extra_str_tmp, freqarr[0]), fontsize = 14)
    xmax = 5000#7100 #lmax+100
    xlim(-100, xmax); 
    ylim(1e-10, 1e4)

    #xlim(None, 200); ylim(1e-6, 1e-2)
#suptitle('Sim set = %s' %(which_sim_set.replace('_','\_')), fontsize = 12, y = 1.02)
show();sys.exit()
plfolder = 'reports/galactic_sims/dust_sync_spectra/%s/' %(extra_str)
os.system('mkdir -p %s' %(plfolder))
#savefig('%s/dust_sync_%s.pdf' %(plfolder, which_spec))


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




