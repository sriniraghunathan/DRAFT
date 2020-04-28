#!/usr/bin/env python
# coding: utf-8

# In[123]:

'''
get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '2')

#%pylab notebook
get_ipython().run_line_magic('matplotlib', 'inline')
'''
from pylab import *
from matplotlib import rc;rc('text', usetex=True);rc('font', weight='bold');matplotlib.rcParams['text.latex.preamble'] = [r'\boldmath']
rcParams['figure.dpi'] = 150
rcParams["figure.facecolor"] = 'white'
try:
    import os
    rc('text.latex',preamble=r'\usepackage{%s/apjfonts}' %(str(os.getcwd())))
except:
    pass


# In[124]:


import argparse, sys, numpy as np, scipy as sc, warnings, os
#sys.path.append('/Users/sraghunathan/Research/SPTPol/analysis/git/tools/')
#or look into https://github.com/sriniraghunathan/tools
import flatsky, tools, misc
import ilc, foregrounds as fg

#import matplotlib.cbook
warnings.filterwarnings('ignore',category=RuntimeWarning)
#warnings.filterwarnings('ignore', category=DeprecationWarning)
#warnings.filterwarnings('ignore', category=matplotlib.cbook.mplDeprecation)


# In[125]:


#some constants
h=6.62607004e-34 #Planck constant in m2 kg / s
k_B=1.38064852e-23 #Boltzmann constant in m2 kg s-2 / K-1
Tcmb = 2.73 #Kelvin


# In[126]:


#params
paramfile = 'params.ini'
#freqarr = [90, 150, 220, 270]

#S4 specs
specs_dic = {
#freq: [beam_arcmins, white_noise_T, elknee_T, alphaknee_T, whitenoise_P, elknee_P, alphaknee_P] 
#20: [10.0, None, None, None, None, None, None],
27: [7.4, 21.8, 471., 3.5, 30.8, 700, 1.4],
39: [5.1, 12.4, 428., 3.5, 17.6, 700, 1.4], 
93: [2.2, 2.0, 2154., 3.5, 2.9, 700, 1.4],
145: [1.4, 2.0, 4364., 3.5, 2.8, 700, 1.4],
225: [1.0, 6.9, 7334., 3.5, 9.8, 700, 1.4],
278: [0.9, 16.7, 7308., 3.5, 23.6, 700, 1.4],
}
'''
specs_dic = {
#freq: [beam_arcmins, white_noise_T, elknee_T, alphaknee_T, whitenoise_P, elknee_P, alphaknee_P] 
#20: [10.0, None, None, None, None, None, None],
#27: [7.4, 21.8, 471., 0., 30.8, 700, 1.4],
#39: [5.1, 12.4, 428., 0., 17.6, 700, 1.4], 
93: [2.2, 2.0, 2154., 0., 2.9, 700, 1.4],
145: [1.4, 2.0, 4364., 0., 2.8, 700, 1.4],
225: [1.0, 6.9, 7334., 0., 9.8, 700, 1.4],
278: [0.9, 16.7, 7308., 0., 23.6, 700, 1.4],
}
'''


freqarr = sorted( specs_dic.keys() )

corr_noise_bands = {20: [20], 27:[39], 39:[27], 93:[145], 145:[93], 225: [278], 278: [225]}
rho = 0.9

freqcalib_fac = None
final_comp = 'cmb'
TParr = ['T', 'P']


# In[127]:


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


# In[128]:


# read and store param dict
param_dict = misc.fn_get_param_dict(paramfile)
el = np.arange(param_dict['lmax'])
include_gal = param_dict['include_gal']
#print(param_dict['which_gal_mask'])



# In[129]:


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


# In[130]:


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


# In[131]:


#get beams
bl_dic = misc.get_beam_dic(freqarr, beam_noise_dic['T'], param_dict['lmax'])
print(bl_dic.keys())
if (1):
    for freq in freqarr:
        plot(bl_dic[freq], label = freq)
    legend(loc = 1)


# In[132]:


#get the CMB, noise, and foreground covriance
try:
    ignore_fg = param_dict['ignore_fg']
except:
    ignore_fg = []

ignore_fg.append(final_comp.lower()) #the required component need not go into the covariance matrix.
print(ignore_fg)

#freqarr = [93]
#param_dict['which_gal_mask'] = 0
cl_dic = {}
for TP in TParr:
    if TP == 'T':
        el, cl_dic[TP] = ilc.get_analytic_covariance(param_dict, freqarr,                 nl_dic = nl_dic['T'], ignore_fg = ignore_fg, include_gal = include_gal, bl_dic = bl_dic)
    elif TP == 'P':
        el, cl_dic[TP] = ilc.get_analytic_covariance                    (param_dict, freqarr, nl_dic = nl_dic['P'], ignore_fg = ignore_fg, pol = 1,                     pol_frac_per_cent_dust = param_dict['pol_frac_per_cent_dust'],                     pol_frac_per_cent_radio = param_dict['pol_frac_per_cent_radio'],                     pol_frac_per_cent_tsz = param_dict['pol_frac_per_cent_tsz'],                     pol_frac_per_cent_ksz = param_dict['pol_frac_per_cent_ksz'],                     include_gal = include_gal, bl_dic = bl_dic)
print(cl_dic.keys(), cl_dic['T'].keys())


# In[133]:


#get the residual power now
cl_residual = {}
for TP in TParr:
    cl_residual[TP] = ilc.residual_power(param_dict, freqarr, el, cl_dic[TP], final_comp = final_comp, freqcalib_fac = freqcalib_fac)


# In[134]:


freq0, lmax = param_dict['freq0'], param_dict['lmax']
foregrounds_to_plot = ['kSZ', 'tSZ', 'DG-Po', 'DG-Cl', 'RG', 'dust']


#CAMB output for plotting
camb_file = param_dict['Dlfile_len']
el_camb = np.loadtxt(camb_file, usecols = [0])
dl_camb = np.loadtxt(camb_file, usecols = [1,2,3,4])

cl_camb = ( Tcmb**2. * dl_camb * 2 * np.pi ) / ( el_camb[:,None] * (el_camb[:,None] + 1) )
cl_camb *= 1e12
cl_TT, cl_EE, cl_BB, cl_TE = cl_camb.T


clf(); 
subplots_adjust(wspace=0.1)
for cntr, TP in enumerate( TParr ):
    ax = subplot(1,2,cntr+1, xscale = 'log', yscale = 'log')
    plot(el, cl_residual[TP], 'black', lw = 2., label = r'Residual')
    if TP == 'T':
        plot(el_camb, cl_TT, 'gray', lw = 1., label = r'TT')
        cl_fg = np.zeros(len(el))
        for curr_fg in foregrounds_to_plot:
            if curr_fg == 'dust':
                el, cl_curr_fg = fg.get_cl_galactic(param_dict, 'dust', 145, 145, 'TT', bl_dic = bl_dic, el = el)
            elif curr_fg == 'sync':
                pass
            else:
                el_, cl_curr_fg = fg.get_foreground_power_spt(curr_fg, freq1 = freq0, lmax = lmax)
            plot(el, cl_curr_fg, lw = 0.5, ls = '--', label = r'150: %s' %(curr_fg), alpha = 0.4)
            cl_fg += cl_curr_fg
        #plot(el, cl_fg, lw = 0.5, ls = '--', label = r'150: All foregrounds', alpha = 1.)
        plot(el, cl_dic['T'][(145,145)], color = 'cyan', lw = 1., ls = '-', label = r'150: Full covariance', alpha = 1.)        
    elif TP == 'P':
        plot(el_camb, cl_EE, 'gray', lw = 0.5, label = r'EE')
        plot(el_camb, cl_TE, 'navy', ls = '-', lw = 0.5, label = r'TE')        
        plot(el_camb, abs( cl_TE ), 'navy', ls = '--', lw = 0.5) 
    for freq in freqarr:
        plot(el, nl_dic[TP][(freq, freq)], lw = 0.5, ls = '-', label = r'Noise: %s' %(freq))#, alpha = 0.5)
    legend(loc=1, fancybox=1, ncol = 2, fontsize = 6);
    xlim(20,lmax);ylim(1e-8,1e6);
    xlabel(r'Multipole $\ell$')
    if cntr == 0: 
        ylabel(r'$C_{\ell}$')
    else:
        setp(ax.get_yticklabels(which = 'both'), visible=False)

suptitle('Galaxy = %s; Mask = %s; Bands = %s' %(include_gal, param_dict['which_gal_mask'], str(freqarr)))
show(); sys.exit()


# In[ ]:


opfname = 'results/S4_ilc_20204015.npy'
opdic = {}
opdic['el'] = el
opdic['cl_residual'] = cl_residual
opdic['freqcalib_fac'] = freqcalib_fac
opdic['param_dict'] = param_dict
#opdic['nl_dic'] = nl_dic
opdic['beam_noise_dic'] = beam_noise_dic
opdic['elknee_dic'] = elknee_dic
####np.save(opfname, opdic)


# In[ ]:




