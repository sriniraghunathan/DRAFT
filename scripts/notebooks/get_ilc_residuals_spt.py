#!/usr/bin/env python
# coding: utf-8

# In[1]:


'''
get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '2')

#%pylab notebook
get_ipython().run_line_magic('matplotlib', 'inline')
'''
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
import flatsky, misc, exp_specs_for_spt as exp_specs
import ilc, foregrounds as fg

#import matplotlib.cbook
warnings.filterwarnings('ignore',category=RuntimeWarning)
#warnings.filterwarnings('ignore', category=DeprecationWarning)
#warnings.filterwarnings('ignore', category=matplotlib.cbook.mplDeprecation)


parser = argparse.ArgumentParser(description='')
parser.add_argument('-expname', dest='expname', action='store', help='expname', type=str, required=True)
parser.add_argument('-final_comp', dest='final_comp', action='store', help='final_comp', type=str, required=True)
parser.add_argument('-null_comp', dest='null_comp', action='store', help='null_comp', nargs='+', default=None)
parser.add_argument('-use_websky_cib', dest='use_websky_cib', action='store', help='use_websky_cib', type = int, default=0)
parser.add_argument('-use_545', dest='use_545', action='store', help='use_545', type = int, default=0)
parser.add_argument('-use_sptspire_for_hfbands', dest='use_sptspire_for_hfbands', action='store', help='use_sptspire_for_hfbands', type = int, default=0)
parser.add_argument('-split_cross', dest='split_cross', action='store', help='split_cross', type = int, default=0)


args = parser.parse_args()
args_keys = args.__dict__
for kargs in args_keys:
    param_value = args_keys[kargs]
    if isinstance(param_value, str):
        cmd = '%s = "%s"' %(kargs, param_value)
    else:
        cmd = '%s = %s' %(kargs, param_value)
    exec(cmd)
# In[4]:

if null_comp is None or null_comp == 'None' or null_comp == ['None']:
    null_comp = None

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
param_dict['include_gal'] = 0
include_gal = param_dict['include_gal'] ##1
if not include_gal:
    param_dict['which_gal_mask'] = 3
which_gal_mask = param_dict['which_gal_mask']
s4like_mask = param_dict['s4like_mask']
try:
    remove_atm = param_dict['remove_atm']
except:
    remove_atm = 0    

#param_dict['ignore_fg'] = ['dust']


# In[6]:

"""
if (0):
    #experiment name
    #expname = 'spt3g_TC'##
    #expname = 'spt3g_summer'#TC'#y12'
    #expname = 'spt3g_ZP'##
    #expname = 'spt3g_y12'##
    #expname = 'sptpol_summer'
    #expname = 'sptpol'
    #expname = 'sptsz'
    '''
    expname = 'sptpolultradeep'
    expname = 'sptpolultradeepplus3g'
    expname = 'sptpolplusultradeep'
    expname = 'sptpolplusultradeepplus3g'
    expname = 'sptpolplusultradeepplus3gfull'
    #expname = 'spt3g'
    '''
    expname = 'spt3g'
    #expname = 'spt4_C4'
    #expname = 'sptpolplusultradeepplus3g'
    #expname = 'sptpolplusultradeepplus3gplusherschel'

    final_comp = 'y' #'cmb' #'tsz'
    null_comp = None

    # In[7]:
"""

remove_atm = 0
specs_dic, corr_noise_bands, rho, corr_noise, cib_corr_coeffs = exp_specs.get_exp_specs(expname, remove_atm = remove_atm)
if use_websky_cib:
    use_545 = 1

if use_545:
    if 600 in specs_dic:
        print('\n\n\t\tModifying 600 GHz to 545 GHz currently when using websky\n\n\n')
        specs_dic[545] = specs_dic[600]
        corr_noise_bands[545] = corr_noise_bands[600]
        specs_dic.pop(600)
        corr_noise_bands.pop(600)
        
freqarr = sorted( specs_dic.keys() )
nc = len( freqarr )
freqcalib_fac = None
#use_websky_cib = 1

if final_comp == 'y':
    TParr = ['T']
    which_spec_arr = ['TT']
elif final_comp == 'cmb':
    TParr = ['T', 'P']
    #which_spec_arr = ['TT', 'EE', 'TE']
    which_spec_arr = ['TT', 'EE']
    ##include_gal = 1

# In[8]:

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

if split_cross:
    noisearr_T = np.asarray(noisearr_T) * np.sqrt(2.)
    noisearr_P = np.asarray(noisearr_P) * np.sqrt(2.)

print(elkneearr_T)


# In[9]:


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

sys.exit()
# In[11]:


#get beams
bl_dic = misc.get_beam_dic(freqarr, beam_noise_dic['T'], param_dict['lmax'])
bl_dic['effective'] = bl_dic[150]
print(bl_dic.keys())
if (0):
    for freq in freqarr:
        plot(bl_dic[freq], label = freq)
    legend(loc = 1)


# In[12]:


#get the CMB, noise, and foreground covriance
try:
    ignore_fg = param_dict['ignore_fg']
except:
    ignore_fg = []

if isinstance(final_comp, str):
    if final_comp.lower() == 'y':
        ignore_fg.append('tsz') #the required component need not go into the covariance matrix.
    else:
        ignore_fg.append(final_comp.lower()) #the required component need not go into the covariance matrix.
print(ignore_fg)

#freqarr = [145]
#param_dict['which_gal_mask'] = 0
cl_dic = {}
for which_spec in which_spec_arr:
    if which_spec == 'TT':
        el, cl_dic[which_spec] = ilc.get_analytic_covariance(param_dict, freqarr, nl_dic = nl_dic['T'], ignore_fg = ignore_fg, include_gal = include_gal, bl_dic = bl_dic, cib_corr_coeffs = cib_corr_coeffs, use_websky_cib = use_websky_cib, use_sptspire_for_hfbands = use_sptspire_for_hfbands)
        ###legend(loc = 3, ncol = 4, fontsize = 8);show(); sys.exit()
    else:
        el, cl_dic[which_spec] = ilc.get_analytic_covariance (param_dict, freqarr, nl_dic = nl_dic['P'], ignore_fg = ignore_fg, which_spec = which_spec, pol_frac_per_cent_dust = param_dict['pol_frac_per_cent_dust'], pol_frac_per_cent_radio = param_dict['pol_frac_per_cent_radio'], pol_frac_per_cent_tsz = param_dict['pol_frac_per_cent_tsz'], pol_frac_per_cent_ksz = param_dict['pol_frac_per_cent_ksz'], include_gal = include_gal, bl_dic = bl_dic, use_websky_cib = use_websky_cib, use_sptspire_for_hfbands = use_sptspire_for_hfbands)

# In[13]:


if expname.find('spt4')>-1:
    colordic = {}
    colordic[90] = 'navy'
    colordic[150] = 'blue'
    colordic[220] = 'darkgreen'
    colordic[225] = 'goldenrod'
    colordic[286] = 'orangered'
    colordic[345] = 'maroon'
elif expname.find('herschel')>-1:
    colordic = {}
    colordic[90] = 'navy'
    colordic[150] = 'green'
    colordic[220] = 'darkred'
    colordic[225] = 'goldenrod'
    colordic[286] = 'orangered'
    colordic[345] = 'maroon'
    colordic[545] = 'darkorchid'
    colordic[600] = 'm'
    colordic[857] = 'lime'
    colordic[1200] = 'cyan'
else:
    colordic = {}
    colordic[90] = 'navy'
    colordic[150] = 'green'
    colordic[220] = 'darkred'


# In[14]:


if (0):
    pass
    '''
    #get the residual power now
    weights_dic, cl_residual = {}, {}
    for which_spec in which_spec_arr:
        print(which_spec)
        cl_residual[which_spec], weights_dic[which_spec] = ilc.residual_power(param_dict, freqarr, el, cl_dic, which_spec, final_comp = final_comp, freqcalib_fac = freqcalib_fac, return_weights = 1)

    print(weights_dic.keys())
    '''


# In[15]:


if null_comp is not None:
    cl_dic_temperature = {}
    cl_dic_temperature['TT'] = cl_dic['TT']
    cl_residual_comp_nulled_arr, weights_comp_nulled_arr = ilc.residual_power_new(param_dict, freqarr, el, cl_dic_temperature, final_comp = final_comp, freqcalib_fac = freqcalib_fac, return_weights = 1, null_comp = null_comp)

    #get weights and residual dic
    cl_residual_comp_nulled, weights_comp_nulled_dic = {}, {}
    for which_spec in ['TT']:
        if which_spec == 'TT':
           cl_res = cl_residual_comp_nulled_arr[0]
           weights = weights_comp_nulled_arr[:nc, 0]
        elif which_spec == 'EE':
           cl_res = cl_residual_comp_nulled_arr[1]
           weights = weights_comp_nulled_arr[nc:, 1]
        elif which_spec == 'TE':
           cl_res = cl_residual_comp_nulled_arr[2]
           weights = np.asarray( [weights_comp_nulled_arr[nc:, 0], weights_comp_nulled_arr[:nc, 1]] )
        cl_residual_comp_nulled[which_spec], weights_comp_nulled_dic[which_spec] = cl_res, weights
    print(cl_residual_comp_nulled.keys())


# In[16]:


cl_residual_arr, weights_arr = ilc.residual_power_new(param_dict, freqarr, el, cl_dic, final_comp = final_comp, freqcalib_fac = freqcalib_fac, return_weights = 1)#, null_comp = null_comp)

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


# In[17]:


#plot weights now
if (0):
    #clf()
    #fig = figure(figsize=(6, 6))
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

                if which_spec == 'TT' and null_comp is not None:
                    weights_comp_nulled_arr = weights_comp_nulled_dic[which_spec][frqcntr+shift]
                    plot(weights_comp_nulled_arr, color = colordic[freq], ls = '-.', lw = lwval)#, label = r'%s' %(freq))

            plot(np.sum(weights_dic[which_spec], axis = 0), 'k--', lw = 0.5, label = r'Sum')    
            if which_spec == 'TT' and null_comp is not None:
                plot(np.sum(weights_comp_nulled_dic[which_spec], axis = 0), 'k:', lw = 2., label = r'Sum')    
            if which_spec == 'TT':
                title(r'Temperature')
            else:
                title(r'Polarisation')
        axhline(lw=0.3);
        if cntr == len(which_spec)-1:
            xlabel(r'Multipole $\ell$');
        else:
            setp(ax.get_xticklabels(which = 'both'), visible=False)
        if cntr == 0:
            legend(loc = 4, ncol = 7, fontsize = 6)

        ylabel(r'Weight $W_{\ell}$')
        if null_comp is not None:
            ylim(-3., 3.);
        else:
            ylim(-2., 2.);
        xlim(100, 1e4)
    show()#;sys.exit()


# In[18]:


#plot and results file name
freqarr_str = '-'.join( np.asarray( freqarr ).astype(str) )
which_spec_arr_str = '-'.join( np.asarray( which_spec_arr ).astype(str) )
parent_folder = 'results/spt/20200708/'
if use_websky_cib:
    parent_folder = 'results/spt/20200708/websky_cib/'
elif use_545:
    parent_folder = 'results/spt/20200708/with_545/'
elif use_sptspire_for_hfbands:
    parent_folder = 'results/spt/20200708/spt_spire_cib_hf/'
if expname.find('spt4')>-1:
    parent_folder = '%s/spt4' %(parent_folder)
opfname = '%s/%s_ilc_%s_%s_%s.npy' %(parent_folder, expname, final_comp, freqarr_str, which_spec_arr_str)

if remove_atm:
    opfname = opfname.replace('.npy', '_noatmnoise.npy')
    
if null_comp is not None:
    if np.ndim(null_comp) == 0:
        null_comp_str = '%s_nulled' %(null_comp)        
    else:
        null_comp_str = '%s_nulled' %('-'.join(null_comp))
    null_comp_str = null_comp_str.replace('nulled_nulled', 'nulled')
        
    opfname = '%s_%s.npy' %(opfname.replace('.npy', ''), null_comp_str)

if split_cross:
    opfname = opfname.replace('.npy', '_splitcross.npy')

plname = opfname.replace('.npy', '.png').replace(parent_folder, '%s/plots/' %(parent_folder))
plfolder = '/'.join(plname.split('/')[:-1])
os.system('mkdir -p %s' %(plfolder))
print(opfname)
print(plname)


# In[19]:


freq0, lmax = param_dict['freq0'], param_dict['lmax']
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
fsval = 8
lwval = 0.75
plot_weights = 1
xmin, xmax = 20, 10000
xmin, xmax = 100, 10000
ymin, ymax = 1e-9, 100000.
if plot_weights:
    tr, tc = 6, len(which_spec_arr)
    subplots_adjust(wspace=0.1, hspace = 0.1)
    #first plot weights
    rspan, cspan = 2, 1
    curr_row = 0
    for cntr, which_spec in enumerate( which_spec_arr ):
        ax = subplot2grid((tr,tc), (curr_row, cntr), rowspan = rspan, colspan = cspan)#, xscale = 'log')#, yscale = 'log')
        if which_spec == 'TE':
            tot_teiter = 1 ##2
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
            legend(loc = 3, fontsize = 5, ncol = 4, handlelength = 2., handletextpad = 0.1)
        else:
            setp(ax.get_yticklabels(which = 'both'), visible=False)
        ylim(-1., 2.);
        xlim(xmin, xmax);
        for label in ax.get_xticklabels(): label.set_fontsize(fsval)
        for label in ax.get_yticklabels(): label.set_fontsize(fsval)        

        title(r'%s' %(which_spec))#, fontsize = 10)

    curr_row = rspan
    rspan = tr - rspan
for cntr, which_spec in enumerate( which_spec_arr ):
    if plot_weights:
        #ax = subplot(1,2,cntr+1, xscale = 'log', yscale = 'log')
        ax = subplot2grid((tr,tc), (curr_row, cntr), rowspan = rspan, colspan = cspan, yscale = 'log')#, xscale = 'log')
    else:
        ax = subplot(1,2,cntr+1, xscale = 'log', yscale = 'log')
    plot(el, cl_residual[which_spec], 'k', lw = 2., label = r'Residual')
    if which_spec == 'TT':
        if final_comp == 'cmb':
            plot(el_camb, cl_TT, 'gray', lw = 1., label = r'TT')
        if null_comp is not None:
            plot(el, cl_residual_comp_nulled[which_spec], 'black', lw = 2., ls = '-.', label = r'Residual: comp-nulled')
        
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
    for freq in freqarr:
        plot(el, cl_dic[which_spec][(freq,freq)], color = colordic[freq], lw = 0.5, ls = '-', label = r'%s' %(freq), alpha = 1.)        
    for freq in freqarr:
        if which_spec == 'TT':
            nl = nl_dic['T'][(freq, freq)]
        elif which_spec == 'EE':
            nl = nl_dic['P'][(freq, freq)]
        elif which_spec == 'TE':
            nl = nl_dic['T'][(freq, freq)] * 0.
        plot(el, nl, color = colordic[freq], lw = lwval, ls = '--', label = r'Noise: %s' %(freq))#, alpha = 0.5)
    #legend(loc=3, fancybox=1, ncol = 4, fontsize = 6);
    
    xlim(xmin, xmax);
    #ylim(1e-8,1e6);

    xlim(xmin, xmax);
    ylim(ymin, ymax); 

    xlabel(r'Multipole $\ell$')
    if cntr == 0: 
        ylabel(r'$C_{\ell}$')
        legend(loc = 1, fontsize = 6, ncol = 2, handlelength = 2., handletextpad = 0.1)
    else:
        setp(ax.get_yticklabels(which = 'both'), visible=False)
    for label in ax.get_xticklabels(): label.set_fontsize(fsval)
    for label in ax.get_yticklabels(): label.set_fontsize(fsval)
        
#tit = 'Galaxy = %s; Mask = %s; Bands = %s' %(include_gal, param_dict['which_gal_mask'], str(freqarr))
if remove_atm:
    tit = 'Bands = %s; no 1/f' %(str(freqarr))
else:
    tit = 'Bands = %s' %(str(freqarr))
suptitle(r'%s: %s' %(expname.replace('_','\_').upper(), final_comp))
#show()
savefig(plname)


# In[ ]:





# In[20]:


opdic = {}
opdic['el'] = el
opdic['cl_residual'] = cl_residual
opdic['cl_dic'] = cl_dic
opdic['freqcalib_fac'] = freqcalib_fac
opdic['param_dict'] = param_dict
opdic['weights'] = weights_dic
if null_comp is not None:
    opdic['null_comp'] = null_comp
    opdic['cl_residual_nulled'] = cl_residual_comp_nulled
    opdic['weights_nulled'] = weights_comp_nulled_dic
opdic['nl_dic'] = nl_dic
opdic['beam_noise_dic'] = beam_noise_dic
opdic['elknee_dic'] = elknee_dic
np.save(opfname, opdic)
print(opfname)
print('\n\n\n')


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




