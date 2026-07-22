import argparse, sys, numpy as np, scipy as sc, warnings, os
sys.path.append('modules')
import misc, exp_specs
import ilc, foregrounds as fg
import pickle, gzip

#import matplotlib.cbook
warnings.filterwarnings('ignore',category=RuntimeWarning)

parser = argparse.ArgumentParser(description='')
parser.add_argument('-expname', dest='expname', action='store', help='expname', type=str, required=True)
parser.add_argument('-total_obs_time', dest='total_obs_time', action='store', help='total_obs_time in years', type=float, default=7.0)
parser.add_argument('-include_gal', dest='include_gal', action='store', help='include_gal', type=int, default=0)
parser.add_argument('-which_gal_mask', dest='which_gal_mask', action='store', help='which_gal_mask', type=int, default=2)
parser.add_argument('-interactive_mode', dest='interactive_mode', action='store', help='interactive_mode', type=int, default=1)
parser.add_argument('-save_fg_res_and_weights', dest='save_fg_res_and_weights', action='store', help='save_fg_res_and_weights', type=int, default=1)
parser.add_argument('-s4_so_joint_configs', dest='s4_so_joint_configs', action='store', help='s4_so_joint_configs', type=int, default=0)
parser.add_argument('-include_fulls4scaledsobaseline', dest='include_fulls4scaledsobaseline', action='store', help='include_fulls4scaledsobaseline', type=int, default=0)

#20230530 - scale noise levels of bands
parser.add_argument('-noise_scalings_for_bands', dest='noise_scalings_for_bands', action='store', help='noise_scalings_for_bands', nargs='+', type=float, default=None)

#20230531 - option to get CMB or y
parser.add_argument('-final_comp', dest='final_comp', action='store', help='final_comp', type=str, default='cmb')
parser.add_argument('-debug', dest='debug', action='store', help='debug', type=int, default=0)


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

# In[24]:


if not debug:
    import matplotlib
    ##matplotlib.use('Agg')
from pylab import *

rcParams['figure.dpi'] = 150
rcParams['font.family'] = 'serif'
rcParams["figure.facecolor"] = 'white'

#some constants
h=6.62607004e-34 #Planck constant in m2 kg / s
k_B=1.38064852e-23 #Boltzmann constant in m2 kg s-2 / K-1
Tcmb = 2.73 #Kelvin

###total_obs_time = float(sys.argv[1]) #in years
# In[25]:


#params
paramfile = 'params.ini'

# read and store param dict
param_dict = misc.get_param_dict(paramfile)

if not os.path.exists(param_dict['data_folder']):
    param_dict['data_folder'] = '/data/spt/sri-data48/git/DRAFT/data/'

el = np.arange(param_dict['lmax'])

#20220112 - moved to argparse
###include_gal = param_dict['include_gal'] ##1  
param_dict['include_gal'] = include_gal
if not include_gal:
    which_gal_mask = -1
param_dict['which_gal_mask'] = which_gal_mask

try:
    remove_atm = param_dict['remove_atm']
except:
    remove_atm = 0

# In[26]:


#S4 specs
if expname == 's4deepv3r025_plus_s4wide':
    specs_dic, corr_noise_bands, rho, corr_noise, Nred_dic = exp_specs.get_exp_specs('s4deepv3r025', remove_atm = remove_atm)
    specs_dic_s4wide, corr_noise_bands_s4wide, rho_s4wide, corr_noise_s4wide, Nred_dic = exp_specs.get_exp_specs('s4wide')
#20220225 - 2028 ASO + single-CHLAT - we will take an inverse variance combination of Nyear single CHLAT and N+1 ASO
elif expname.find('s4wide_single_chlat_plus_2028aso')>-1: 
    specs_dic_aso, corr_noise_bands_aso, rho_aso, corr_noise_aso, Nred_dic = exp_specs.get_exp_specs('s4wide_scaled_aso')
    specs_dic, corr_noise_bands, rho, corr_noise, Nred_dic = exp_specs.get_exp_specs('s4wide_single_chlat')
else:
    specs_dic, corr_noise_bands, rho, corr_noise, Nred_dic = exp_specs.get_exp_specs(expname, remove_atm = remove_atm)

#20220223 - include full SO-baseline, if requested
if include_fulls4scaledsobaseline:
    specs_dic_fulls4scaledsobaseline, corr_noise_bands_fulls4scaledsobaseline, rho_fulls4scaledsobaseline, corr_noise_fulls4scaledsobaseline = exp_specs.get_exp_specs('s4wide_scaled_sobaseline')

freqarr = sorted( specs_dic.keys() )
nc = len( freqarr )
freqcalib_fac = None
null_comp = None
TParr = ['T', 'P']
which_spec_arr = ['TT', 'EE']
reduce_cib_power = None
total_obs_time_default = 7. ###10. #years
if expname.find('s4_all_chile_config_lat_')>-1: #20250504
    total_obs_time_default = 10.
if expname.find('cmbhd')>-1:
    reduce_cib_power = 17. #150 GHz power reduction after removing sources above 0.04 mJy

#cl multipler - multiply a given spectra by some amount to perform partial ILC. similar to https://arxiv.org/abs/2102.05033
cl_multiplier_dic = {}

if (1): #20230530
    freqarr_str = '-'.join( np.asarray( freqarr ).astype(str) )
    which_spec_arr_str = '-'.join( np.asarray( which_spec_arr ).astype(str) )
    parent_folder = 'results/20210506_with202102designtoolinputforpySM3sims_sedscalingfordust'
    parent_folder = '%s/202310xx_modified_PBDR_config_for_Neff_paper/' %(parent_folder) #20231025 - modified PBDR config: https://docs.google.com/spreadsheets/d/10fL76XTzhgP_B_GKsEW4nqNTkRgvp2dh4zYh6Y-G2AE/edit#gid=0

    if expname.find('s4_all_chile_config_lat_')>-1 or expname.find('advanced_so')>-1:
        #parent_folder = 'results/s4_all_chile_config'
        parent_folder = 'results/s4_all_chile_config/report/'

    if final_comp != 'cmb':
        parent_folder = '%s/%s/' %(parent_folder, final_comp)

    if noise_scalings_for_bands is not None and len(np.unique(noise_scalings_for_bands))>1: #20230530 - scale noise levels of bands
        parent_folder = '%s/noise_scalings/' %(parent_folder)
        noise_scalings_for_bands_str = '-'.join([str(n) for n in noise_scalings_for_bands])
        noise_scalings_for_bands_str = '_noisescalings%s' %(noise_scalings_for_bands_str)

    if s4_so_joint_configs:
        parent_folder = '%s/s4_so_joint_configs/' %(parent_folder)
    if null_comp is not None:
        null_comp_str = 'nulled_%s' %('-'.join(null_comp))
        parent_folder = '%s/%s/' %(parent_folder, null_comp_str)

    if param_dict['lmax']>5002:
        parent_folder = '%s/lmax_%s/' %(parent_folder, param_dict['lmax'])    
    parent_folder = '%s/%s/' %(parent_folder, expname)

    opfname = '%s/%s_ilc_galaxy%s_%s_%s.npy' %(parent_folder, expname, include_gal, freqarr_str, which_spec_arr_str)
    if null_comp is not None:
        opfname = '%s_%s.npy' %(opfname.replace('.npy', ''), null_comp_str)

    if not corr_noise:
        opfname = opfname.replace('.npy', '_nocorrnoise.npy')

    if expname.find('s4')>-1 or expname.find('cmbhd')>-1:
        opfname = opfname.replace(parent_folder, '%s/planck_mask/TT-EE/baseline/' %(parent_folder))

    if include_gal:
        opfname = opfname.replace('.npy', '_galmask%s.npy' %(which_gal_mask))

    if remove_atm:
        opfname = opfname.replace('.npy', '_noatmnoise.npy')

    if cl_multiplier_dic is not None:
        if len(cl_multiplier_dic) > 1:
            cl_multiplier_str = 'pilcmultfacs'
            for kkk in cl_multiplier_dic:
                cl_multiplier_str = '%s-%s%s'%(cl_multiplier_str, kkk, cl_multiplier_dic[kkk])
            cl_multiplier_str = cl_multiplier_str.strip('-')
            opfname = opfname.replace('.npy', '_%s.npy' %(cl_multiplier_str))

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

    if param_dict['lmax']!=5000:
        opfname = opfname.replace('.npy', '_lmax%s.npy' %(param_dict['lmax']))

    if noise_scalings_for_bands is not None and len(np.unique(noise_scalings_for_bands))>1: #20230530 - scale noise levels of bands
        opfname = opfname.replace('.npy', '%s.npy' %(noise_scalings_for_bands_str))


    opfolder = '/'.join(opfname.split('/')[:-1])
    if not os.path.exists(opfolder): os.system('mkdir -p %s' %(opfolder))

    if expname.lower().find('s4')>-1:#total_obs_time_default != total_obs_time:
        opfname = opfname.replace('.npy', '_for%gyears.npy' %(total_obs_time))
        
    print(opfname)

#beam and noise arr
beamarr = []
noisearr_T, elkneearr_T, alphakneearr_T = [], [], []
noisearr_P, elkneearr_P, alphakneearr_P = [], [], []
for fcntr, freq in enumerate( freqarr ):
    beam_arcmins, white_noise_T, elknee_T, alphaknee_T, white_noise_P, elknee_P, alphaknee_P = specs_dic[freq]

    if (1): #20230530 - noise scaling of bands
        if noise_scalings_for_bands is not None:
            assert len(noise_scalings_for_bands) == len(freqarr)
            white_noise_T = white_noise_T * noise_scalings_for_bands[fcntr]
            white_noise_P = white_noise_P * noise_scalings_for_bands[fcntr]

    if (1): #noise scaling based on total_obs_time
        ###print(total_obs_time, total_obs_time_default)
        noise_scaling_fac = (total_obs_time_default / total_obs_time)**0.5
        ###print(noise_scaling_fac); sys.exit()
        white_noise_T = white_noise_T * noise_scaling_fac
        white_noise_P = white_noise_P * noise_scaling_fac

    #add N+1 year aso via inverse variance combination now
    if expname.find('s4wide_single_chlat_plus_2028aso')>-1:
        #the above noise numbers for N year single CHLAT
        #now we will add N+1 year ASO
        total_obs_time_2028_aso = total_obs_time + 1
        white_noise_T_aso, white_noise_P_aso = specs_dic_aso[freq][1], specs_dic_aso[freq][4]
        noise_scaling_fac_2028_aso = (total_obs_time_default / (total_obs_time_2028_aso))**0.5
        white_noise_T_2028_aso = white_noise_T_aso * noise_scaling_fac_2028_aso
        white_noise_P_2028_aso = white_noise_P_aso * noise_scaling_fac_2028_aso

        white_noise_T = ( (1./white_noise_T**2.) + (1./white_noise_T_2028_aso**2.) )**-0.5
        white_noise_P = ( (1./white_noise_P**2.) + (1./white_noise_P_2028_aso**2.) )**-0.5


    #20220223 - include full SO-baseline, if requested
    if include_fulls4scaledsobaseline:
        beam_arcmins_2, white_noise_T_2, elknee_T_2, alphaknee_T_2, white_noise_P_2, elknee_P_2, alphaknee_P_2 = specs_dic_fulls4scaledsobaseline[freq]
        white_noise_T_2 = white_noise_T_2 * np.sqrt(7./5.)
        white_noise_P_2 = white_noise_P_2 * np.sqrt(7./5.)
        white_noise_T = (1./white_noise_T**2. + 1./white_noise_T_2**2.)**-0.5
        white_noise_P = (1./white_noise_P**2. + 1./white_noise_P_2**2.)**-0.5

    beamarr.append(beam_arcmins)
    noisearr_T.append(white_noise_T)
    noisearr_P.append(white_noise_P)
    elkneearr_T.append(elknee_T)
    elkneearr_P.append(elknee_P)
    alphakneearr_T.append(alphaknee_T)
    alphakneearr_P.append(alphaknee_P)    

print('\n')
#print(elkneearr_T)
print('Delta T =', noisearr_T)
print('Delta P =', noisearr_P)
#print(beamarr)
print('\n')
###sys.exit()

# In[28]:


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


# In[29]:


#get beam deconvolved noise nls
nl_dic = {}
for TPcntr, TP in enumerate( TParr ):
    nl_dic[TP]={}
    for freq1 in freqarr:
        beamval1, noiseval1 = beam_noise_dic[TP][freq1]
        elknee1, alphaknee1 = elknee_dic[TP][freq1]
        for freq2 in freqarr:        
            beamval2, noiseval2 = beam_noise_dic[TP][freq2]
            elknee2, alphaknee2 = elknee_dic[TP][freq2]

            ##elknee1, elknee2 = -1, -1.
            ##alphaknee1, alphaknee2 = -1., -1.

            Nred1, Nred2 = -1., -1.
            if freq1 in Nred_dic:
                Nred1 = Nred_dic[freq1][TPcntr]
            if freq2 in Nred_dic:
                Nred2 = Nred_dic[freq2][TPcntr]
            
            if freq1 == freq2:
                nl = misc.get_nl(noiseval1, el, beamval1, elknee = elknee1, alphaknee = alphaknee1, Nred1 = Nred1)
            else:
                if freq2 in corr_noise_bands[freq1]:
                    nl = misc.get_nl(noiseval1, el, beamval1, elknee = elknee1, alphaknee = alphaknee1, beamval2 = beamval2, noiseval2 = noiseval2, elknee2 = elknee2, alphaknee2 = alphaknee2, rho = rho, Nred1 = Nred1, Nred2 = Nred2)
                else:
                    nl = np.zeros( len(el) )
            nl[el<=param_dict['lmin']] = 0.
            ##nl[nl == 0.] = np.min(nl[nl!=0.])/1e3

            if expname == 's4deepv3r025_plus_s4wide' and np.sum(nl)>0.: #20210428: inv. var. combination of S4-Wide and S4-Ultra deep nl
                beam_arcmins_s4wide_f1, white_noise_T_s4wide_f1, elknee_T_s4wide_f1, alphaknee_T_s4wide_f1, white_noise_P_s4wide_f1, elknee_P_s4wide_f1, alphaknee_P_s4wide_f1 = np.copy(specs_dic_s4wide[freq1])
                beam_arcmins_s4wide_f2, white_noise_T_s4wide_f2, elknee_T_s4wide_f2, alphaknee_T_s4wide_f2, white_noise_P_s4wide_f2, elknee_P_s4wide_f2, alphaknee_P_s4wide_f2 = np.copy(specs_dic_s4wide[freq2])
                if TP == 'T':
                    white_noise_s4wide_f1, elknee_s4wide_f1, alphaknee_s4wide_f1 = white_noise_T_s4wide_f1, elknee_T_s4wide_f1, alphaknee_T_s4wide_f1
                    white_noise_s4wide_f2, elknee_s4wide_f2, alphaknee_s4wide_f2 = white_noise_T_s4wide_f2, elknee_T_s4wide_f2, alphaknee_T_s4wide_f2
                elif TP == 'P':
                    white_noise_s4wide_f1, elknee_s4wide_f1, alphaknee_s4wide_f1 = white_noise_P_s4wide_f1, elknee_P_s4wide_f1, alphaknee_P_s4wide_f1
                    white_noise_s4wide_f2, elknee_s4wide_f2, alphaknee_s4wide_f2 = white_noise_P_s4wide_f2, elknee_P_s4wide_f2, alphaknee_P_s4wide_f2

                if (0): #noise scaling based on total_obs_time
                    noise_scaling_fac = (total_obs_time_default / total_obs_time)**0.5
                    white_noise_s4wide_f1 = white_noise_s4wide_f1 * noise_scaling_fac
                    white_noise_s4wide_f2 = white_noise_s4wide_f2 * noise_scaling_fac        

                #print('s4wide', white_noise_s4wide_f1, elknee_s4wide_f1, alphaknee_s4wide_f1, white_noise_s4wide_f2, elknee_s4wide_f2, alphaknee_s4wide_f2)
                #print('s4deep', noiseval1, elknee1, alphaknee1, noiseval2, elknee2, alphaknee2)
                if freq1 == freq2:
                    nl_s4_wide = misc.get_nl(white_noise_s4wide_f1, el, beam_arcmins_s4wide_f1, elknee = elknee_s4wide_f1, alphaknee = alphaknee_s4wide_f1)
                else:
                    if freq2 in corr_noise_bands_s4wide[freq1]:
                        nl_s4_wide = misc.get_nl(white_noise_s4wide_f1, el, beam_arcmins_s4wide_f1, elknee = elknee_s4wide_f1, alphaknee = alphaknee_s4wide_f1, beamval2 = beam_arcmins_s4wide_f2, noiseval2 = white_noise_s4wide_f2, elknee2 = elknee_s4wide_f2, alphaknee2 = alphaknee_s4wide_f2, rho = rho_s4wide)
                    else:
                        nl_s4_wide = np.zeros( len(el) )
                nl_s4_wide[el<=param_dict['lmin']] = 0.

                #perform inverse variance combination now
                nl_s4_ultradeep = np.copy(nl)
                nl = 1./ ( (1./nl) + (1./nl_s4_wide) )

                if (0):#freq1 == 93:# and freq2 == 145:
                    loglog(el, nl_s4_ultradeep, color = 'black', label = 'S4-Ultra deep'); loglog(el, nl_s4_wide, color = 'red', label = 'S4-Wide'); 
                    loglog(el, nl, color = 'darkgreen', label = 'S4'); 
                    title('%s: (%s, %s)' %(TP, freq1, freq2)); legend(loc = 3); show()

            nl_dic[TP][(freq1, freq2)] = nl

print(nl_dic['T'].keys())
##sys.exit()

# In[30]:


#get beams
bl_dic = misc.get_beam_dic(freqarr, beam_noise_dic['T'], param_dict['lmax'])
print(bl_dic.keys())
if (0):
    for freq in freqarr:
        plot(bl_dic[freq], label = freq)
    legend(loc = 1)

if (0):
    use_dls = True
    beam_decon = True
    color_arr = ['navy', 'blue', 'darkgreen', 'goldenrod', 'orangered', 'darkred']
    ax=subplot(111, yscale = 'log')
    for fcntr, freq in enumerate( freqarr ):
        currnl = nl_dic['T'][(freq,freq)]*bl_dic[freq]**2.
        #currnl = nl_dic['P'][(freq,freq)]*bl_dic[freq]**2.
        tmpinds = np.where( (el>=3000) & (el<=5000) )[0]
        meannl = np.median(currnl[tmpinds])
        noise_uk_arcmin = np.sqrt(meannl)/np.radians(1./60.)
        print(noise_uk_arcmin, beam_noise_dic['T'][freq])

        if beam_decon:
            currnl = currnl / bl_dic[freq]**2.
        
        if use_dls:
            dl_fac = el * (el+1)/2/np.pi
            plot(el, dl_fac * currnl, label = '%s GHz ' %(freq), ls = '-', color = color_arr[fcntr])
        else:
            plot(el, currnl, label = '%s GHz (%.2f $\mu$K-arcmin)' %(freq, noise_uk_arcmin), ls = '-', color = color_arr[fcntr])
        #plot(nl_dic_actual['T'][(freq,freq)], lw = 2., color = colordic[freq])
    legend(loc = 1)
    if not use_dls:
        xlim(0, 5000); ylim(1e-8, 1.)
        ylabel(r'N$_{\ell}$ [$\mu$K$^{2}$]', fontsize = 14)
    else:
        ylabel(r'$\ell(\ell+1)/(2\pi)$ N$_{\ell}$ [$\mu$K$^{2}$]', fontsize = 14)
        xlim(0, 5000); ylim(.1, 1e5)
    xlabel(r'Multipole $\ell$', fontsize = 14)
    expname_str = expname.replace('spt3g_', 'SPT-3G: ').replace('summer', 'Summer').replace('_', '\_')
    ##expname_str = 'S4-Wide'
    title(r'%s' %(expname_str), fontsize = 14)
    ##savefig('s4_wide_nl.png', dpi= 200.); sys.exit()
    show(); sys.exit()

# In[31]:


#get the CMB, noise, and foreground covriance
try:
    ignore_fg = param_dict['ignore_fg']
except:
    ignore_fg = []

ignore_fg.append(final_comp.lower()) #the required component need not go into the covariance matrix.
ignore_fg.append('tsz_cib')
print(ignore_fg)

#freqarr = [145]
#param_dict['which_gal_mask'] = 0
cl_dic = {}
fg_cl_dic = {}
for which_spec in which_spec_arr:
    if which_spec == 'TT':
        el, cl_dic[which_spec], fg_cl_dic[which_spec] = ilc.get_analytic_covariance(param_dict, freqarr, el = el,nl_dic = nl_dic['T'], ignore_fg = ignore_fg, include_gal = include_gal, bl_dic = bl_dic, reduce_cib_power = reduce_cib_power, cl_multiplier_dic = cl_multiplier_dic)
    else:
        el, cl_dic[which_spec], fg_cl_dic[which_spec] = ilc.get_analytic_covariance(param_dict, freqarr, el = el, nl_dic = nl_dic['P'], ignore_fg = ignore_fg, which_spec = which_spec,                     pol_frac_per_cent_dust = param_dict['pol_frac_per_cent_dust'],                     pol_frac_per_cent_radio = param_dict['pol_frac_per_cent_radio'],                     pol_frac_per_cent_tsz = param_dict['pol_frac_per_cent_tsz'],                     pol_frac_per_cent_ksz = param_dict['pol_frac_per_cent_ksz'],                     include_gal = include_gal, bl_dic = bl_dic, reduce_cib_power = reduce_cib_power)
    ###print(which_spec, cl_dic[which_spec])
print(cl_dic.keys())
print(el)
###sys.exit()
op_dic_for_fg_noise = {}
op_dic_for_fg_noise['fg_cl_dic'] = fg_cl_dic
op_dic_for_fg_noise['nl_dic'] = nl_dic
##np.save('data/cmbs4_nl_and_fgcl_dict.npy', op_dic_for_fg_noise); ##sys.exit()

#get the residual power now
#null_comp = None
cl_residual, weights_dic = {}, {}
if null_comp is None:
    cl_residual_arr, weights_arr = ilc.residual_power(param_dict, freqarr, el, cl_dic, final_comp = final_comp, freqcalib_fac = freqcalib_fac, return_weights = 1)
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
else:
    cl_residual, weights_dic = {}, {}
    cl_dic_TT, cl_dic_EE = {}, {}
    cl_dic_TT['TT'] = cl_dic['TT']
    cl_dic_EE['EE'] = cl_dic['EE']
    cl_residual_TT_arr, weights_TT_arr = ilc.residual_power(param_dict, freqarr, el, cl_dic_TT, final_comp = final_comp, freqcalib_fac = freqcalib_fac, return_weights = 1, null_comp = null_comp)
    cl_residual_EE_arr, weights_EE_arr = ilc.residual_power(param_dict, freqarr, el, cl_dic_EE, final_comp = final_comp, freqcalib_fac = freqcalib_fac, return_weights = 1, null_comp = null_comp)
    cl_residual['TT'], weights_dic['TT'] = cl_residual_TT_arr[0], weights_TT_arr[:nc, 0]
    cl_residual['EE'], weights_dic['EE'] = cl_residual_EE_arr[0], weights_EE_arr[:nc, 0]

print(cl_residual.keys())
# In[38]:

fg_res_dic = {}
#signal_arr = ['galdust', 'galsync']
if include_gal:
    signal_arr = ['tsz', 'cib', 'radio', 'galdust', 'galsync', 'noise']#, 'tsz-cib']
else:
    signal_arr = ['tsz', 'cib', 'radio', 'noise']#, 'tsz-cib']
for which_spec in ['TT', 'EE']:
    fg_res_dic[which_spec] = {}
    for elcnt, currel in enumerate(el):
        if (elcnt%2500) == 0: print(which_spec, elcnt)
        for s in signal_arr:
            #print(s)
            if s == 'galdust':
                #curr_cl_dic = cl_gal_dust_dic[which_spec]
                curr_cl_dic = fg_cl_dic[which_spec]['galdust']
            elif s == 'galsync':
                #curr_cl_dic = cl_gal_sync_dic[which_spec]
                curr_cl_dic = fg_cl_dic[which_spec]['galsync']
            elif s == 'ksz':
                #curr_cl_dic = cl_ksz_dic[which_spec]
                curr_cl_dic = fg_cl_dic[which_spec]['ksz']
            elif s == 'cmb':
                #curr_cl_dic = cl_cmb_dic[which_spec]
                curr_cl_dic = fg_cl_dic[which_spec]['cmb']
            elif s == 'tsz':
                #curr_cl_dic = cl_tsz_dic[which_spec]
                curr_cl_dic = fg_cl_dic[which_spec]['tsz']
            elif s == 'cib':
                #curr_cl_dic = cl_dust_dic[which_spec]
                curr_cl_dic = fg_cl_dic[which_spec]['cib']
            elif s == 'radio':
                #curr_cl_dic = cl_radio_dic[which_spec]
                curr_cl_dic = fg_cl_dic[which_spec]['radio']
            elif s == 'tsz-cib':
                #curr_cl_dic = cl_tsz_cib_dic[which_spec]
                curr_cl_dic = fg_cl_dic[which_spec]['tsz_cib']
            elif s == 'noise':
                curr_cl_dic = fg_cl_dic[which_spec]['noise']

            tmp_cl_dic = {which_spec: curr_cl_dic}
            clmat = np.mat( ilc.create_clmat(freqarr, elcnt, tmp_cl_dic) )
            currw_ilc = np.mat( weights_dic[which_spec][:, elcnt] )
            
            curr_res_ilc = np.asarray(np.dot(currw_ilc, np.dot(clmat, currw_ilc.T)))[0][0]
            if s not in fg_res_dic[which_spec]:
                fg_res_dic[which_spec][s] = []
            fg_res_dic[which_spec][s].append( curr_res_ilc )
    
    for s in signal_arr:
        fg_res_dic[which_spec][s] = np.asarray(fg_res_dic[which_spec][s])
print(fg_res_dic.keys())


freq0, lmax = param_dict['freq0'], param_dict['lmax']

#save residual files
if noise_scalings_for_bands is not None:
    save_fg_res_and_weights = 0
if include_gal:
    cl_gal_dic_dust_fname = param_dict['cl_gal_dic_dust_fname']
    try:
        cl_gal_folder = param_dict['cl_gal_folder']
        cl_gal_dic_dust_fname = '%s/%s' %(cl_gal_folder, cl_gal_dic_dust_fname)
    except:
        pass
    galdustsims_cl = np.load(cl_gal_dic_dust_fname, allow_pickle=1, encoding = 'latin1').item()
    if not include_gal:
        fsky_val = 0.68
    else:
        fsky_val = galdustsims_cl['fsky_arr'][param_dict['which_gal_mask']]
opdic = {}
opdic['el'] = el
if cl_multiplier_dic is not None:
    opdic['cl_multiplier_dic'] = cl_multiplier_dic
opdic['cl_residual'] = cl_residual
if save_fg_res_and_weights:
    opdic['fg_res_dic'] = fg_res_dic
    opdic['weights'] = weights_dic
opdic['freqcalib_fac'] = freqcalib_fac
opdic['param_dict'] = param_dict
if include_gal:
    opdic['fsky_val'] = fsky_val
opdic['which_gal_mask'] = which_gal_mask
#opdic['nl_dic'] = nl_dic
opdic['beam_noise_dic'] = beam_noise_dic
opdic['elknee_dic'] = elknee_dic
np.save(opfname, opdic)
print(opfname)

print('\nDone.\n')



