import pickle, gzip, numpy as np, glob, sys, os
#fname = 's4deepv3r025_202310xx_pbdr_config/s4deepv3r025_202310xx_pbdr_config_ilc_galaxy0_20-27-39-93-145-225-278_TT-EE_lmax6500_for7years.npy'
fname = 's4wide_202310xx_pbdr_config/s4wide_202310xx_pbdr_config_ilc_galaxy1_27-39-93-145-225-278_TT-EE_galmask2_AZ_lmax6500_for7years.npy'

res_dic = np.load(fname, allow_pickle = 1, encoding = 'latin1').item()
print(fname, res_dic.keys())
cl_residual_dic = res_dic['cl_residual'] #ILC residuals
els = res_dic['el'] #ells

'''
#push the noise curves to an ascii file.
ascii_fname = fname.replace('.npy', '.txt')
oparr = np.asarray( [els, cl_residual_dic['TT'], cl_residual_dic['EE']] ).T
np.savetxt(ascii_fname, oparr, header = 'els cl_ilc_tt cl_ilc_ee')
print(ascii_fname)
'''
