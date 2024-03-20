import pickle, gzip, numpy as np, glob, sys, os
fname = 's4wide_202310xx_pbdr_config/s4wide_202310xx_pbdr_config_ilc_galaxy1_27-39-93-145-225-278_TT-EE_galmask2_AZ_lmax6500_for7years.npy'
res_dic = np.load(fname, allow_pickle = 1, encoding = 'latin1').item()
print(fname, res_dic.keys())
cl_residual_dic = res_dic['cl_residual'] #ILC residuals
els = res_dic['el'] #ells

