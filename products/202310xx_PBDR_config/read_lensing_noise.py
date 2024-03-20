import pickle, gzip, numpy as np, glob, sys, os
fname = 's4wide_202310xx_pbdr_config/lensing_noise_curves//s4wide_202310xx_pbdr_config_lmin100_lmax5000_lmaxtt3000.npy'
dic = np.load(fname, allow_pickle = 1, encoding='latin1').item()
#print(dic.keys()); sys.exit()
els, cl_kk, nl_tt, nl_eb, nl_mv, nl_mvpol = dic['els'], dic['cl_kk'], dic['Nl_TT'], dic['Nl_EB'], dic['Nl_MV'], dic['Nl_MVpol']
nl_te, nl_tb, nl_ee = dic['Nl_ET'], dic['Nl_TB'], dic['Nl_EE']
