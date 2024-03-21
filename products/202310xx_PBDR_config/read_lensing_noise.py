import pickle, gzip, numpy as np, glob, sys, os
fname = 's4deepv3r025_202310xx_pbdr_config/lensing_noise_curves/s4deepv3r025_202310xx_pbdr_config_lmin100_lmax5000_lmaxtt3000.npy'
#fname = 's4wide_202310xx_pbdr_config/lensing_noise_curves//s4wide_202310xx_pbdr_config_lmin100_lmax5000_lmaxtt3000.npy'
#fname = 's4wide_acheived_performance_pbdr_202312xx/lensing_noise_curves/s4wide_acheived_performance_pbdr_202312xx_lmin100_lmax5000_lmaxtt3000.npy'
dic = np.load(fname, allow_pickle = 1, encoding='latin1').item()
#print(dic.keys()); sys.exit()
els, cl_kk, nl_tt, nl_eb, nl_mv, nl_mvpol = dic['els'], dic['cl_kk'], dic['Nl_TT'], dic['Nl_EB'], dic['Nl_MV'], dic['Nl_MVpol']
nl_te, nl_tb, nl_ee = dic['Nl_ET'], dic['Nl_TB'], dic['Nl_EE']


'''
#push the noise curves to an ascii file.
ascii_fname = fname.replace('.npy', '.txt')
oparr = np.asarray( [els, cl_kk.real, nl_tt.real, nl_ee.real, nl_te.real, nl_tb.real, nl_eb.real, nl_mv.real, nl_mvpol.real] ).T
np.savetxt(ascii_fname, oparr, header = 'els cl_kk n0_tt n0_ee n0_te n0_tb n0_eb n0_mv n0_mvpol')
print(ascii_fname)
'''
