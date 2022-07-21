import numpy as np 
fname = 's4wide_ilc_galaxy1_27-39-93-145-225-278_TT-EE_galmask2_AZ_lmin100_lmax5000_lmaxtt3000.npy'
lensing_noise_dic = np.load(fname, allow_pickle=True, encoding='latin1').item()

print(lensing_noise_dic.keys())

which_lensing_estimator = 'Nl_MV' #'Nl_MVpol' #'Nl_EB' #'Nl_TT'
els=lensing_noise_dic['els']
nl=lensing_noise_dic[which_lensing_estimator]

print(els, nl)