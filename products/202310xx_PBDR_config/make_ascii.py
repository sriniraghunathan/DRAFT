import pickle, gzip, numpy as np, glob, sys, os

#fd = 's4wide/'
#fd = 's4deepv3r025_202310xx_pbdr_config/'
#fd = 's4wide_202310xx_pbdr_config/'
#fd = 's4wide_acheived_performance_pbdr_202312xx/'

#fd = 'lmax_12000/s4wide/'
fd = 'lmax_12000/s4wide_202310xx_pbdr_config/'
#fd = 'lmax_12000/s4wide_acheived_performance_pbdr_202312xx/'

flist = glob.glob( '%s/*.npy' %(fd) )
for fname in flist:
    res_dic = np.load(fname, allow_pickle = 1, encoding = 'latin1').item()
    ##print(fname, res_dic.keys())
    cl_residual_dic = res_dic['cl_residual'] #ILC residuals
    els = res_dic['el'] #ells

    #push the noise curves to an ascii file.
    ascii_fname = fname.replace('.npy', '.txt')
    oparr = np.asarray( [els, cl_residual_dic['TT'], cl_residual_dic['EE']] ).T
    np.savetxt(ascii_fname, oparr, header = 'els cl_ilc_tt cl_ilc_ee')
    print(ascii_fname)
