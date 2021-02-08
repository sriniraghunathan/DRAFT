import numpy as np
fname = 's4wide_ilc_galaxy1_27-39-93-145-225-278_TT-EE_galmask2_AZ_galresiduals.npy'
fg_res_dic = np.load(fname, allow_pickle = 1, encoding = 'latin1').item()
print('\n')
print(fg_res_dic.keys())
for which_spec in fg_res_dic:
    print('\tSpectra = %s' %which_spec)
    for gal_res_name in fg_res_dic[which_spec]:
        gal_res_cl = fg_res_dic[which_spec][gal_res_name]
        print('\t\t Signal = %s; length = %s' %(gal_res_name, len(gal_res_cl)))
