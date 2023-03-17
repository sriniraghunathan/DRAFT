import numpy as np
#fname = 's4like_mask_v2/TT-EE/baseline/s4wide_ilc_galaxy1_27-39-93-145-225-278_TT-EE_galmask2_AZ_for7years.npy'
#fname = '20230317/s4wide/s4like_mask_v2/TT-EE/baseline/s4wide_ilc_galaxy1_27-39-93-145-225-278_TT-EE_galmask2_AZ_for7years.npy'
#fname = '20230317/s4wide/planck_mask/TT-EE/baseline/s4wide_ilc_galaxy1_27-39-93-145-225-278_TT-EE_galmask2_AZ_for7years.npy'
fname = '20230317/s4wide/planck_mask/TT-EE/baseline/s4wide_ilc_galaxy1_27-39-93-145-225-278_TT-EE_galmask1_AZ_for7years.npy'

res_dic = np.load(fname, allow_pickle = 1, encoding = 'latin1').item()
cl_residual = res_dic['cl_residual']
print('\nILC residuals: %s' %(cl_residual.keys()))
for which_spec in cl_residual:
    ilc_nl = cl_residual[which_spec] # ILC residuals
    print('\tSpectra = %s; length = %s' %(which_spec,len(ilc_nl) ))



fg_res_dic = res_dic['fg_res_dic']
print('\nNow residual foreground signals: %s' %(fg_res_dic.keys()))
for which_spec in fg_res_dic:
    print('\tSpectra = %s' %which_spec)
    for gal_res_name in fg_res_dic[which_spec]:
        gal_res_cl = fg_res_dic[which_spec][gal_res_name]
        print('\t\t Signal = %s; length = %s' %(gal_res_name, len(gal_res_cl)))
