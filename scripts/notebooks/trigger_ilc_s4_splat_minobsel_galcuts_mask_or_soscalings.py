import numpy as np, sys, os

pgmname = 'get_ilc_residuals_withcorrnoise_withgal_TTEETE_partialilcusingmultipliers_cilc_0.5.py'

if (0):
    include_gal = 1
    which_gal_mask_arr = [0, 1, 2, 3, 4, 5]
    which_gal_mask_arr = [3]
    '''
    total_obs_time_arr_1 = [0.1, 0.5]
    total_obs_time_arr_2 = np.arange(1.0, 20., 1.0)
    total_obs_time_arr = np.concatenate((total_obs_time_arr_1, total_obs_time_arr_2))
    '''
    total_obs_time_arr_1 = np.arange(1.0, 10.1, 1.0)
    total_obs_time_arr_2 = np.arange(12.0, 20.1, 2.0)
    total_obs_time_arr = np.concatenate((total_obs_time_arr_1, total_obs_time_arr_2))

    expname = 's4deepv3r025'
    interactive_mode = 0
    save_fg_res_and_weights = 0

    for which_gal_mask in which_gal_mask_arr:
        for total_obs_time in total_obs_time_arr:
            cmd = 'python3 %s -expname %s -include_gal %s -which_gal_mask %s -total_obs_time %s -interactive_mode %s -save_fg_res_and_weights %s' %(pgmname, expname, include_gal, which_gal_mask, total_obs_time, interactive_mode, save_fg_res_and_weights)
            print('\n###############\n%s\n' %(cmd))
            #os.system(cmd)
            sys.exit()

if (1): #SO scalings
    #expname_arr = ['s4wide', 's4wide_scaled_sobaseline', 's4wide_scaled_aso', 's4wide_single_chlat_plus_aso']
    #expname_arr = ['s4wide_scaled_sobaseline', 's4wide_scaled_aso', 's4wide_single_chlat_plus_aso']
    #expname_arr = ['s4wide']

    #expname_arr = ['s4wide', 's4wide_scaled_aso_plus_fulls4scaledsobaseline', 's4wide_single_chlat_plus_fulls4scaledsobaseline', 's4wide_single_chlat_plus_aso_plus_fulls4scaledsobaseline']
    #expname_arr = ['s4wide_scaled_aso_plus_fulls4scaledsobaseline', 's4wide_plus_fulls4scaledsobaseline', 's4wide_single_chlat_plus_aso_plus_fulls4scaledsobaseline']
    #expname_arr = ['s4wide_single_chlat_plus_2028aso_plus_fulls4scaledsobaseline', 's4wide_single_chlat_plus_2028aso']
    expname_arr = ['s4deepv3r025']

    if (1): #20220726 - all experiments.
        expname_arr = ['s4wide', 's4deepv3r025', 'sobaseline', 'sogoal', 'spt3g']#, 'spt4_C4', 'ccat_prime_so']

    if (1): #20230517 - for S4 Neff paper
        expname_arr = ['s4wide', 's4deepv3r025']

    s4_so_joint_configs = 1
    include_gal = 1
    which_gal_mask = 2
    interactive_mode = 0
    save_fg_res_and_weights = 0
    total = 0
    for expname in expname_arr:
        if expname == 's4deepv3r025':
            include_gal = 0
            which_gal_mask = -1
        if expname.find('plus_fulls4scaledsobaseline')>-1:
            include_fulls4scaledsobaseline = 1
        else:
            include_fulls4scaledsobaseline = 0
        if expname == 's4wide_scaled_sobaseline' or expname == 's4wide_scaled_aso' or expname == 's4wide_scaled_aso_plus_fulls4scaledsobaseline':
            total_obs_time_arr = np.arange(1., 5.1, 1.)
        else:
            total_obs_time_arr = np.arange(1., 10.1, 1.)

        if (0):
            total_obs_time_arr = [7.]
            include_fulls4scaledsobaseline = 0
            include_gal = 0
            s4_so_joint_configs = 0
            save_fg_res_and_weights = 1
            which_gal_mask = -1

        if (1): #Planck masks
            ##total_obs_time_arr = [7.]
            #total_obs_time_arr = np.arange(1., 7.1, 1.)
            total_obs_time_arr = np.arange(8., 10.1, 1.)
            include_fulls4scaledsobaseline = 0
            include_gal = 1 ##[0, 1]
            s4_so_joint_configs = 0
            save_fg_res_and_weights = 1
            which_gal_mask = 2 ##-1
            if expname == 's4wide':
                include_gal = 1 
                which_gal_mask = 2
            elif expname == 's4deepv3r025':
                include_gal = 0
                which_gal_mask = -1


        for total_obs_time in total_obs_time_arr:
            cmd = 'python3 %s -expname %s -include_gal %s -which_gal_mask %s -total_obs_time %s -s4_so_joint_configs %s -include_fulls4scaledsobaseline %s -interactive_mode %s -save_fg_res_and_weights %s' %(pgmname, expname, include_gal, which_gal_mask, total_obs_time, s4_so_joint_configs, include_fulls4scaledsobaseline, interactive_mode, save_fg_res_and_weights)
            print('\n###############\n%s\n' %(cmd)); ##sys.exit()
            os.system(cmd)
            ###sys.exit()
            total += 1
    print(total)
sys.exit()


