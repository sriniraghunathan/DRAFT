#20250504
#All Chile configuration

import numpy as np, os, sys


#cmd = 'python3 get_ilc_residuals.py -expname exppatchval -include_gal galval -which_gal_mask galmaskval -total_obs_time totalyearval -save_fg_res_and_weights 0 -final_comp cmb -interactive_mode 0'
cmd = 'python3 get_ilc_residuals.py -expname exppatchval -include_gal galval -which_gal_mask galmaskval -save_fg_res_and_weights 0 -final_comp cmb -interactive_mode 0'
total_year_arr = np.arange(1, 10.1, 1)
survey_arr = ['lat_wide', 'lat_delensing', 'advanced_so_goal']

aso_s4_start_year_diff = (2033 - 2028)
total = 0
for survey in survey_arr:
    if survey in ['lat_wide']:
        patch_arr = [2]
        galval, galmaskval = 1, 2
    elif survey in ['lat_delensing']:
        patch_arr = [1]
        galval, galmaskval = 0, 0
    elif survey in ['advanced_so_goal']:
        patch_arr = [-1]
        galval, galmaskval = 0, 0

    for patch in patch_arr:
        for total_year in total_year_arr:
            if survey == 'lat_wide': #combine with ASO
                total_year_for_aso = total_year + aso_s4_start_year_diff
                exppatchval = 's4_all_chile_config_%s---patch%s---year%s+advanced_so_goal---year%s' %(survey, patch, total_year, total_year_for_aso)
            elif survey == 'lat_delensing': #delensing
                exppatchval = 's4_all_chile_config_%s---patch%s---year%s' %(survey, patch, total_year)
            elif survey == 'advanced_so_goal':
                exppatchval = '%s---year%s' %(survey, total_year)

            curr_cmd = cmd.replace('exppatchval', exppatchval).replace('galval', str(galval)).replace('galmaskval', str(galmaskval))#.replace('totalyearval', str(total_year))
            print('\n', curr_cmd); ##sys.exit()
            ##os.system(curr_cmd)
            total += 1

print('\nTotal = %s.\n' %(total))    
sys.exit()

#------------- older runs below.

cmd = 'python3 get_ilc_residuals.py -expname exppatchval -include_gal galval -which_gal_mask galmaskval -total_obs_time totalyearval -save_fg_res_and_weights 0 -final_comp cmb -interactive_mode 0'
total_year_arr = np.arange(1, 10.1, 1)
survey_arr = ['lat_wide', 'lat_wide_phase2', 'lat_wide_dc0', 'lat_delensing', 'lat_roman']
#survey_arr = ['lat_wide_phase2', 'lat_wide_dc0']
survey_arr = ['lat_wide', 'lat_wide_dc0']
galval = 0
galmaskval = 0
combined_with_advanced_so_arr = [None]




if (0): #with galaxy
    #survey_arr = ['lat_wide']
    #survey_arr = ['lat_wide_phase2', 'lat_wide_dc0']
    survey_arr = ['lat_wide', 'lat_wide_dc0']
    galval = 1
    galmaskval = 2

if (1): #Combined with advanced SO
    combined_with_advanced_so_arr = ['advanced_so_baseline', 'advanced_so_goal']

total = 0
for survey in survey_arr:
    if survey in ['lat_wide', 'lat_delensing', 'lat_roman']:
        patch_arr = [1, 2, 3, 4]
    elif survey in ['lat_wide_dc0']:
        patch_arr = [1, 2]
    elif survey in ['lat_wide_phase2']:
        patch_arr = [1, 2, 3]

    for patch in patch_arr:
        for total_year in total_year_arr:
            for curr_advanced_so_config in combined_with_advanced_so_arr:
                if curr_advanced_so_config is None:
                    ###if total_year == 10.: continue
                    exppatchval = 's4_all_chile_config_%s---patch%s' %(survey, patch)
                else:
                    exppatchval = 's4_all_chile_config_%s---patch%s+%s' %(survey, patch, curr_advanced_so_config)
                curr_cmd = cmd.replace('exppatchval', exppatchval).replace('galval', str(galval)).replace('galmaskval', str(galmaskval)).replace('totalyearval', str(total_year))
                print('\n', curr_cmd); ###sys.exit()
                os.system(curr_cmd)
                total += 1

print('\nTotal = %s.\n' %(total))

"""
if (1): #Advanced-SO
    expname_arr = ['advanced_so_baseline', 'advanced_so_goal']
    galval = 0
    galmaskval = 0

    for expname in expname_arr:

        curr_cmd = cmd.replace('exppatchval', expname).replace('galval', str(galval)).replace('galmaskval', str(galmaskval)).replace(' -total_obs_time totalyearval', '')
        print('\n', curr_cmd)
        os.system(curr_cmd)
"""

"""
if (1): #SPT-3G+
    expname_arr = ['spt3g_plus_spt3g+_WG2']
    galval = 0
    galmaskval = 0

    for expname in expname_arr:

        curr_cmd = cmd.replace('exppatchval', expname).replace('galval', str(galval)).replace('galmaskval', str(galmaskval)).replace(' -total_obs_time totalyearval', '')
        print('\n', curr_cmd)
        os.system(curr_cmd)
"""

