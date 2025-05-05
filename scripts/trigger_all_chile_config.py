#20250504
#All Chile configuration

import numpy as np, os, sys

cmd = 'python3 get_ilc_residuals.py -expname exppatchval -include_gal galval -which_gal_mask galmaskval -total_obs_time totalyearval -save_fg_res_and_weights 0 -final_comp cmb -interactive_mode 0'
total_year_arr = np.arange(1, 10.1, 1)
survey_arr = ['lat_wide', 'lat_delensing', 'lat_roman']
patch_arr = [1, 2, 3, 4]
galval = 0
galmaskval = 0

if (0): #with galaxy
    survey_arr = ['lat_wide']
    galval = 1
    galmaskval = 2

for survey in survey_arr:
    for patch in patch_arr:
        for total_year in total_year_arr:
            if total_year == 10.: continue
            exppatchval = 's4_all_chile_config_%s---patch%s' %(survey, patch)
            curr_cmd = cmd.replace('exppatchval', exppatchval).replace('galval', str(galval)).replace('galmaskval', str(galmaskval)).replace('totalyearval', str(total_year))
            print('\n', curr_cmd)
            os.system(curr_cmd)
