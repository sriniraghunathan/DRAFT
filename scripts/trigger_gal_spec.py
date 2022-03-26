import numpy as np, sys, os

pgmname = 'get_pspec_galactic_sims.py'
use_planck_mask = 0
use_lat_step_mask = 0
use_s4like_mask = 0
use_s4like_mask_v2 = 0
use_s4like_mask_v3 = 0
use_spt3g_mask = 0
use_s4delensing_mask = 0
##s,e = 0, 4
zonca_sims = 1
pySM_yomori = 0
mask_arr_default = None
min_obs_el = 30.0
use_splat_minobsel_galcuts =0

if (0):
    use_lat_step_mask = 1
    t_only = 1 ##0
    nside = 4096 ##2048
    lmax = 10000 ##7000
    s,e = 3, 4


if (0):
    use_s4like_mask = 1
    t_only = 0
    nside = 4096 ##2048
    lmax = 7000
    s,e = 0, 3

if (1):
    use_s4delensing_mask = 0
    #use_s4like_mask_v2 = 1
    #s,e = 0, 6#3, 5

    use_s4like_mask_v2 = 1
    s,e = 0, 6#3, 5
    #mask_arr_default = [2, 5]

    #delensing LAT
    use_s4like_mask_v2 = 0
    use_s4delensing_mask = 1
    s,e = 0, 1

    use_s4like_mask_v3 = 0
    #s,e = 0, 3

    t_only = 0
    #nside = 4096 ##2048
    #lmax = 7000
    nside = 2048
    lmax = 5000

if (0): ### 20210512
    use_spt3g_mask = 1
    t_only = 0
    nside = 2048
    lmax = 6000
    s,e = 0, 4

    use_planck_mask = 0
    use_lat_step_mask = 0
    use_s4like_mask = 0
    use_s4like_mask_v2 = 0
    use_s4like_mask_v3 = 0
    use_s4delensing_mask = 0
    splat_minobsel_galcuts = 1
    min_obs_el = 30.0

dust_sync_arr = ['dust', 'sync']

if (0): #pySM YOmori for SPT-3G + Planck
    use_spt3g_mask = 1
    t_only = 0
    nside = 2048
    lmax = 2000 ##5000 ###6000
    s,e = 0, 1 #only the main winter field
    zonca_sims = 1 ###0
    pySM_yomori = 0 ##1
    dust_sync_arr = ['dust', 'sync', 'freefree']
    #dust_sync_arr = ['both']

if (0): #20210111 - CMB-S4 SP-LAT - different min obs el and gal cuts
    use_spt3g_mask = 0
    t_only = 0
    min_obs_el = 30.
    use_splat_minobsel_galcuts = 1
    nside = 2048
    lmax = 5000
    s,e = 0, 6 #different gal cuts #check https://github.com/sriniraghunathan/cmbs4_minobsel_galcuts/blob/main/analyse_results.ipynb
    zonca_sims = 1
    dust_sync_arr = ['dust', 'sync']
    
mask_arr = np.arange(s, e)
#if mask_arr_default is not None:
#    mask_arr = mask_arr_default

template_fname = 'batch_jobs/template_hoff.sh'

for which_mask in mask_arr:
    for dust_or_sync in dust_sync_arr:
        opfname = 'batch_jobs/gspec_mask%s_%s_%s.sh' %(which_mask, dust_or_sync, lmax)

        opf = open(opfname, 'w')
        template = open(template_fname, 'r')

        for line in template:
            opf.writelines('%s\n' %(line.strip()))

        #opline = 'python %s -dust_or_sync %s -zonca_sims %s -pySM_yomori %s -which_mask %s -use_planck_mask %s -use_lat_step_mask %s -use_s4like_mask %s -use_s4like_mask_v2 %s -use_s4like_mask_v3 %s -use_spt3g_mask %s -t_only %s -nside %s -lmax %s ' %(pgmname, dust_or_sync, zonca_sims, pySM_yomori, int(which_mask), use_planck_mask, use_lat_step_mask, use_s4like_mask, use_s4like_mask_v2, use_s4like_mask_v3, use_spt3g_mask, t_only, nside, lmax)
        #opline = 'python %s -dust_or_sync %s -zonca_sims %s -pySM_yomori %s -which_mask %s -use_planck_mask %s -use_lat_step_mask %s -use_s4like_mask %s -use_s4like_mask_v2 %s -use_s4like_mask_v3 %s -use_s4delensing_mask %s -use_spt3g_mask %s -t_only %s -nside %s -lmax %s ' %(pgmname, dust_or_sync, zonca_sims, pySM_yomori, int(which_mask), use_planck_mask, use_lat_step_mask, use_s4like_mask, use_s4like_mask_v2, use_s4like_mask_v3, use_s4delensing_mask, use_spt3g_mask, t_only, nside, lmax)
        #20210111 - CMB-S4 SP-LAT - different min obs el and gal cuts
        opline = 'python %s -dust_or_sync %s -zonca_sims %s -pySM_yomori %s -which_mask %s -use_planck_mask %s -use_lat_step_mask %s -use_s4like_mask %s -use_s4like_mask_v2 %s -use_s4like_mask_v3 %s -use_s4delensing_mask %s -use_spt3g_mask %s -t_only %s -nside %s -lmax %s -min_obs_el %s -use_splat_minobsel_galcuts %s ' %(pgmname, dust_or_sync, zonca_sims, pySM_yomori, int(which_mask), use_planck_mask, use_lat_step_mask, use_s4like_mask, use_s4like_mask_v2, use_s4like_mask_v3, use_s4delensing_mask, use_spt3g_mask, t_only, nside, lmax, min_obs_el, use_splat_minobsel_galcuts)
        #print(opline); sys.exit()
        opf.writelines('%s\n\n' %(opline))
        opf.close()
        template.close()

        cmd = 'qsub -V %s' %(opfname)
        os.system(cmd)
        print('%s\n' %(cmd))

