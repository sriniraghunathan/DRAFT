import numpy as np, sys, os

pgmname = 'get_pspec_galactic_sims.py'
use_planck_mask = 0
use_lat_step_mask = 0
use_s4like_mask = 0
use_s4like_mask_v2 = 0
use_spt3g_mask = 0
##s,e = 0, 4
zonca_sims = 1
pySM_yomori = 0
mask_arr_default = None

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
    #use_s4like_mask_v2 = 1
    #s,e = 0, 6#3, 5

    use_s4like_mask_v2 = 1
    s,e = 0, 6#3, 5
    #mask_arr_default = [2, 5]

    use_s4like_mask_v3 = 0
    #s,e = 0, 3

    t_only = 0
    #nside = 4096 ##2048
    #lmax = 7000
    nside = 2048
    lmax = 5000

if (0):
    use_spt3g_mask = 1
    t_only = 0
    nside = 2048
    lmax = 6000
    s,e = 0, 4

dust_sync_arr = ['dust', 'sync']

if (0): #pySM YOmori for SPT-3G + Planck
    use_spt3g_mask = 1
    t_only = 0
    nside = 2048
    lmax = 6000
    s,e = 0, 1 #only the main winter field
    zonca_sims = 0
    pySM_yomori = 1
    #dust_sync_arr = ['both']

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

        opline = 'python %s -dust_or_sync %s -zonca_sims %s -pySM_yomori %s -which_mask %s -use_planck_mask %s -use_lat_step_mask %s -use_s4like_mask %s -use_s4like_mask_v2 %s -use_s4like_mask_v3 %s -use_spt3g_mask %s -t_only %s -nside %s -lmax %s ' %(pgmname, dust_or_sync, zonca_sims, pySM_yomori, int(which_mask), use_planck_mask, use_lat_step_mask, use_s4like_mask, use_s4like_mask_v2, use_s4like_mask_v3, use_spt3g_mask, t_only, nside, lmax)
        opf.writelines('%s\n\n' %(opline))
        opf.close()
        template.close()

        cmd = 'qsub -V %s' %(opfname)
        os.system(cmd)
        print('%s\n' %(cmd))

