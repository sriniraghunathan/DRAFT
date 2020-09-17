import numpy as np, sys, os

s,e = 0, 7
mask_arr = np.arange(s, e+1)
dust_sync_arr = ['dust', 'sync']

pgmname = 'get_pspec_galactic_sims.py'
use_planck_mask = 0
use_lat_step_mask = 0
use_spt3g_mask = 0
use_s4like_mask_v2 = 0
if (1): #spt 3G summer/ winter field
    use_spt3g_mask = 1
    mask_arr = [0, 1]

template_fname = 'template_hoff.sh'

for which_mask in mask_arr:
    for dust_or_sync in dust_sync_arr:
        opfname = 'batch_jobs/gspec_mask%s_%s.sh' %(which_mask, dust_or_sync)

        opf = open(opfname, 'w')
        template = open(template_fname, 'r')

        for line in template:
            opf.writelines('%s\n' %(line))

        opline = 'python %s -dust_or_sync %s -which_mask %s -use_planck_mask %s -use_lat_step_mask %s -use_s4like_mask_v2 %s -use_spt3g_mask %s' %(pgmname, which_mask, dust_or_sync, use_planck_mask, use_lat_step_mask, use_s4like_mask_v2, use_spt3g_mask)
        opf.writelines('%s\n' %(opline))
        opf.close()
        template.close()

        cmd = 'qsub -V %s' %(opfname)
        os.system(cmd)
        print(cmd)
