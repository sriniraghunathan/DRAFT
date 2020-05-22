import numpy as np, sys, os

pgmname = 'get_pspec_galactic_sims.py'
use_planck_mask = 0
use_lat_step_mask = 0
use_s4like_mask = 0
##s,e = 0, 4

use_lat_step_mask = 1
t_only = 1 ##0
nside = 4096 ##2048
lmax = 10000 ##7000
s,e = 3, 4

use_s4like_mask = 1
t_only = 0
nside = 4096 ##2048
lmax = 7000
s,e = 0, 3

mask_arr = np.arange(s, e+0.1)
dust_sync_arr = ['dust', 'sync']


template_fname = 'batch_jobs/template_hoff.sh'

for which_mask in mask_arr:
    for dust_or_sync in dust_sync_arr:
        opfname = 'batch_jobs/gspec_mask%s_%s.sh' %(which_mask, dust_or_sync)

        opf = open(opfname, 'w')
        template = open(template_fname, 'r')

        for line in template:
            opf.writelines('%s\n' %(line.strip()))

        opline = 'python %s -dust_or_sync %s -which_mask %s -use_planck_mask %s -use_lat_step_mask %s -use_s4like_mask %s -t_only %s -nside %s -lmax %s ' %(pgmname, dust_or_sync, int(which_mask), use_planck_mask, use_lat_step_mask, use_s4like_mask, t_only, nside, lmax)
        opf.writelines('%s\n\n' %(opline))
        opf.close()
        template.close()

        cmd = 'qsub -V %s' %(opfname)
        os.system(cmd)
        print('%s\n' %(cmd))


