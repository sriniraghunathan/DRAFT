import numpy as np, sys, os

pgmname = 'get_ilc_residuals_spt.py'


if (1):
    exparr = ['spt3g', 'spt4_C1', 'spt4_C2', 'spt4_C3', 'spt4_C4', 'spt4_C5']

    exparr = ['spt3gplusherschel', 'spt4_C1plusherschel', 'spt4_C5']
    null_comp_arr = [None, 'cib', 'cmb', 'cib cmb']
    null_comp_arr = ['cib cmb']
    final_comp = 'y'

    for expname in exparr:
        for null_comp in null_comp_arr:
            cmd = 'python %s -expname %s -null_comp %s -final_comp %s ' %(pgmname, expname, null_comp, final_comp)
            print('\n###############\n%s\n' %(cmd))
            os.system(cmd)
            #sys.exit()
sys.exit()


