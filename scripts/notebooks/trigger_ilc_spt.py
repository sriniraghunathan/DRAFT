import numpy as np, sys, os

pgmname = 'get_ilc_residuals_spt.py'


if (0):
    #exparr = ['spt3g', 'spt4_C1', 'spt4_C2', 'spt4_C3', 'spt4_C4', 'spt4_C5']
    exparr = ['spt3g', 'spt4_C5_v2']
    #exparr = ['spt3gplusherschel', 'spt4_C1plusherschel', 'spt4_C5']
    null_comp_arr = [None, 'cib', 'cmb', 'cib cmb']
    #null_comp_arr = ['cib cmb']
    use_websky_cib_arr = [0]##, 1]
    final_comp = 'y'

    for expname in exparr:
        for use_websky_cib in use_websky_cib_arr:
            for null_comp in null_comp_arr:
                cmd = 'python3 %s -expname %s -null_comp %s -final_comp %s -use_websky_cib %s' %(pgmname, expname, null_comp, final_comp, use_websky_cib)
                print('\n###############\n%s\n' %(cmd))
                os.system(cmd)
                #sys.exit()

if (1):
    #exparr = ['spt3g', 'spt4_C1', 'spt4_C2', 'spt4_C3', 'spt4_C4', 'spt4_C5']
    exparr = ['spt3g', 'spt4_C3']
    #exparr = ['spt3g', 'sptpolplusultradeepplus3gplusherschel']#, 'spt4_C3'] #CIB at spt4 bands are not computed for websky
    exparr = ['sptpolplusultradeepplus3gplusherschel']#, 'spt4_C3'] #CIB at spt4 bands are not computed for websky
    #null_comp_arr = [None, 'y', 'cib', 'cib y']
    #null_comp_arr = ['cibpo cibclus y']
    null_comp_arr = ['cibpo cibclus']
    use_websky_cib_arr = [1]#0, 1]
    final_comp = 'cmb'

    for expname in exparr:
        for use_websky_cib in use_websky_cib_arr:
            for null_comp in null_comp_arr:
                cmd = 'python3 %s -expname %s -null_comp %s -final_comp %s -use_websky_cib %s' %(pgmname, expname, null_comp, final_comp, use_websky_cib)
                print('\n###############\n%s\n' %(cmd))
                os.system(cmd)
                #sys.exit()

sys.exit()


