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
    exparr = ['spt3g']#, 'sptpolplusultradeepplus3gplusherschel']#, 'spt4_C3'] #CIB at spt4 bands are not computed for websky
    #exparr = ['sptpolplusultradeepplus3gplusherschel']#, 'spt4_C3'] #CIB at spt4 bands are not computed for websky
    #exparr = ['spt4_C3'] #CIB at spt4 bands are not computed for websky
    #exparr = ['spt4_C3_high220noise']
    #exparr = ['spt4_C3plusplanckHF'] #CIB at spt4 bands are not computed for websky

    #exparr = ['sptpolplusultradeepplus3gplusherschel_v2', 'sptpolplusultradeepplus3gplusherschel_v3', 'sptpolplusultradeepplus3gplusherschel_v4'] #CIB at spt4 bands are not computed for websky

    #null_comp_arr = [None, 'y', 'cib', 'cib y']
    #null_comp_arr = ['cibpo cibclus y']
    #null_comp_arr = ['cibpo cibclus']
    #null_comp_arr = ['misc_cib_tcib18_beta2.0']
    #null_comp_arr = ['misc_cib_tcib18_beta2.0 misc_cib_tcib20_beta1.5']

    #null_comp_arr = ['misc_cib_tcib20_beta2.0 misc_cib_tcib20_beta1.5']
    #null_comp_arr = ['misc_cib_tcib20_beta2.0 misc_cib_tcib20_beta1.5']
    #null_comp_arr = ['misc_cib_tcib18_beta1.5 misc_cib_tcib18_beta2.5']
    #null_comp_arr = ['misc_cib_tcib25_beta1.5 misc_cib_tcib25_beta2.5']
    null_comp_arr = ['misc_cib_tcib20_beta1.5 misc_cib_tcib20_beta2.5']

    #null_comp_arr = ['misc_cib_tcib20_beta1.5 misc_cib_tcib20_beta2.5 radio']
    null_comp_arr = [ null_comp_arr[0], 'y %s' %(null_comp_arr[0]) ]


    if exparr[0].find('spt4')>-1:
        null_comp_arr = ['y misc_cib_tcib20_beta1.5', 'y misc_cib_tcib25_beta1.5']


    #null_comp_arr = ['y', 'radio']
    split_cross = 0 ###1

    #null_comp_arr = ['cibpo cibclus y']

    use_websky_cib = 1 ## ###0 ##1
    use_sptspire_for_hfbands = 0 ###1 ##1 ##0
    use_mdpl2_cib = 0
    final_comp = 'cmb'

    final_comp = 'y'
    null_comp_arr = ['cmb misc_cib_tcib20_beta1.5 misc_cib_tcib20_beta2.5', 'misc_cib_tcib20_beta1.5 misc_cib_tcib20_beta2.5']

    if (1):
        final_comp = 'cmb'
        null_comp_arr = ['y']
        exparr = ['actd56', 'sptsz']

    for expname in exparr:
        for null_comp in null_comp_arr:
            cmd = 'python3 %s -expname %s -null_comp %s -final_comp %s -use_websky_cib %s -use_sptspire_for_hfbands %s -use_mdpl2_cib %s -split_cross %s' %(pgmname, expname, null_comp, final_comp, use_websky_cib, use_sptspire_for_hfbands, use_mdpl2_cib, split_cross)
            print('\n###############\n%s\n' %(cmd))
            #os.system(cmd)
            sys.exit()

sys.exit()


