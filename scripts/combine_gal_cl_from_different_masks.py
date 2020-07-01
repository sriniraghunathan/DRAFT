import numpy as np, glob, sys, os, re

local = 1
if str(os.getcwd()).find('sri')>-1:
    local =0

s4like_mask =0
lat_steps_mask=0

nside, lmax = 2048, 3500
t_only = 0

zonca_sims = 1

'''
lat_steps_mask = 1
nside, lmax = 4096, 10000 ##7000
'''

#s4like_mask = 1
#nside, lmax = 4096, 7000
s4like_mask_v2 = 1
nside, lmax = 4096, 7000
nside, lmax = 2048, 5000
t_only = 0


if zonca_sims:
    data_folder = '/Volumes/data_PHD_WD_babbloo/s4/cmbs4/map_based_simulations/202002_foregrounds_extragalactic_cmb_tophat/4096/'
    if not local:
        data_folder = '//u/home/s/srinirag/project-nwhiteho/cmbs4/map_based_simulations/202002_foregrounds_extragalactic_cmb_tophat/4096/'
    comp_arr = ['dust', 'synchrotron']
    #comp_arr = ['synchrotron']
else:
    comp_arr = ['dust', 'sync']

for which_comp in comp_arr:

    if not zonca_sims:
        searchstr = 'cls_galactic_sims_xxxx_CUmilta_20200319_maskplanck_nside%s_lmax%s_mask?.npy' %(nside, lmax)
    else:
        if lat_steps_mask:
            searchstr = '%s/%s/0000/lat_steps/cls_galactic_sims_xxxx_nside%s_lmax%s_TTonly_mask?.npy' %(data_folder, which_comp, nside, lmax)
        elif s4like_mask:
            if not t_only:
                searchstr = '%s/%s/0000/s4like_mask/cls_galactic_sims_xxxx_nside%s_lmax%s_mask?.npy' %(data_folder, which_comp, nside, lmax)
            else:
                searchstr = '%s/%s/0000/s4like_mask/cls_galactic_sims_xxxx_nside%s_lmax%s_TTonly_mask?.npy' %(data_folder, which_comp, nside, lmax)
	elif s4like_mask_v2:
	    searchstr = '%s/%s/0000/s4like_mask_v2/cls_galactic_sims_xxxx_nside%s_lmax%s_mask?.npy' %(data_folder, which_comp, nside, lmax)
            searchstr = searchstr.replace('.npy', '_cos_el_40.npy')
        else:
            searchstr = '%s/%s/0000/lat_steps/cls_galactic_sims_xxxx_maskplanck_nside%s_lmax%s_TTonly_mask?.npy' %(data_folder, which_comp, nside, lmax)

    #print(which_comp, searchstr)
    curr_searchstr = searchstr.replace('xxxx', which_comp)
    curr_searchstr = curr_searchstr.replace('_synchrotron_', '_sync_')
    flist = sorted( glob.glob(curr_searchstr) )

    opdic = {}
    opdic['cl_dic'] = {}
    fksy_arr = []
    print('\n')
    for fcntr, f in enumerate( flist ):
        ##print(f)
        #which_mask = int( f.split('_')[-1].replace('mask','').replace('.npy','') )
        mask_str = re.findall(r'mask[0-9]?', f.split('/')[-1])[0]
        which_mask = int(mask_str.replace('mask',''))
        curr_dic = np.load(f, allow_pickle = 1).item()

        opdic['cl_dic'][which_mask] = curr_dic['cl_dic'][0]

        lmax = curr_dic['lmax']
        fsky = curr_dic['fsky_arr'][0]
        fksy_arr.append( fsky )

        print(f.split('/')[-1],which_comp, which_mask, len( opdic['cl_dic'][which_mask].keys()), len( curr_dic['cl_dic'][0].keys()) )

    opdic['lmax'] = lmax
    opdic['fsky_arr'] = np.asarray( fksy_arr )

    opfname = f.replace('%s_' %mask_str, '')
    print(opfname)
    print('\n')
    np.save(opfname, opdic)

sys.exit()


