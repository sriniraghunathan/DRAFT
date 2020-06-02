import numpy as np, sys, os

pgmname = 'get_ilc_residuals_withcorrnoise_withgal_TTEETE_0.3.py'

totlf, totmf, tothf = 2, 12, 5
lf_to_mf_hf_arr = [(0,0), (0,1), (1,0), (1,1), (2,0), (0,2)]
mf_to_lf_hf_arr = [(0,0)]
hf_to_lf_mf_arr = [(0,0), (0,1), (1,0), (1,1), (2,0), (0,2)]#, (3,0), (4,0), (5,0)]
#also_te_arr = [ 0, 1]
also_te_arr = [ 1 ]
include_gal_arr = [0, 1]

for include_gal in include_gal_arr:
    if include_gal == 0:
        which_gal_mask_arr = [3]
    else:
        which_gal_mask_arr = [0, 1, 2, 3]
    for which_gal_mask in which_gal_mask_arr:
        for also_te in also_te_arr:
            tube_comb_arr = []
            for lmh in lf_to_mf_hf_arr:
                for mlh in mf_to_lf_hf_arr:
                    for hlm in hf_to_lf_mf_arr:
                        lf_to_mf, lf_to_hf = lmh
                        mf_to_lf, mf_to_hf = mlh
                        hf_to_lf, hf_to_mf = hlm

                        totlfmoved = lf_to_mf + lf_to_hf
                        totmfmoved = mf_to_lf + mf_to_hf
                        tothfmoved = hf_to_lf + hf_to_mf

                        newlf = totlf - totlfmoved + mf_to_lf + hf_to_lf
                        newmf = totmf - totmfmoved + lf_to_mf + hf_to_mf
                        newhf = tothf - tothfmoved + lf_to_hf + mf_to_hf

                        if newlf<0 or newmf<0 or newhf<0: continue

                        tube_comb = [newlf, newmf, newhf]
                        if tube_comb in tube_comb_arr: continue

                        tube_comb_arr.append( tube_comb )

                        #print(lf_to_mf, lf_to_hf, mf_to_lf, mf_to_hf, hf_to_lf, hf_to_mf)
                        #print(tube_comb)

                        cmd = 'python %s -lf_to_mf %s -lf_to_hf %s -mf_to_lf %s -mf_to_hf %s -hf_to_lf %s -hf_to_mf %s -also_te %s -include_gal %s -which_gal_mask %s ' %(pgmname, lf_to_mf, lf_to_hf, mf_to_lf, mf_to_hf, hf_to_lf, hf_to_mf, also_te, include_gal, which_gal_mask)
                        print('\n###############\n%s\n' %(cmd))
                        os.system(cmd)
                        #sys.exit()
sys.exit()


