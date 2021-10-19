import numpy as np

def get_exp_specs(expname, remove_atm = 0, corr_noise_for_spt = 1):

    if expname.find('s4')>-1 or expname.find('cmbhd')>-1 or expname.find('cmb-hd')>-1:
        if expname == 's4wide':

            specs_dic = {
            #freq: [beam_arcmins, white_noise_T, elknee_T, alphaknee_T, whitenoise_P, elknee_P, alphaknee_P] 
            27: [7.4, 21.34, 415., 3.5, 30.23, 700, 1.4],
            39: [5.1, 11.67, 391., 3.5, 16.53, 700, 1.4], 
            93: [2.2, 1.89, 1932., 3.5, 2.68, 700, 1.4],
            145: [1.4, 2.09, 3917., 3.5, 2.96, 700, 1.4],
            225: [1.0, 6.90, 6740., 3.5, 9.78, 700, 1.4],
            278: [0.9, 16.88, 6792., 3.5, 23.93, 700, 1.4],
            }

            if remove_atm:
                specs_dic = {
                #freq: [beam_arcmins, white_noise_T, elknee_T, alphaknee_T, whitenoise_P, elknee_P, alphaknee_P] 
                27: [7.4, 21.34, 415., 0., 30.23, 700, 0.],
                39: [5.1, 11.67, 391., 0., 16.53, 700, 0.], 
                93: [2.2, 1.89, 1932., 0., 2.68, 700, 0.],
                145: [1.4, 2.09, 3917., 0., 2.96, 700, 0.],
                225: [1.0, 6.90, 6740., 0., 9.78, 700, 0.],
                278: [0.9, 16.88, 6792., 0., 23.93, 700, 0.],
                }

            '''
            #DSR specs
            specs_dic = {
            #freq: [beam_arcmins, white_noise_T, elknee_T, alphaknee_T, whitenoise_P, elknee_P, alphaknee_P] 
            #20: [10.0, None, None, None, None, None, None],
            27: [7.4, 21.8, 471., 3.5, 30.8, 700, 1.4],
            39: [5.1, 12.4, 428., 3.5, 17.6, 700, 1.4], 
            93: [2.2, 2.0, 2154., 3.5, 2.9, 700, 1.4],
            145: [1.4, 2.0, 4364., 3.5, 2.8, 700, 1.4],
            225: [1.0, 6.9, 7334., 3.5, 9.8, 700, 1.4],
            278: [0.9, 16.7, 7308., 3.5, 23.6, 700, 1.4],
            #225: [1.0, 100., 7334., 3.5, 9.8, 700, 1.4],
            #278: [0.9, 100., 7308., 3.5, 23.6, 700, 1.4],
            }

            if remove_atm:
                specs_dic = {
                #freq: [beam_arcmins, white_noise_T, elknee_T, alphaknee_T, whitenoise_P, elknee_P, alphaknee_P] 
                #20: [10.0, None, None, None, None, None, None],
                27: [7.4, 21.8, 471., 0., 30.8, 700, 0.],
                39: [5.1, 12.4, 428., 0., 17.6, 700, 0.], 
                93: [2.2, 2.0, 2154., 0., 2.9, 700, 0.],
                145: [1.4, 2.0, 4364., 0., 2.8, 700, 0.],
                225: [1.0, 6.9, 7334., 0., 9.8, 700, 0.],
                278: [0.9, 16.7, 7308., 0., 23.6, 700, 0.],
                }
            '''

        elif expname == 's4wide_chlat_el40':
            #https://cmb-s4.atlassian.net/wiki/spaces/XC/pages/680853505/Neff+forecasts+for+CHLAT+for+different+observing+elevations
            specs_dic = {
            #freq: [beam_arcmins, white_noise_T, elknee_T, alphaknee_T, whitenoise_P, elknee_P, alphaknee_P] 
            27: [7.3, 12.96, 415., 3.5, 18.33, 700, 1.4],
            39: [5.5, 9.83, 391., 3.5, 13.90, 700, 1.4], 
            93: [2.3, 1.49, 1932., 3.5, 2.11, 700, 1.4],
            145: [1.5, 1.39, 3917., 3.5, 1.97, 700, 1.4],
            225: [1.0, 3.44, 6740., 3.5, 4.87, 700, 1.4],
            278: [0.8, 5.44, 6792., 3.5, 7.69, 700, 1.4],
            }


        elif expname.find('cmbhd')>-1 or expname.find('cmb-hd')>-1:

            specs_dic = {
            #freq: [beam_arcmins, white_noise_T, elknee_T, alphaknee_T, whitenoise_P, elknee_P, alphaknee_P] 
            27: [1.4, 6.5, 415., 3.5, 9.2, 700, 1.4],
            39: [1.05, 3.4, 391., 3.5, 4.8, 700, 1.4], 
            93: [0.45, 0.73, 1932., 3.5, 1.03, 700, 1.4],
            145: [0.25, 0.79, 3917., 3.5, 1.12, 700, 1.4],
            225: [0.2, 2.0, 6740., 3.5, 2.828, 700, 1.4],
            278: [0.15, 4.6, 6792., 3.5, 6.5, 700, 1.4],
            350: [0.12, 4.6, 6792., 3.5, 6.5, 700, 1.4],
            }

        elif expname == 's4deep':
            '''
            specs_dic = {
            #freq: [beam_arcmins, white_noise_T, elknee_T, alphaknee_T, whitenoise_P, elknee_P, alphaknee_P] 
            #20: [10.0, None, None, None, None, None, None],
            27: [7.4, 21.8, 471., 3.5, 30.8, 700, 1.4],
            39: [5.1, 12.4, 428., 3.5, 17.6, 700, 1.4], 
            93: [2.2, 0.48, 2154., 3.5, 0.68, 700, 1.4],
            145: [1.4, 0.67, 4364., 3.5, 0.96, 700, 1.4],
            225: [1.0, 4.04, 7334., 3.5, 5.72, 700, 1.4],
            278: [0.9, 6.92, 7308., 3.5, 9.8, 700, 1.4],
            }
            '''
            #different 1/f noise definitions from S4V3R0, I think.
            specs_dic = {
            #freq: [beam_arcmins, white_noise_T, elknee_T, alphaknee_T, whitenoise_P, elknee_P, alphaknee_P] 
            #20: [10.0, None, None, None, None, None, None],
            27: [7.4, 21.8, 1200., 4.2, 30.8, 700, 1.4],
            39: [5.1, 12.4, 1200., 4.2, 17.6, 700, 1.4], 
            93: [2.2, 0.48, 1200., 4.2, 0.68, 700, 1.4],
            145: [1.4, 0.67, 1900., 4.1, 0.96, 700, 1.4],
            225: [1.0, 4.04, 2100., 4.1, 5.72, 700, 1.4],
            278: [0.9, 6.92, 2100., 3.9, 9.8, 700, 1.4],
            }

            '''
            specs_dic = {
            #freq: [beam_arcmins, white_noise_T, elknee_T, alphaknee_T, whitenoise_P, elknee_P, alphaknee_P] 
            #20: [10.0, None, None, None, None, None, None],
            #27: [7.4, 21.8, 471., 3.5, 30.8, 700, 1.4],
            #39: [5.1, 12.4, 428., 3.5, 17.6, 700, 1.4], 
            93: [2.2, 0.48, 2154., 0., 0.68, 700, 1.4],
            145: [1.4, 0.67, 4364., 0., 0.96, 700, 1.4],
            225: [1.0, 4.04, 7334., 0., 5.72, 700, 1.4],
            278: [0.9, 6.92, 7308., 0., 9.8, 700, 1.4],
            }
            '''

        elif expname == 's4deepv3r025':
            #https://cmb-s4.org/wiki/index.php/Delensing_sensitivity_-_updated_sensitivities,_beams,_TT_noise
            specs_dic = {
            #freq: [beam_arcmins, white_noise_T, elknee_T, alphaknee_T, whitenoise_P, elknee_P, alphaknee_P] 
            #20: [11.0, None, None, None, None, None, None],
            27: [8.4, 4.6, 1200., 4.2, 6.5, 150, 2.7],
            39: [5.8, 2.94, 1200., 4.2, 4.15, 150, 2.7], 
            93: [2.5, 0.45, 1200., 4.2, 0.63, 150, 2.6],
            145: [1.6, 0.41, 1900., 4.1, 0.59, 200, 2.2],
            225: [1.1, 1.29, 2100., 4.1, 1.83, 200, 2.2],
            278: [1.0, 3.07, 2100., 3.9, 4.34, 200, 2.2],
            }

        freqarr = sorted( specs_dic.keys() )
        nc = len( freqarr )

        corr_noise = 1
        if corr_noise:
            corr_noise_bands = {20: [20], 27:[39], 39:[27], 93:[145], 145:[93], 225: [278], 278: [225], 350: [350]}
        else:
            corr_noise_bands = {20: [20], 27:[39], 39:[27], 93:[93], 145:[145], 225: [225], 278: [278], 350: [350]}
        rho = 0.9

    elif expname.find('spt')>-1:

        elknee_T_90, elknee_T_150, elknee_T_220 = 1200., 2200., 2300.
        alphaknee_T_90, alphaknee_T_150, alphaknee_T_220 = 3., 4., 4.
        elknee_P_90, elknee_P_150, elknee_P_220 = 300., 300., 300.
        alphaknee_P_90, alphaknee_P_150, alphaknee_P_220 = 1., 1., 1.

        if remove_atm:
            alphaknee_T_90, alphaknee_T_150, alphaknee_T_220 = 0., 0., 0.
            alphaknee_P_90, alphaknee_P_150, alphaknee_P_220 = 0., 0., 0.

        beam_90, beam_150, beam_220 = 1.7, 1.2, 1.0

        if expname == 'sptsz':
            #slide 3 of Lindsey https://pole.uchicago.edu/sptpol/images/Spt_des_clusters_july17_f2f.pdf
            white_noise_T_90, white_noise_T_150, white_noise_T_220 = 40., 17., 80.
            white_noise_P_90, white_noise_P_150, white_noise_P_220 = 1e6, 1e6, 1e6
        elif expname == 'sptpolultradeep':
            #slide 3 of Lindsey https://pole.uchicago.edu/sptpol/images/Spt_des_clusters_july17_f2f.pdf
            white_noise_T_90, white_noise_T_150, white_noise_T_220 = 10., 3.8, 1e6

            white_noise_P_90 = white_noise_T_90 * np.sqrt(2.)
            white_noise_P_150 = white_noise_T_150 * np.sqrt(2.)
            white_noise_P_220 = white_noise_T_220 * np.sqrt(2.)
        elif expname == 'sptpolultradeepplus3g':
            #slide 3 of Lindsey https://pole.uchicago.edu/sptpol/images/Spt_des_clusters_july17_f2f.pdf
            white_noise_T_90, white_noise_T_150, white_noise_T_220 = 10., 3.8, 30.

            white_noise_P_90 = white_noise_T_90 * np.sqrt(2.)
            white_noise_P_150 = white_noise_T_150 * np.sqrt(2.)
            white_noise_P_220 = white_noise_T_220 * np.sqrt(2.)
        elif expname == 'sptpolplusultradeep':
            #slide 3 of Lindsey https://pole.uchicago.edu/sptpol/images/Spt_des_clusters_july17_f2f.pdf
            white_noise_T_90, white_noise_T_220 = 10., 1e6
            white_noise_T_150 = ( (1./3.8**2.) + (1./5**2.) )**-0.5

            white_noise_P_90 = white_noise_T_90 * np.sqrt(2.)
            white_noise_P_150 = white_noise_T_150 * np.sqrt(2.)
            white_noise_P_220 = white_noise_T_220 * np.sqrt(2.)
        elif expname == 'sptpolplusultradeepplus3g':
            #slide 3 of Lindsey https://pole.uchicago.edu/sptpol/images/Spt_des_clusters_july17_f2f.pdf
            white_noise_T_90, white_noise_T_220 = 10., 30.
            white_noise_T_150 = ( (1./3.8**2.) + (1./5**2.) )**-0.5

            white_noise_P_90 = white_noise_T_90 * np.sqrt(2.)
            white_noise_P_150 = white_noise_T_150 * np.sqrt(2.)
            white_noise_P_220 = white_noise_T_220 * np.sqrt(2.)
        elif expname == 'sptpolplusultradeepplus3gfull':
            #slide 3 of Lindsey https://pole.uchicago.edu/sptpol/images/Spt_des_clusters_july17_f2f.pdf
            white_noise_T_90, white_noise_T_220 = 10., 10.
            white_noise_T_150 = ( (1./3.8**2.) + (1./5**2.) )**-0.5

            white_noise_P_90 = white_noise_T_90 * np.sqrt(2.)
            white_noise_P_150 = white_noise_T_150 * np.sqrt(2.)
            white_noise_P_220 = white_noise_T_220 * np.sqrt(2.)

        elif expname == 'sptpolsummer':
            white_noise_T_150 = 28.
            white_noise_T_90 = white_noise_T_150 * np.sqrt(2.)
            white_noise_T_220 = 1e6

            white_noise_P_90 = white_noise_T_90 * np.sqrt(2.)
            white_noise_P_150 = white_noise_T_150 * np.sqrt(2.)
            white_noise_P_220 = white_noise_T_220 * np.sqrt(2.)

        elif expname == 'spt3gsummer':
            white_noise_T_90 = 15.
            white_noise_T_150 = 15.
            white_noise_T_220 = 40.

            white_noise_P_90 = white_noise_T_90 * np.sqrt(2.)
            white_noise_P_150 = white_noise_T_150 * np.sqrt(2.)
            white_noise_P_220 = white_noise_T_220 * np.sqrt(2.)

        elif expname == 'spt3g':
            white_noise_T_90 = 3.0
            white_noise_T_150 = 2.2
            white_noise_T_220 = 8.8

            white_noise_P_90 = np.sqrt(2.) * white_noise_T_90
            white_noise_P_150 = np.sqrt(2.) * white_noise_T_150
            white_noise_P_220 = np.sqrt(2.) * white_noise_T_220

        elif expname == 'spt3g_ZP':
            white_noise_T_90 = 14.3
            white_noise_T_150 = 14.0
            white_noise_T_220 = 48.

            white_noise_P_90 = 25.
            white_noise_P_150 = 17.
            white_noise_P_220 = 56.

        elif expname == 'spt3g_y12':
            white_noise_T_90 = 9.
            white_noise_T_150 = 10.
            white_noise_T_220 = 40.

            white_noise_P_90 = np.sqrt(2.) * white_noise_T_90
            white_noise_P_150 = np.sqrt(2.) * white_noise_T_150
            white_noise_P_220 = np.sqrt(2.) * white_noise_T_220

        elif expname == 'spt3g_summer':
            white_noise_T_90 = 13.
            white_noise_T_150 = 13.0
            white_noise_T_220 = 40.

            white_noise_P_90 = np.sqrt(2.) * white_noise_T_90
            white_noise_P_150 = np.sqrt(2.) * white_noise_T_150
            white_noise_P_220 = np.sqrt(2.) * white_noise_T_220

        elif expname == 'sptpol_summer':
            white_noise_T_90 = 60.
            white_noise_T_150 = 28.0
            white_noise_T_220 = 10000.

            white_noise_P_90 = np.sqrt(2.) * white_noise_T_90
            white_noise_P_150 = np.sqrt(2.) * white_noise_T_150
            white_noise_P_220 = np.sqrt(2.) * white_noise_T_220

        elif expname == 'sptpol':
            white_noise_T_90 = 14.
            white_noise_T_150 = 6.0
            white_noise_T_220 = 10000.

            white_noise_P_90 = np.sqrt(2.) * white_noise_T_90
            white_noise_P_150 = np.sqrt(2.) * white_noise_T_150
            white_noise_P_220 = np.sqrt(2.) * white_noise_T_220

        elif expname == 'spt3g_TC':
            white_noise_T_90 = 7.960
            white_noise_T_150 = 6.330
            white_noise_T_220 = 23.39

            white_noise_P_90 = 11.26
            white_noise_P_150 = 8.950
            white_noise_P_220 = 33.08

        specs_dic = {
        90: [beam_90, white_noise_T_90, elknee_T_90, alphaknee_T_90, white_noise_P_90, elknee_P_90, alphaknee_P_90],
        150: [beam_150, white_noise_T_150, elknee_T_150, alphaknee_T_150, white_noise_P_150, elknee_P_150, alphaknee_P_150],
        220: [beam_220, white_noise_T_220, elknee_T_220, alphaknee_T_220, white_noise_P_220, elknee_P_220, alphaknee_P_220],
        }
        
        corr_noise = corr_noise_for_spt
        if corr_noise:
            corr_noise_bands = {90:[90], 150:[220], 220: [150]}
        else:
            corr_noise_bands = {90:[90], 150:[150], 220: [220]}

        rho = 1.0

    return specs_dic, corr_noise_bands, rho, corr_noise
