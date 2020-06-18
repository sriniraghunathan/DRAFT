import numpy as np

def get_exp_specs(expname, remove_atm = 0):

    if expname == 's4':
        specs_dic = {
        #freq: [beam_arcmins, white_noise_T, elknee_T, alphaknee_T, whitenoise_P, elknee_P, alphaknee_P] 
        #20: [10.0, None, None, None, None, None, None],
        #27: [7.4, 21.8, 471., 3.5, 30.8, 700, 1.4],
        #39: [5.1, 12.4, 428., 3.5, 17.6, 700, 1.4], 
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
            #27: [7.4, 21.8, 471., 0., 30.8, 700, 0.],
            #39: [5.1, 12.4, 428., 0., 17.6, 700, 0.], 

            93: [2.2, 2.0, 2154., 0., 2.9, 700, 0.],
            145: [1.4, 2.0, 4364., 0., 2.8, 700, 0.],
            225: [1.0, 6.9, 7334., 0., 9.8, 700, 0.],
            278: [0.9, 16.7, 7308., 0., 23.6, 700, 0.],
        #    225: [1.0, 100., 7334., 0., 9.8, 700, 0.],
        #    278: [0.9, 100., 7308., 0., 23.6, 700, 0.],        
            }

        corr_noise = 1
        if corr_noise:
            corr_noise_bands = {20: [20], 27:[39], 39:[27], 93:[145], 145:[93], 225: [278], 278: [225]}
        else:
            corr_noise_bands = {20: [20], 27:[39], 39:[27], 93:[93], 145:[145], 225: [225], 278: [278]}
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
            white_noise_T_90, white_noise_T_150, white_noise_T_220 = 40., 18., 55.
            white_noise_P_90, white_noise_P_150, white_noise_P_220 = 1e6, 1e6, 1e6
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

        elif expname == 'spt3g_ZP':
            white_noise_T_90 = 14.3
            white_noise_T_150 = 14.0
            white_noise_T_220 = 48.

            white_noise_P_90 = 25.
            white_noise_P_150 = 17.
            white_noise_P_220 = 56.

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

        corr_noise_bands = {90:[90], 150:[150], 220: [220]}
        rho = 1.

    return specs_dic, corr_noise_bands, rho
