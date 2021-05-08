import numpy as np

def get_exp_specs(expname, remove_atm = False):

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

    elif expname == 's4satEE':
        specs_dic = {
        #Table 2.1 of PBDR. Note that we are using noise in E maps (and not B maps) here for P.
        #freq: [beam_arcmins, white_noise_T, elknee_T, alphaknee_T, whitenoise_B, elknee_B, alphaknee_B] 
        30: [72.8, None, -1., 0., 3.74, 60., 2.2],
        40: [72.8, None, -1., 0., 4.73, 60., 2.2], 
        85: [25.5, None, -1., 0., 0.93, 60., 2.2],
        95: [22.7, None, -1., 0., 0.82, 60., 2.2],
        145: [25.5, None, -1., 0., 1.25, 65., 3.1],
        155: [22.7, None, -1., 0., 1.34, 65., 3.1],
        220: [13.0, None, -1., 0., 3.48, 65., 3.1],
        270: [13.0, None, -1., 0., 8.08, 65., 3.1],
        }

    elif expname == 's4satBB':
        specs_dic = {
        #Table 2.1 of PBDR. Note that we are using noise in B maps (and not E maps) here for P.
        #freq: [beam_arcmins, white_noise_T, elknee_T, alphaknee_T, whitenoise_B, elknee_B, alphaknee_B] 
        30: [72.8, None, -1., 0., 3.53, 60., 1.7],
        40: [72.8, None, -1., 0., 4.46, 60., 1.7], 
        85: [25.5, None, -1., 0., 0.88, 60., 1.7],
        95: [22.7, None, -1., 0., 0.78, 60., 1.7],
        145: [25.5, None, -1., 0., 1.23, 60., 3.],
        155: [22.7, None, -1., 0., 1.34, 60., 3.],
        220: [13.0, None, -1., 0., 3.48, 60., 3.],
        270: [13.0, None, -1., 0., 5.97, 60., 3.],
        }

    elif expname == 's4BB_DC06':
        specs_dic = {
        #Table 1 of http://bicep.rc.fas.harvard.edu/CMB-S4/analysis_logbook/20200208_06_sims_details/. Note that we are using noise in B maps (and not E maps) here for P.
        #freq: [beam_arcmins, white_noise_T, elknee_T, alphaknee_T, whitenoise_B, elknee_B, alphaknee_B] 
        20: [11.0, None, -1., 0., 1.69, 200., 1.5], #I think this is LAT
        30: [72.8, None, -1., 0., 1.4, 75., 1.5],
        40: [72.8, None, -1., 0., 1.76, 75., 1.5], 
        85: [25.5, None, -1., 0., 2.13, 75., 1.5],
        95: [22.7, None, -1., 0., 1.88, 75., 1.5],
        145: [25.5, None, -1., 0., 1.4, 60., 2.8],
        155: [22.7, None, -1., 0., 1.5, 60., 2.8],
        220: [13.0, None, -1., 0., 11.53, 60., 2.9],
        270: [13.0, None, -1., 0., 19.77, 60., 2.9],
        #now adding other LATS
        31: [7.3, None, -1., 0., 1.66, 20., 1.5],
        41: [5.5, None, -1., 0., 1.5, 200., 1.5], 
        96: [2.3, None, -1., 0., 0.226, 200., 1.5],
        146: [1.5, None, -1., 0., 0.32, 200., 2.8],
        221: [1.0, None, -1., 0., 1.9, 200., 2.9],
        271: [0.8, None, -1., 0., 3.27, 200., 2.9],
        }
    freqarr = sorted( specs_dic.keys() )
    nc = len( freqarr )

    if expname.find('s4sat')==-1:
        rho = 0.9
        corr_noise = 1
        if corr_noise:
            corr_noise_bands = {20: [20], 27:[39], 39:[27], 93:[145], 145:[93], 225: [278], 278: [225], 350: [350]}
        else:
            corr_noise_bands = {20: [20], 27:[39], 39:[27], 93:[93], 145:[145], 225: [225], 278: [278], 350: [350]}
    else:
        rho = 0.
        corr_noise = 0
        corr_noise_bands = {}
        for freq in freqarr:
            corr_noise_bands[freq] = [freq]

    return specs_dic, corr_noise_bands, rho, corr_noise
