import numpy as np, os, sys
import foregrounds as fg

def get_exp_specs(expname, remove_atm = 0):

    if expname.find('s4')>-1:
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

        elif expname == 's4deepv3r025' or expname == 's4deepv3r025plusspt4HF':
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

            if expname.find('spt4')>-1:
                spt4_specs_dic = get_spt4_specs(expname, which_spt4 = 'spt4_C3')
                for nu in spt4_specs_dic:
                    if nu in specs_dic:
                        nu_mod = nu + 1
                    else:
                        nu_mod = nu
                    specs_dic[nu_mod] = spt4_specs_dic[nu]

        freqarr = sorted( specs_dic.keys() )
        nc = len( freqarr )

        corr_noise = 1
        if corr_noise:
            corr_noise_bands = {20: [20], 27:[39], 39:[27], 93:[145], 145:[93], 225: [278], 278: [225], 350: [350]}
        else:
            corr_noise_bands = {20: [20], 27:[39], 39:[27], 93:[93], 145:[145], 225: [225], 278: [278], 350: [350]}
        rho = 0.9
        cib_corr_coeffs = None

        for freq in freqarr:
            if freq not in sorted(corr_noise_bands):
                corr_noise_bands[freq] = [freq]

    elif expname.find('actd56')>-1:

        if not remove_atm:
            elknee_T_90, elknee_T_150 = 2000., 4000
            alphaknee_T_90, alphaknee_T_150 = 3.5, 3.5
            elknee_P_90, elknee_P_150 = 700, 700
            alphaknee_P_90, alphaknee_P_150 = 1.4, 1.4
        else:
            elknee_T_90, elknee_T_150 = -1., -1.
            alphaknee_T_90, alphaknee_T_150 = 0., 0.
            elknee_P_90, elknee_P_150 = -1., -1.
            alphaknee_P_90, alphaknee_P_150 = 0., 0.

        beam_90, beam_150 = 2.2, 1.4
        white_noise_T_90, white_noise_T_150 = 19., 12.
        white_noise_P_90 = white_noise_T_90 * np.sqrt(2.)
        white_noise_P_150 = white_noise_T_150 * np.sqrt(2.)

        specs_dic = {
        90: [beam_90, white_noise_T_90, elknee_T_90, alphaknee_T_90, white_noise_P_90, elknee_P_90, alphaknee_P_90],
        150: [beam_150, white_noise_T_150, elknee_T_150, alphaknee_T_150, white_noise_P_150, elknee_P_150, alphaknee_P_150],
        }
        corr_noise = 0
        corr_noise_bands = {90:[90], 150:[150]}
        rho = 1.
        cib_corr_coeffs = None

    elif expname.find('actBN')>-1:

        if not remove_atm:
            elknee_T_90, elknee_T_150 = 2000., 4000
            alphaknee_T_90, alphaknee_T_150 = 3.5, 3.5
            elknee_P_90, elknee_P_150 = 700, 700
            alphaknee_P_90, alphaknee_P_150 = 1.4, 1.4
        else:
            elknee_T_90, elknee_T_150 = -1., -1.
            alphaknee_T_90, alphaknee_T_150 = 0., 0.
            elknee_P_90, elknee_P_150 = -1., -1.
            alphaknee_P_90, alphaknee_P_150 = 0., 0.

        beam_90, beam_150 = 2.2, 1.4
        white_noise_T_90, white_noise_T_150 = 35., 35.
        white_noise_P_90 = white_noise_T_90 * np.sqrt(2.)
        white_noise_P_150 = white_noise_T_150 * np.sqrt(2.)

        specs_dic = {
        90: [beam_90, white_noise_T_90, elknee_T_90, alphaknee_T_90, white_noise_P_90, elknee_P_90, alphaknee_P_90],
        150: [beam_150, white_noise_T_150, elknee_T_150, alphaknee_T_150, white_noise_P_150, elknee_P_150, alphaknee_P_150],
        }
        corr_noise = 0
        corr_noise_bands = {90:[90], 150:[150]}
        rho = 1.
        cib_corr_coeffs = None

    elif expname.find('advactC14')>-1:

        if not remove_atm:
            elknee_T_90, elknee_T_150, elknee_T_220 = 2000., 4000, 6740
            alphaknee_T_90, alphaknee_T_150, alphaknee_T_220 = 3.5, 3.5, 3.5
            elknee_P_90, elknee_P_150, elknee_P_220 = 700, 700, 700
            alphaknee_P_90, alphaknee_P_150, alphaknee_P_220 = 1.4, 1.4, 1.4
        else:
            elknee_T_90, elknee_T_150, elknee_T_220 = -1., -1., -1.
            alphaknee_T_90, alphaknee_T_150, alphaknee_T_220 = 0., 0., 0.
            elknee_P_90, elknee_P_150, elknee_P_220 = -1., -1., -1.
            alphaknee_P_90, alphaknee_P_150, alphaknee_P_220 = 0., 0., 0.
        '''

        #check CCAT-prime below for how Igot these numbers
        elknee_t_arr = [1000*2., 1000.*3, 1000.*6]
        alpha_knee_t_arr = [3.5, 3.5, 3.5]

        elknee_p_arr = [700., 700., 700.]
        alpha_knee_p_arr = [1.4, 1.4, 1.4]

        elknee_T_90, elknee_T_150, elknee_T_220 = elknee_t_arr
        alphaknee_T_90, alphaknee_T_150, alphaknee_T_220 = alpha_knee_t_arr
        elknee_P_90, elknee_P_150, elknee_P_220 = alpha_knee_t_arr
        alphaknee_P_90, alphaknee_P_150, alphaknee_P_220 = alpha_knee_p_arr
        '''


        #table 1 of https://arxiv.org/pdf/1406.4794.pdf
        beam_90, beam_150, beam_220 = 2.2, 1.3, 0.9
        white_noise_T_90, white_noise_T_150, white_noise_T_220 = 7.8, 6.9, 25.
        white_noise_P_90, white_noise_P_150, white_noise_P_220 = white_noise_T_90 * np.sqrt(2.), white_noise_T_150 * np.sqrt(2.), white_noise_T_220 * np.sqrt(2.)

        specs_dic = {
        90: [beam_90, white_noise_T_90, elknee_T_90, alphaknee_T_90, white_noise_P_90, elknee_P_90, alphaknee_P_90],
        150: [beam_150, white_noise_T_150, elknee_T_150, alphaknee_T_150, white_noise_P_150, elknee_P_150, alphaknee_P_150],
        220: [beam_220, white_noise_T_220, elknee_T_220, alphaknee_T_220, white_noise_P_220, elknee_P_220, alphaknee_P_220],
        }
        corr_noise = 0
        corr_noise_bands = {90:[90], 150:[150], 220: [220]}
        rho = 1.
        cib_corr_coeffs = None

    elif expname.find('ccat_prime_so')>-1:

        #SO and CCAT-prime combinations

        noise_220 = (1./22.**2. + 1./14.5**2.)**-0.5
        noise_278 = (1./54.**2. + 1./27.5**2.)**-0.5
        noise_350, noise_410 = 105., 376.5

        #SO/CCAT_prime beams. HF (>220 GHz) beams are CCAT-prime
        beam_size_arr = [7.4, 5.1, 2.2, 1.4, 1.0, 45./60., 35./60., 30./60.]
        noiselevel_arr = [52., 27., 5.8, 6.3, noise_220, noise_278, noise_350, noise_410]
        nu_arr = [30, 39, 93, 145, 220, 278, 353, 410]


        '''
        if (0): #this is how I modified SO Nwhite/Nred into normal nl definitons.

            elknee_t_93, alpha_knee_93 = 1000., 3.5
            elknee_t_145, alpha_knee_145 = 1000., 3.5
            elknee_t_220, alpha_knee_220 = 1000., 3.5
            elknee_t_278, alpha_knee_278 = 1000., 3.5
            elknee_t_353, alpha_knee_353 = 1000., 3.5
            elknee_t_410, alpha_knee_410 = 1000., 3.5
            Nred_93, Nred_145, Nred_220, Nred_278 = 230., 1450., 17000., 31000. #uK^2 seconds

            NTube_years_MF = 5.
            f_sky = 0.35
            survey_time = 1. #years--- given we are using "tube-years above, this should not be changed.
            t = survey_time * 365.25 * 24. * 3600.    ## convert years to seconds
            #t = t * 0.2 *0.5  ## retention after observing efficiency and cuts, daytime only
            t = t * 0.2  ## retention after observing efficiency and cuts, both day and night as i cannot find this in SO forecast
            t = t * 0.85  ## a kludge for the noise non-uniformity of the map edges
            A_SR = 4. * np.pi * f_sky  ## sky areas in Steradians
            A_deg =  A_SR * (180/np.pi)**2  ## sky area in square degrees
            A_arcmin = A_deg * 3600.

            Nred_93  = Nred_93  * A_SR / t / NTube_years_MF
            Nred_145  = Nred_145  * A_SR / t / NTube_years_MF
            Nred_220  = Nred_220  * A_SR / t / NTube_years_MF
            Nred_278  = Nred_278  * A_SR / t / NTube_years_MF

            #CCAT-prime from table VII of https://arxiv.org/pdf/2008.11688.pdf
            Nred_220, Nred_278 = 1.6e-2, 1.1e-1
            Nred_353, Nred_410 = 2.7, 17.

            w93, w145, w220, w278, w353, w410 = noiselevel_arr
            ell = np.arange(10000)

            nl_93 = (w93*np.pi/180./60.)**2 + Nred_93 * (ell/elknee_t_93)**-alpha_knee_93
            nl_93_v2 = (w93*np.pi/180./60.)**2 * ( 1. + (ell/2000.)**-alpha_knee_93 )
            nl_145 = (w145*np.pi/180./60.)**2 + Nred_145 * (ell/elknee_t_145)**-alpha_knee_145
            nl_145_v2 = (w145*np.pi/180./60.)**2 * ( 1. + (ell/elknee_t_145/3.)**-alpha_knee_145 )
            nl_220 = (w220*np.pi/180./60.)**2 + Nred_220 * (ell/elknee_t_220)**-alpha_knee_220
            nl_220_v2 = (w220*np.pi/180./60.)**2 * ( 1. + (ell/elknee_t_220/6.)**-alpha_knee_220 )
            nl_278 = (w278*np.pi/180./60.)**2 + Nred_278 * (ell/elknee_t_278)**-alpha_knee_278
            nl_278_v2 = (w278*np.pi/180./60.)**2 * (1. + (ell/elknee_t_278/8.)**-alpha_knee_278 )
            nl_353 = (w353*np.pi/180./60.)**2 + Nred_353 * (ell/elknee_t_353)**-alpha_knee_353
            nl_353_v2 = (w353*np.pi/180./60.)**2 * (1. + (ell/elknee_t_353/8.)**-alpha_knee_353 )
            nl_410 = (w410*np.pi/180./60.)**2 + Nred_410 * (ell/elknee_t_410)**-alpha_knee_410
            nl_410_v2 = (w410*np.pi/180./60.)**2 * (1. + (ell/elknee_t_410/8.)**-alpha_knee_410 )

            ax = subplot(111, yscale = 'log')
            plot(nl_39)
            plot(nl_39_v2)
            show()

        elknee_t_arr = [1000./2., 1000./2., 1000*2., 1000.*3, 1000.*6, 1000.*8, 1000.*8, 1000.*8]
        alpha_knee_t_arr = [3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5]

        elknee_p_arr = [700., 700., 700., 700., 700., 700., 700., 700.]
        alpha_knee_p_arr = [1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4]
        '''
        if not remove_atm:
            elknee_T_27, elknee_T_39, elknee_T_90, elknee_T_150, elknee_T_225, elknee_T_278, elknee_T_353, elknee_T_410 = 415., 391., 1932., 3917., 6740., 6792., 6792., 6792.
            alphaknee_T_27, alphaknee_T_39, alphaknee_T_90, alphaknee_T_150, alphaknee_T_225, alphaknee_T_278, alphaknee_T_353, alphaknee_T_410 = 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5
            elknee_P_27, elknee_P_39, elknee_P_90, elknee_P_150, elknee_P_225, elknee_P_278, elknee_P_353, elknee_P_410 = 700., 700., 700., 700., 700., 700., 700., 700.
            alphaknee_P_27, alphaknee_P_39, alphaknee_P_90, alphaknee_P_150, alphaknee_P_225, alphaknee_P_278, alphaknee_P_353, alphaknee_P_410 = 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4
        else:
            elknee_T_27, elknee_T_39, elknee_T_90, elknee_T_150, elknee_T_225, elknee_T_278, elknee_T_353, elknee_T_410 = -1., -1., -1., -1., -1., -1., -1., -1.
            alphaknee_T_27, alphaknee_T_39, alphaknee_T_90, alphaknee_T_150, alphaknee_T_225, alphaknee_T_278, alphaknee_T_353, alphaknee_T_410 = 0.0, 0., 0., 0., 0., 0., 0., 0.
            elknee_P_27, elknee_P_39, elknee_P_90, elknee_P_150, elknee_P_225, elknee_P_278, elknee_P_353, elknee_P_410 = -1., -1., -1., -1., -1., -1., -1., -1.
            alphaknee_P_27, alphaknee_P_39, alphaknee_P_90, alphaknee_P_150, alphaknee_P_225, alphaknee_P_278, alphaknee_P_353, alphaknee_P_410 = 0.0, 0., 0., 0., 0., 0., 0., 0.

        elknee_t_arr = [elknee_T_27, elknee_T_39, elknee_T_90, elknee_T_150, elknee_T_225, elknee_T_278, elknee_T_353, elknee_T_410]
        alpha_knee_t_arr = [alphaknee_T_27, alphaknee_T_39, alphaknee_T_90, alphaknee_T_150, alphaknee_T_225, alphaknee_T_278, alphaknee_T_353, alphaknee_T_410]
        elknee_p_arr = [elknee_P_27, elknee_P_39, elknee_P_90, elknee_P_150, elknee_P_225, elknee_P_278, elknee_P_353, elknee_P_410 ]
        alpha_knee_p_arr = [alphaknee_P_27, alphaknee_P_39, alphaknee_P_90, alphaknee_P_150, alphaknee_P_225, alphaknee_P_278, alphaknee_P_353, alphaknee_P_410]

        specs_dic = {}
        for nucntr, nu in enumerate( nu_arr ):
            specs_dic[nu] = [beam_size_arr[nucntr], noiselevel_arr[nucntr], elknee_t_arr[nucntr], alpha_knee_t_arr[nucntr], noiselevel_arr[nucntr]*1.414, elknee_p_arr[nucntr], alpha_knee_p_arr[nucntr]]

        freqarr = sorted( specs_dic.keys() )
        nc = len( freqarr )

        corr_noise = 1
        corr_noise_bands = {30:[39], 39:[30], 93:[145], 145:[93], 220: [278], 278: [220], 353: [410], 410: [353]}
        rho = 0.9
        cib_corr_coeffs = None

    elif expname.lower().find('cmb_bharat')>-1:

        draft_code_path = '/Users/sraghunathan/Research/SPTPol/analysis/git/DRAFT/'
        if not os.path.exists(draft_code_path):
            draft_code_path = '/data48/sri/git/DRAFT/'

        cmb_bharat_fname = '%s/data/cmb_bharat_config.txt' %(draft_code_path)
        cmb_bharat_specs = np.loadtxt(cmb_bharat_fname, usecols = [0, 1, 3, 4], delimiter = '&')

        v1_v2_line_splitter = 22
        if expname.lower() == 'cmb_bharat_1':
            cmb_bharat_specs = cmb_bharat_specs[: v1_v2_line_splitter]
        elif expname.lower() == 'cmb_bharat_2':
            cmb_bharat_specs = cmb_bharat_specs[v1_v2_line_splitter: ]

        freq_arr, beam_arr, white_noise_T_arr, white_noise_P_arr = cmb_bharat_specs.T

        specs_dic = {}
        corr_noise = 0
        corr_noise_bands = {}
        rho = 1.
        for freq, beam, white_noise_T, white_noise_P in zip( freq_arr, beam_arr, white_noise_T_arr, white_noise_P_arr):
            elknee_T, alphaknee_T, elknee_P, alphaknee_P = -1., 0., -1., 0.
            specs_dic[freq] = [beam, white_noise_T, elknee_T, alphaknee_T, white_noise_P, elknee_P, alphaknee_P]
            corr_noise_bands[freq] = [freq]

        cib_corr_coeffs = fg.get_cib_decorrelations_via_intrp(freq_arr)

    elif expname.find('so')>-1:

        if not remove_atm:
            elknee_T_27, elknee_T_39, elknee_T_90, elknee_T_150, elknee_T_225, elknee_T_278 = 415., 391., 1932., 3917., 6740., 6792.
            alphaknee_T_27, alphaknee_T_39, alphaknee_T_90, alphaknee_T_150, alphaknee_T_225, alphaknee_T_278 = 3.5, 3.5, 3.5, 3.5, 3.5, 3.5
            elknee_P_27, elknee_P_39, elknee_P_90, elknee_P_150, elknee_P_225, elknee_P_278 = 700., 700., 700., 700., 700., 700.
            alphaknee_P_27, alphaknee_P_39, alphaknee_P_90, alphaknee_P_150, alphaknee_P_225, alphaknee_P_278 = 1.4, 1.4, 1.4, 1.4, 1.4, 1.4
        else:
            elknee_T_27, elknee_T_39, elknee_T_90, elknee_T_150, elknee_T_225, elknee_T_278 = -1., -1., -1., -1., -1., -1.
            alphaknee_T_27, alphaknee_T_39, alphaknee_T_90, alphaknee_T_150, alphaknee_T_225, alphaknee_T_278 = 0.0, 0., 0., 0., 0., 0.
            elknee_P_27, elknee_P_39, elknee_P_90, elknee_P_150, elknee_P_225, elknee_P_278 = -1., -1., -1., -1., -1., -1.
            alphaknee_P_27, alphaknee_P_39, alphaknee_P_90, alphaknee_P_150, alphaknee_P_225, alphaknee_P_278 = 0.0, 0., 0., 0., 0., 0.
        '''
        #check CCAT-prime below for how Igot these numbers
        elknee_t_arr = [1000./2., 1000./2., 1000*2., 1000.*3, 1000.*6, 1000.*8]
        alpha_knee_t_arr = [3.5, 3.5, 3.5, 3.5, 3.5, 3.5]

        elknee_p_arr = [700., 700., 700., 700., 700., 700.]
        alpha_knee_p_arr = [1.4, 1.4, 1.4, 1.4, 1.4, 1.4]

        elknee_T_27, elknee_T_39, elknee_T_90, elknee_T_150, elknee_T_225, elknee_T_278 = elknee_t_arr
        alphaknee_T_27, alphaknee_T_39, alphaknee_T_90, alphaknee_T_150, alphaknee_T_225, alphaknee_T_278 = alpha_knee_t_arr
        elknee_P_27, elknee_P_39, elknee_P_90, elknee_P_150, elknee_P_225, elknee_P_278 = alpha_knee_t_arr
        alphaknee_P_27, alphaknee_P_39, alphaknee_P_90, alphaknee_P_150, alphaknee_P_225, alphaknee_P_278 = alpha_knee_p_arr
        '''

        if expname == 'sobaseline':
            white_noise_T_27, white_noise_T_39 = 71., 36.
            white_noise_T_90, white_noise_T_150 = 8., 10.
            white_noise_T_225, white_noise_T_278 = 22., 54.
        elif expname == 'sogoal':
            white_noise_T_27, white_noise_T_39 = 52., 27.
            white_noise_T_90, white_noise_T_150 = 5.8, 6.3
            white_noise_T_225, white_noise_T_278 = 15., 37.

        white_noise_P_27, white_noise_P_39 = white_noise_T_27* np.sqrt(2.), white_noise_T_39 * np.sqrt(2.)
        white_noise_P_90, white_noise_P_150 = white_noise_T_90* np.sqrt(2.), white_noise_T_150 * np.sqrt(2.)
        white_noise_P_225, white_noise_P_278 = white_noise_T_225* np.sqrt(2.), white_noise_T_278 * np.sqrt(2.)

        beam_27, beam_39, beam_90, beam_150, beam_225, beam_278 = 7.4, 5.1, 2.2, 1.4, 1.0, 0.9

        specs_dic = {
        27: [beam_27, white_noise_T_27, elknee_T_27, alphaknee_T_27, white_noise_P_27, elknee_P_27, alphaknee_P_27],
        39: [beam_39, white_noise_T_39, elknee_T_39, alphaknee_T_39, white_noise_P_39, elknee_P_39, alphaknee_P_39],
        93: [beam_90, white_noise_T_90, elknee_T_90, alphaknee_T_90, white_noise_P_90, elknee_P_90, alphaknee_P_90],
        145: [beam_150, white_noise_T_150, elknee_T_150, alphaknee_T_150, white_noise_P_150, elknee_P_150, alphaknee_P_150],
        225: [beam_225, white_noise_T_225, elknee_T_225, alphaknee_T_225, white_noise_P_225, elknee_P_225, alphaknee_P_225],
        278: [beam_278, white_noise_T_278, elknee_T_278, alphaknee_T_278, white_noise_P_278, elknee_P_278, alphaknee_P_278],
        }

        freqarr = sorted( specs_dic.keys() )
        nc = len( freqarr )

        corr_noise = 1
        if corr_noise:
            corr_noise_bands = {27:[39], 39:[27], 93:[145], 145:[93], 225: [278], 278: [225]}
        else:
            corr_noise_bands = {27:[39], 39:[27], 93:[93], 145:[145], 225: [225], 278: [278]}
        rho = 0.9
        cib_corr_coeffs = None

    elif expname.find('spt')>-1:

        elknee_T_90, elknee_T_150, elknee_T_220 = 1200., 2200., 2300.
        alphaknee_T_90, alphaknee_T_150, alphaknee_T_220 = 3., 4., 4.
        elknee_P_90, elknee_P_150, elknee_P_220 = 300., 300., 300.
        alphaknee_P_90, alphaknee_P_150, alphaknee_P_220 = 1., 1., 1.

        if remove_atm:
            elknee_T_90, elknee_T_150, elknee_T_220 = -1., -1., -1.
            alphaknee_T_90, alphaknee_T_150, alphaknee_T_220 = 0., 0., 0.
            elknee_P_90, elknee_P_150, elknee_P_220 = -1., -1., -1.
            alphaknee_P_90, alphaknee_P_150, alphaknee_P_220 = 0., 0., 0.

        beam_90, beam_150, beam_220 = 1.7, 1.2, 1.0
        white_noise_T_90, white_noise_T_150, white_noise_T_220 = 1e6, 1e6, 1e6
        white_noise_P_90, white_noise_P_150, white_noise_P_220 = 1e6, 1e6, 1e6

        if expname == 'sptsz':
            #slide 3 of Lindsey https://pole.uchicago.edu/sptpol/images/Spt_des_clusters_july17_f2f.pdf
            #white_noise_T_90, white_noise_T_150, white_noise_T_220 = 40., 17., 80.

            white_noise_T_90, white_noise_T_150, white_noise_T_220 = 35., 16., 61. # from TC slack spt_ymap chat on 4 Sep 2020.
            white_noise_P_90, white_noise_P_150, white_noise_P_220 = 1e6, 1e6, 1e6
        elif expname == 'sptszplanck':
            white_noise_T_90, white_noise_T_150, white_noise_T_220 = 35., 16., 61. # from TC slack spt_ymap chat on 4 Sep 2020.
            white_noise_P_90, white_noise_P_150, white_noise_P_220 = 1e6, 1e6, 1e6
        elif expname == 'sptpolultradeep':
            #slide 3 of Lindsey https://pole.uchicago.edu/sptpol/images/Spt_des_clusters_july17_f2f.pdf
            #white_noise_T_90, white_noise_T_150, white_noise_T_220 = 10., 3.8, 1e6
            white_noise_T_90, white_noise_T_150, white_noise_T_220 = 10., 6.0, 1e6

            white_noise_P_90 = white_noise_T_90 * np.sqrt(2.)
            white_noise_P_150 = white_noise_T_150 * np.sqrt(2.)
            white_noise_P_220 = white_noise_T_220 * np.sqrt(2.)
        elif expname == 'sptpolultradeepplus3g':
            #slide 3 of Lindsey https://pole.uchicago.edu/sptpol/images/Spt_des_clusters_july17_f2f.pdf
            #white_noise_T_90, white_noise_T_150, white_noise_T_220 = 10., 3.8, 30.
            white_noise_T_90, white_noise_T_150, white_noise_T_220 = 10., 6.0, 30.

            white_noise_P_90 = white_noise_T_90 * np.sqrt(2.)
            white_noise_P_150 = white_noise_T_150 * np.sqrt(2.)
            white_noise_P_220 = white_noise_T_220 * np.sqrt(2.)
        elif expname == 'sptpolplusultradeep':
            #slide 3 of Lindsey https://pole.uchicago.edu/sptpol/images/Spt_des_clusters_july17_f2f.pdf
            white_noise_T_90, white_noise_T_220 = 10., 1e6
            #white_noise_T_150 = ( (1./3.8**2.) + (1./5**2.) )**-0.5
            white_noise_T_150 = ( (1./6.**2.) + (1./6.**2.) )**-0.5

            white_noise_P_90 = white_noise_T_90 * np.sqrt(2.)
            white_noise_P_150 = white_noise_T_150 * np.sqrt(2.)
            white_noise_P_220 = white_noise_T_220 * np.sqrt(2.)
        elif expname == 'sptpolplusultradeepplus3g':
            #slide 3 of Lindsey https://pole.uchicago.edu/sptpol/images/Spt_des_clusters_july17_f2f.pdf
            #white_noise_T_90, white_noise_T_220 = 10., 30.
            #white_noise_T_150 = ( (1./3.8**2.) + (1./5**2.) )**-0.5

            """
            #3G noise levels come from Wei's presentation during SPTf2f 20200707. slide 11 from the below link.
            #https://pole.uchicago.edu/spt3g/images/20200707_Face-to-Face_Meeting_White_Noise.pdf
            #100d noise levels come from this: https://pole.uchicago.edu/sptpol/images/Spt_des_clusters_july17_f2f.pdf
            """
            white_noise_T_90 = ( (1./10.**2.) + (1./12.**2.) + (1./7.2**2.) )**-0.5 #100d + 500d + 3G
            #white_noise_T_150 = ( (1./3.8**2.) + (1./5.**2.) + (1./5.8**2.) )**-0.5 #100d + 500d + 3G
            white_noise_T_150 = ( (1./6.**2.) + (1./6.**2.) + (1./5.8**2.) )**-0.5 #100d + 500d + 3G
            white_noise_T_220 = 21. #3G-only


            white_noise_P_90 = white_noise_T_90 * np.sqrt(2.)
            white_noise_P_150 = white_noise_T_150 * np.sqrt(2.)
            white_noise_P_220 = white_noise_T_220 * np.sqrt(2.)
        elif expname == 'sptpolplusultradeepplus3gfull':
            #slide 3 of Lindsey https://pole.uchicago.edu/sptpol/images/Spt_des_clusters_july17_f2f.pdf
            #white_noise_T_90, white_noise_T_220 = 10., 10.
            #white_noise_T_150 = ( (1./3.8**2.) + (1./5**2.) )**-0.5

            """
            #100d + 500d + full 3g_depth
            #100d noise levels come from this: https://pole.uchicago.edu/sptpol/images/Spt_des_clusters_july17_f2f.pdf
            """
            white_noise_T_90 = ( (1./10.**2.) + (1./12.**2.) + (1./3.0**2.) )**-0.5 #100d + 500d + 3G
            #white_noise_T_150 = ( (1./3.8**2.) + (1./5.**2.) + (1./2.2**2.) )**-0.5 #100d + 500d + 3G
            white_noise_T_150 = ( (1./6.**2.) + (1./6.**2.) + (1./2.2**2.) )**-0.5 #100d + 500d + 3G
            white_noise_T_220 = 8. #3G-only            

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

        elif expname == 'spt3g' or expname == 'spt3gplusherschel' or expname == 'spt3g_plusplanck' or expname == 'spt3gplusherschel_plusplanck' or expname == 'spt3gplancksevenbands':
            white_noise_T_90 = 3.0
            white_noise_T_150 = 2.2
            white_noise_T_220 = 8.8

            white_noise_P_90 = np.sqrt(2.) * white_noise_T_90
            white_noise_P_150 = np.sqrt(2.) * white_noise_T_150
            white_noise_P_220 = np.sqrt(2.) * white_noise_T_220

        ######################################################################
        ######################################################################
        ######################################################################
        elif expname == 'spt3g_winter_2020':

            white_noise_T_90 = 5.0
            white_noise_T_150 = 4.0
            white_noise_T_220 = 14.0

            white_noise_P_90, white_noise_P_150, white_noise_P_220 = np.sqrt(2.) * white_noise_T_90, np.sqrt(2.) * white_noise_T_150, np.sqrt(2.) * white_noise_T_220

        elif expname == 'spt3g_summer_el1c_el2c_2020':
            #https://docs.google.com/spreadsheets/d/1GwF8kJjl3czbJcy68ft-LpbStAWtkq3RPNX5aA1MaAY/edit#gid=0
            noise_90_arr = [10.9, 13.1]
            noise_150_arr = [10.8, 12.]
            noise_220_arr = [38.9, 44.2]

            white_noise_T_90 = np.mean(noise_90_arr)
            white_noise_T_150 = np.mean(noise_150_arr)
            white_noise_T_220 = np.mean(noise_220_arr)

            white_noise_P_90, white_noise_P_150, white_noise_P_220 = np.sqrt(2.) * white_noise_T_90, np.sqrt(2.) * white_noise_T_150, np.sqrt(2.) * white_noise_T_220

        elif expname == 'spt3g_summer_el1b_el2b_2020':
            #https://docs.google.com/spreadsheets/d/1GwF8kJjl3czbJcy68ft-LpbStAWtkq3RPNX5aA1MaAY/edit#gid=0
            noise_90_arr = [10.9, 13.1]
            noise_150_arr = [10.8, 12.]
            noise_220_arr = [38.9, 44.2]

            white_noise_T_90 = np.mean(noise_90_arr)
            white_noise_T_150 = np.mean(noise_150_arr)
            white_noise_T_220 = np.mean(noise_220_arr)

            white_noise_P_90, white_noise_P_150, white_noise_P_220 = np.sqrt(2.) * white_noise_T_90, np.sqrt(2.) * white_noise_T_150, np.sqrt(2.) * white_noise_T_220

        elif expname == 'spt3g_summer_el1_e5_2020':
            #https://docs.google.com/spreadsheets/d/1GwF8kJjl3czbJcy68ft-LpbStAWtkq3RPNX5aA1MaAY/edit#gid=0
            noise_90_arr = [10.4, 10.6, 10.4, 10.6, 10.5]
            noise_150_arr = [10.3, 9.7, 9.5, 9.6, 9.5]
            noise_220_arr = [37.1, 35.7, 34.5, 35.1, 34.7]

            white_noise_T_90 = np.mean(noise_90_arr)
            white_noise_T_150 = np.mean(noise_150_arr)
            white_noise_T_220 = np.mean(noise_220_arr)

            white_noise_P_90, white_noise_P_150, white_noise_P_220 = np.sqrt(2.) * white_noise_T_90, np.sqrt(2.) * white_noise_T_150, np.sqrt(2.) * white_noise_T_220

        ######################################################################
        ######################################################################
        ######################################################################

        elif expname == 'spt3g20192020':
            white_noise_T_90 = 6.
            white_noise_T_150 = 5.
            white_noise_T_220 = 17.

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

        elif expname == 'spt3g_201920' or expname == 'spt3g_201920_plusplanck' or expname.lower() == 'spt3g_201920_plancksevenbands':
            white_noise_T_90 = 6.2
            white_noise_T_150 = 4.6
            white_noise_T_220 = 16.

            white_noise_P_90 = np.sqrt(2.) * white_noise_T_90
            white_noise_P_150 = np.sqrt(2.) * white_noise_T_150
            white_noise_P_220 = np.sqrt(2.) * white_noise_T_220

        elif expname == 'sptpolplusultradeepplus3gplusherschel' or expname == 'sptpolplusultradeepplus3gplusherschel_202009'\
            or expname == 'sptpolplusultradeepplus3gplusherschel_plusplanck' or expname == 'sptpolplusultradeepplus3gplusherschel_202009_plusplanck'\
            or expname == 'sptpolplusultradeepplus3gplusherschelplus353planck_202009'\
            or expname == 'sptpolplusultradeepplus3gplusherschelplus353planck_202009_remove3g150'\
            or expname == 'sptpolplusultradeepplus3gplusherschel_202009_highnoiseherschel'\
            or expname.lower() == 'allspt_herschel_planck':
            #slide 3 of Lindsey https://pole.uchicago.edu/sptpol/images/Spt_des_clusters_july17_f2f.pdf
            #white_noise_T_90, white_noise_T_220 = 10., 30.

            #white_noise_T_90 = ( (1./10.**2.) + (1./12.**2.) + (1./10.2**2.) )**-0.5
            #white_noise_T_150 = ( (1./3.8**2.) + (1./5.**2.) + (1./8.2**2.) )**-0.5
            #white_noise_T_220 = 30. ##8.

            """
            #3G noise levels come from Wei's presentation during SPTf2f 20200707. slide 11 from the below link.
            #https://pole.uchicago.edu/spt3g/images/20200707_Face-to-Face_Meeting_White_Noise.pdf
            #100d noise levels come from this: https://pole.uchicago.edu/sptpol/images/Spt_des_clusters_july17_f2f.pdf
            """
            white_noise_T_90 = ( (1./10.**2.) + (1./12.**2.) + (1./7.2**2.) )**-0.5 #100d + 500d + 3G
            #white_noise_T_150 = ( (1./3.8**2.) + (1./5.**2.) + (1./5.8**2.) )**-0.5 #100d + 500d + 3G
            white_noise_T_150 = ( (1./6.**2.) + (1./6.**2.) + (1./5.8**2.) )**-0.5 #100d + 500d + 3G
            if expname.lower().find('remove3g150')>-1: #removing 3G 150
                white_noise_T_150 = ( (1./6.**2.) + (1./6.**2.) )**-0.5 #100d + 500d
                #white_noise_T_150 = 6.

            if expname.lower().find('202009')>-1 or expname.lower() == 'allspt_herschel_planck':
                white_noise_T_220 = 17. #3G-only
                ###white_noise_T_220 = 1e6 #3G-only
            else:
                white_noise_T_220 = 21. #3G-only

            if expname.find('plus353planck')>-1:
                beam_353 = 4.818
                white_noise_T_353 = 150.
                white_noise_P_353 = 1e6
                elknee_T_353 = elknee_P_353 = -1.
                alphaknee_T_353 = alphaknee_P_353 = 0.

            #white_noise_T_600, white_noise_T_857, white_noise_T_1200 = 1., 1., 1.
            #20200707
            ###white_noise_T_600, white_noise_T_857, white_noise_T_1200 = 2.1e4/2., 3.2e5/2., 2.5e7/2. #Based on Gil's 2013 paper
            white_noise_T_600, white_noise_T_857, white_noise_T_1200 = 2.07769677e+03, 3.29858731e+04, 2.51686637e+06 #based on SPTxSPIRE paper

            if expname.find('highnoiseherschel')>-1:
                white_noise_T_600 = white_noise_T_600 * 2.
                white_noise_T_857 = white_noise_T_857 * 2.
                white_noise_T_1200 = white_noise_T_1200 * 2.
        
            white_noise_P_90 = white_noise_T_90 * np.sqrt(2.)
            white_noise_P_150 = white_noise_T_150 * np.sqrt(2.)
            white_noise_P_220 = white_noise_T_220 * np.sqrt(2.)

            white_noise_P_600 = 1e6
            white_noise_P_857 = 1e6
            white_noise_P_1200 = 1e6

            beam_600, beam_857, beam_1200 = 36.6/60., 25.2/60., 18.1/60.
            elknee_T_600 = elknee_T_857 = elknee_T_1200 = -1. ##2500
            alphaknee_T_600 = alphaknee_T_857 = alphaknee_T_1200 = 0. ##4.

            elknee_P_600 = elknee_P_857 = elknee_P_1200 = -1. ##300.
            alphaknee_P_600 = alphaknee_P_857 = alphaknee_P_1200 = 0. ##1.

        elif expname == 'sptpolplusultradeepplus3gplusherschel_v2' or expname == 'sptpolplusultradeepplus3gplusherschel_v2_plusplanck':

            #same as above but does not use SPT-3G 220 GHz data
            white_noise_T_90 = ( (1./10.**2.) + (1./12.**2.) + (1./7.2**2.) )**-0.5 #100d + 500d + 3G
            #white_noise_T_150 = ( (1./3.8**2.) + (1./5.**2.) + (1./5.8**2.) )**-0.5 #100d + 500d + 3G
            white_noise_T_150 = ( (1./6.**2.) + (1./6.**2.) + (1./5.8**2.) )**-0.5 #100d + 500d + 3G
            white_noise_T_220 = 1e6
        
            white_noise_P_90 = white_noise_T_90 * np.sqrt(2.)
            white_noise_P_150 = white_noise_T_150 * np.sqrt(2.)
            white_noise_P_220 = white_noise_T_220 * np.sqrt(2.)

        elif expname == 'sptpolplusultradeepplus3gplusherschel_v3':

            #same as above but does not use SPTpol 90 GHz data
            white_noise_T_90 = 7.2 #3G-only
            #white_noise_T_150 = ( (1./3.8**2.) + (1./5.**2.) + (1./5.8**2.) )**-0.5 #100d + 500d + 3G
            white_noise_T_150 = ( (1./6.**2.) + (1./6.**2.) + (1./5.8**2.) )**-0.5 #100d + 500d + 3G
            white_noise_T_220 = 21. #3G-only
        
            white_noise_P_90 = white_noise_T_90 * np.sqrt(2.)
            white_noise_P_150 = white_noise_T_150 * np.sqrt(2.)
            white_noise_P_220 = white_noise_T_220 * np.sqrt(2.)

        elif expname == 'sptpolplusultradeepplus3gplusherschel_v4':

            #same as above but does not use SPTpol 90 GHz or SPT-3G 220 data
            white_noise_T_90 = 7.2 #3G-only
            #white_noise_T_150 = ( (1./3.8**2.) + (1./5.**2.) + (1./5.8**2.) )**-0.5 #100d + 500d + 3G
            white_noise_T_150 = ( (1./6.**2.) + (1./6.**2.) + (1./5.8**2.) )**-0.5 #100d + 500d + 3G
            white_noise_T_220 = 1e6
        
            white_noise_P_90 = white_noise_T_90 * np.sqrt(2.)
            white_noise_P_150 = white_noise_T_150 * np.sqrt(2.)
            white_noise_P_220 = white_noise_T_220 * np.sqrt(2.)

        #remove 1/f if Planck is included
        if expname.find('plusplanck')>-1:
            elknee_T_90, elknee_T_150, elknee_T_220 = -1., -1., -1.
            alphaknee_T_90, alphaknee_T_150, alphaknee_T_220 = 0., 0., 0.
            elknee_P_90, elknee_P_150, elknee_P_220 = -1., -1., -1.
            alphaknee_P_90, alphaknee_P_150, alphaknee_P_220 = 0., 0., 0.            

        specs_dic = {
        90: [beam_90, white_noise_T_90, elknee_T_90, alphaknee_T_90, white_noise_P_90, elknee_P_90, alphaknee_P_90],
        150: [beam_150, white_noise_T_150, elknee_T_150, alphaknee_T_150, white_noise_P_150, elknee_P_150, alphaknee_P_150],
        220: [beam_220, white_noise_T_220, elknee_T_220, alphaknee_T_220, white_noise_P_220, elknee_P_220, alphaknee_P_220],
        }
        corr_noise = 0
        corr_noise_bands = {90:[90], 150:[150], 220: [220]}
        rho = 1.

        cib_corr_coeffs = None

        if expname.lower() == 'spt3gplancksevenbands' or expname.lower() == 'spt3g_201920_plancksevenbands' or expname.lower() == 'allspt_herschel_planck':

            planck_nu_arr = [100, 143, 217, 353]
            planck_beam_arr = [9.651, 7.248, 4.990, 4.818]
            planck_white_noise_T_arr = [70., 35., 47., 150.]
            planck_white_noise_P_arr = None
            planck_elknee_T_arr = None
            planck_elknee_P_arr = None
            planck_alphaknee_T_arr = None
            planck_alphaknee_P_arr = None

            for nucntr, nu in enumerate( planck_nu_arr ):
                specs_dic[nu] = [planck_beam_arr[nucntr], planck_white_noise_T_arr[nucntr], -1., 0., 1e6, -1., 0.]
                corr_noise_bands[nu] = [nu]

        if expname.find('plus353planck')>-1:
            specs_dic[353] = [beam_353, white_noise_T_353, elknee_T_353, alphaknee_T_353, white_noise_P_353, elknee_P_353, alphaknee_P_353]
            corr_noise_bands[353] = [353]

        if expname == 'sptszplanck':
            white_noise_T_353 = 150.
            white_noise_P_353 = 1e6

            elknee_T_90, elknee_T_150, elknee_T_220 = -1., -1., -1.
            alphaknee_T_90, alphaknee_T_150, alphaknee_T_220 = 0., 0., 0.
            elknee_P_90, elknee_P_150, elknee_P_220 = -1., -1., -1.
            alphaknee_P_90, alphaknee_P_150, alphaknee_P_220 = 0., 0., 0.

            elknee_T_353 = -1.
            alphaknee_T_353 = 0.
            elknee_P_353 = -1.
            alphaknee_P_353 = 0.

            beam_353 = 4.818

            specs_dic = {
            90: [beam_90, white_noise_T_90, elknee_T_90, alphaknee_T_90, white_noise_P_90, elknee_P_90, alphaknee_P_90],
            150: [beam_150, white_noise_T_150, elknee_T_150, alphaknee_T_150, white_noise_P_150, elknee_P_150, alphaknee_P_150],
            220: [beam_220, white_noise_T_220, elknee_T_220, alphaknee_T_220, white_noise_P_220, elknee_P_220, alphaknee_P_220],
            353: [beam_353, white_noise_T_353, elknee_T_353, alphaknee_T_353, white_noise_P_353, elknee_P_353, alphaknee_P_353],
            }

            corr_noise = 0
            corr_noise_bands = {90:[90], 150:[150], 220: [220], 353:[353]}
            rho = 1.

        #### SPT4
        if expname.find('spt4')>-1:
            white_noise_T_90 = 3.0
            white_noise_T_150 = 2.2
            white_noise_T_220 = 8.8

            white_noise_P_90 = np.sqrt(2.) * white_noise_T_90
            white_noise_P_150 = np.sqrt(2.) * white_noise_T_150
            white_noise_P_220 = np.sqrt(2.) * white_noise_T_220

            spt4_white_noise_dic = {'spt4_C1': [4.1, 4.8, 28.4], 'spt4_C2': [4.1, 4.3, 40.2], 'spt4_C3': [2.4, 5.6, 40.2], 'spt4_C4': [2.9, 5.6, 28.4], 'spt4_C5': [2.9, 6.8, 23.2], 'spt4_C3_v2': [2.4, 5.6], 'spt4_C5_v2': [2.9, 6.8]}
            spt4_white_noise_dic['spt4_C3_high220noise'] =  [2.4, 5.6, 40.2]
            if expname.find('herschel')>-1:
                spt4_white_noise_dic = {'spt4_C1plusherschel': [4.1, 4.8, 28.4], 'spt4_C2plusherschel': [4.1, 4.3, 40.2], 'spt4_C3plusherschel': [2.4, 5.6, 40.2], 'spt4_C4plusherschel': [2.9, 5.6, 28.4], 'spt4_C5plusherschel': [2.9, 6.8, 23.2]}
            elif expname.find('planckHF')>-1:
                spt4_white_noise_dic = {'spt4_C1_plusplanckHF': [4.1, 4.8, 28.4], 'spt4_C2_plusplanckHF': [4.1, 4.3, 40.2], 'spt4_C3_plusplanckHF': [2.4, 5.6, 40.2], 'spt4_C4_plusplanckHF': [2.9, 5.6, 28.4], 'spt4_C5_plusplanckHF': [2.9, 6.8, 23.2]}
            elif expname.find('planck')>-1:
                spt4_white_noise_dic = {'spt4_C1_plusplanck': [4.1, 4.8, 28.4], 'spt4_C2_plusplanck': [4.1, 4.3, 40.2], 'spt4_C3_plusplanck': [2.4, 5.6, 40.2], 'spt4_C4_plusplanck': [2.9, 5.6, 28.4], 'spt4_C5_plusplanck': [2.9, 6.8, 23.2]}

            if len(spt4_white_noise_dic[expname]) == 3:
                white_noise_T_225, white_noise_T_286, white_noise_T_345 = spt4_white_noise_dic[expname]
            else:
                white_noise_T_225, white_noise_T_286 = spt4_white_noise_dic[expname]
                white_noise_T_345 = 0.

            if expname == 'spt4_C3_high220noise':
                white_noise_T_220 = 21.
                white_noise_T_225 = 21.


            white_noise_P_225 = np.sqrt(2.) * white_noise_T_225
            white_noise_P_286 = np.sqrt(2.) * white_noise_T_286
            white_noise_P_345 = np.sqrt(2.) * white_noise_T_345

            elknee_T_225, elknee_T_286, elknee_T_345 = 2100., 2100., 2600.
            alphaknee_T_225, alphaknee_T_286, alphaknee_T_345 = 3.9, 3.9, 3.9
            elknee_P_225, elknee_P_286, elknee_P_345 = 200., 200., 200.
            alphaknee_P_225, alphaknee_P_286, alphaknee_P_345 = 2.2, 2.2, 2.2

            if expname.find('plusplanck')>-1:
                elknee_T_345 = -1.
                alphaknee_T_345 = 0.
                elknee_T_345 = -1.
                alphaknee_P_345 = 0.

            '''
            freq_arr = np.asarray( [90., 150., 220., 225., 286., 345.] )
            c = 3e8
            lamb_arr = c/(freq_arr * 1e9)
            D = 8. #10.metres
            beam_arr = np.degrees( 1.22 * lamb_arr / D ) * 60.
            '''
            
            beam_225, beam_286, beam_345 = beam_220, 0.55, 0.45

            specs_dic = {
            90: [beam_90, white_noise_T_90, elknee_T_90, alphaknee_T_90, white_noise_P_90, elknee_P_90, alphaknee_P_90],
            150: [beam_150, white_noise_T_150, elknee_T_150, alphaknee_T_150, white_noise_P_150, elknee_P_150, alphaknee_P_150],
            220: [beam_220, white_noise_T_220, elknee_T_220, alphaknee_T_220, white_noise_P_220, elknee_P_220, alphaknee_P_220],
            225: [beam_225, white_noise_T_225, elknee_T_225, alphaknee_T_225, white_noise_P_225, elknee_P_225, alphaknee_P_225],
            286: [beam_286, white_noise_T_286, elknee_T_286, alphaknee_T_286, white_noise_P_286, elknee_P_286, alphaknee_P_286],
            }
            if white_noise_T_345 != 0.:
                specs_dic[345] = [beam_345, white_noise_T_345, elknee_T_345, alphaknee_T_345, white_noise_P_345, elknee_P_345, alphaknee_P_345]

            corr_noise = 0
            corr_noise_bands = {90:[90], 150:[150], 220: [220], 225:[225], 286:[286], 345:[345]}
            rho = 1.

            cib_corr_coeffs = {}
            cib_corr_coeffs[(90, 150)] = 1.
            cib_corr_coeffs[(90, 220)] = 1.
            cib_corr_coeffs[(90, 225)] = 1.
            cib_corr_coeffs[(90, 286)] = 1.
            cib_corr_coeffs[(90, 345)] = 1.
            cib_corr_coeffs[(150, 220)] = 1.
            cib_corr_coeffs[(150, 225)] = 1.
            cib_corr_coeffs[(150, 286)] = 0.95
            cib_corr_coeffs[(150, 345)] = 0.93
            cib_corr_coeffs[(220, 225)] = 1.
            cib_corr_coeffs[(220, 286)] = 0.98
            cib_corr_coeffs[(220, 345)] = 0.95
            cib_corr_coeffs[(225, 286)] = 0.98
            cib_corr_coeffs[(225, 345)] = 0.95
            cib_corr_coeffs[(286, 345)] = 0.98

        if expname.find('herschel')>-1:

            white_noise_T_600, white_noise_T_857, white_noise_T_1200 = 2.07769677e+03, 3.29858731e+04, 2.51686637e+06 #based on SPTxSPIRE paper

            if expname.find('highnoiseherschel')>-1:
                white_noise_T_600 = white_noise_T_600 * 2.
                white_noise_T_857 = white_noise_T_857 * 2.
                white_noise_T_1200 = white_noise_T_1200 * 2.
        
            white_noise_P_90 = white_noise_T_90 * np.sqrt(2.)
            white_noise_P_150 = white_noise_T_150 * np.sqrt(2.)
            white_noise_P_220 = white_noise_T_220 * np.sqrt(2.)

            white_noise_P_600 = 1e6
            white_noise_P_857 = 1e6
            white_noise_P_1200 = 1e6

            beam_600, beam_857, beam_1200 = 36.6/60., 25.2/60., 18.1/60.
            elknee_T_600 = elknee_T_857 = elknee_T_1200 = -1. ##2500
            alphaknee_T_600 = alphaknee_T_857 = alphaknee_T_1200 = 0. ##4.

            elknee_P_600 = elknee_P_857 = elknee_P_1200 = -1. ##300.
            alphaknee_P_600 = alphaknee_P_857 = alphaknee_P_1200 = 0. ##1.

            specs_dic[600] = [beam_600, white_noise_T_600, elknee_T_600, alphaknee_T_600, white_noise_P_600, elknee_P_600, alphaknee_P_600]
            specs_dic[857] = [beam_857, white_noise_T_857, elknee_T_857, alphaknee_T_857, white_noise_P_857, elknee_P_857, alphaknee_P_857]
            ###specs_dic[1200] = [beam_1200, white_noise_T_1200, elknee_T_1200, alphaknee_T_1200, white_noise_P_1200, elknee_P_1200, alphaknee_P_1200]

            try:
                cib_corr_coeffs
            except:
                cib_corr_coeffs = {}

            if cib_corr_coeffs is None:
                cib_corr_coeffs = {}

            cib_corr_coeffs[(90, 150)] = cib_corr_coeffs[(150, 90)] = 1.
            cib_corr_coeffs[(90, 220)] = cib_corr_coeffs[(220, 90)] = 1.
            cib_corr_coeffs[(150, 220)] = cib_corr_coeffs[(220, 150)] = 1.
            corr_noise_bands[600] = [600]
            corr_noise_bands[857] = [857]
            corr_noise_bands[1200] = [1200]
            cib_corr_coeffs[(90, 600)] = 0.46
            cib_corr_coeffs[(90, 857)] = 0.40
            cib_corr_coeffs[(90, 1200)] = 0.36
            cib_corr_coeffs[(150, 600)] = 0.879
            cib_corr_coeffs[(150, 857)] = 0.789
            cib_corr_coeffs[(150, 1200)] = 0.650
            cib_corr_coeffs[(220, 600)] = 0.879
            cib_corr_coeffs[(220, 857)] = 0.786
            cib_corr_coeffs[(220, 1200)] = 0.643

            cib_corr_coeffs[(225, 600)] = 0.879
            cib_corr_coeffs[(225, 857)] = 0.786
            cib_corr_coeffs[(225, 1200)] = 0.643

            cib_corr_coeffs[(286, 600)] = 0.9
            cib_corr_coeffs[(286, 857)] = 0.8
            cib_corr_coeffs[(286, 1200)] = 0.7

            cib_corr_coeffs[(345, 600)] = 0.93
            cib_corr_coeffs[(345, 857)] = 0.83
            cib_corr_coeffs[(345, 1200)] = 0.73

            cib_corr_coeffs[(600, 857)] = 0.970
            cib_corr_coeffs[(600, 1200)] = 0.861
            cib_corr_coeffs[(857, 1200)] = 0.9551            

        elif expname.find('planckHF')>-1:

            white_noise_T_545, white_noise_T_857 = 1e3, 3e4
        
            white_noise_P_90 = white_noise_T_90 * np.sqrt(2.)
            white_noise_P_150 = white_noise_T_150 * np.sqrt(2.)

            white_noise_P_545 = 1e6
            white_noise_P_857 = 1e6

            beam_545, beam_857 = 4.682, 4.325
            elknee_T_545 = elknee_T_857 = -1
            alphaknee_T_545 = alphaknee_T_857 = 0.

            elknee_P_545 = elknee_P_857 = -1. ##300.
            alphaknee_P_545 = alphaknee_P_857 = 0. ##1.

            specs_dic[545] = [beam_545, white_noise_T_545, elknee_T_545, alphaknee_T_545, white_noise_P_545, elknee_P_545, alphaknee_P_545]
            specs_dic[857] = [beam_857, white_noise_T_857, elknee_T_857, alphaknee_T_857, white_noise_P_857, elknee_P_857, alphaknee_P_857]        

            corr_noise_bands[545] = [545]
            corr_noise_bands[857] = [857]

    return specs_dic, corr_noise_bands, rho, corr_noise, cib_corr_coeffs

def get_spt4_specs(expname, which_spt4 = 'spt4_C3'):
    spt4_white_noise_dic = {'spt4_C1': [4.1, 4.8, 28.4], 'spt4_C2': [4.1, 4.3, 40.2], 'spt4_C3': [2.4, 5.6, 40.2], 'spt4_C4': [2.9, 5.6, 28.4], 'spt4_C5': [2.9, 6.8, 23.2], 'spt4_C3_v2': [2.4, 5.6], 'spt4_C5_v2': [2.9, 6.8]}
    spt4_white_noise_dic['spt4_C3_high220noise'] =  [2.4, 5.6, 40.2]
    if expname.find('herschel')>-1:
        spt4_white_noise_dic = {'spt4_C1plusherschel': [4.1, 4.8, 28.4], 'spt4_C2plusherschel': [4.1, 4.3, 40.2], 'spt4_C3plusherschel': [2.4, 5.6, 40.2], 'spt4_C4plusherschel': [2.9, 5.6, 28.4], 'spt4_C5plusherschel': [2.9, 6.8, 23.2]}
        which_spt4 = '%splusherschel' %(which_spt4)
    elif expname.find('planckHF')>-1:
        spt4_white_noise_dic = {'spt4_C1_plusplanckHF': [4.1, 4.8, 28.4], 'spt4_C2_plusplanckHF': [4.1, 4.3, 40.2], 'spt4_C3_plusplanckHF': [2.4, 5.6, 40.2], 'spt4_C4_plusplanckHF': [2.9, 5.6, 28.4], 'spt4_C5_plusplanckHF': [2.9, 6.8, 23.2]}
        which_spt4 = '%splusplanckHF' %(which_spt4)
    elif expname.find('planck')>-1:
        spt4_white_noise_dic = {'spt4_C1_plusplanck': [4.1, 4.8, 28.4], 'spt4_C2_plusplanck': [4.1, 4.3, 40.2], 'spt4_C3_plusplanck': [2.4, 5.6, 40.2], 'spt4_C4_plusplanck': [2.9, 5.6, 28.4], 'spt4_C5_plusplanck': [2.9, 6.8, 23.2]}
        which_spt4 = '%splusplanck' %(which_spt4)

    if len(spt4_white_noise_dic[which_spt4]) == 3:
        white_noise_T_225, white_noise_T_286, white_noise_T_345 = spt4_white_noise_dic[which_spt4]
    else:
        white_noise_T_225, white_noise_T_286 = spt4_white_noise_dic[which_spt4]
        white_noise_T_345 = 0.

    if expname == 'spt4_C3_high220noise':
        white_noise_T_220 = 21.
        white_noise_T_225 = 21.

    white_noise_P_225 = np.sqrt(2.) * white_noise_T_225
    white_noise_P_286 = np.sqrt(2.) * white_noise_T_286
    white_noise_P_345 = np.sqrt(2.) * white_noise_T_345

    elknee_T_225, elknee_T_286, elknee_T_345 = 2100., 2100., 2600.
    alphaknee_T_225, alphaknee_T_286, alphaknee_T_345 = 3.9, 3.9, 3.9
    elknee_P_225, elknee_P_286, elknee_P_345 = 200., 200., 200.
    alphaknee_P_225, alphaknee_P_286, alphaknee_P_345 = 2.2, 2.2, 2.2

    if expname.find('plusplanck')>-1:
        elknee_T_345 = -1.
        alphaknee_T_345 = 0.
        elknee_T_345 = -1.
        alphaknee_P_345 = 0.

    '''
    freq_arr = np.asarray( [90., 150., 220., 225., 286., 345.] )
    c = 3e8
    lamb_arr = c/(freq_arr * 1e9)
    D = 8. #10.metres
    beam_arr = np.degrees( 1.22 * lamb_arr / D ) * 60.
    '''
    
    beam_225, beam_286, beam_345 = 1.0, 0.55, 0.45

    spt4_specs_dic = {
    225: [beam_225, white_noise_T_225, elknee_T_225, alphaknee_T_225, white_noise_P_225, elknee_P_225, alphaknee_P_225],
    286: [beam_286, white_noise_T_286, elknee_T_286, alphaknee_T_286, white_noise_P_286, elknee_P_286, alphaknee_P_286],
    }
    if white_noise_T_345 != 0.:
        spt4_specs_dic[345] = [beam_345, white_noise_T_345, elknee_T_345, alphaknee_T_345, white_noise_P_345, elknee_P_345, alphaknee_P_345]

    return spt4_specs_dic