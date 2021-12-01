import numpy as np, os, sys
import foregrounds as fg

def get_exp_specs(expname, remove_atm = 0, include_planck = 0):

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

        elif expname == 'spt3g':
            white_noise_T_90 = 3.0
            white_noise_T_150 = 2.2
            white_noise_T_220 = 8.8

            white_noise_P_90 = np.sqrt(2.) * white_noise_T_90
            white_noise_P_150 = np.sqrt(2.) * white_noise_T_150
            white_noise_P_220 = np.sqrt(2.) * white_noise_T_220

        specs_dic = {
        90: [beam_90, white_noise_T_90, elknee_T_90, alphaknee_T_90, white_noise_P_90, elknee_P_90, alphaknee_P_90],
        150: [beam_150, white_noise_T_150, elknee_T_150, alphaknee_T_150, white_noise_P_150, elknee_P_150, alphaknee_P_150],
        220: [beam_220, white_noise_T_220, elknee_T_220, alphaknee_T_220, white_noise_P_220, elknee_P_220, alphaknee_P_220],
        }
        corr_noise = 0
        corr_noise_bands = {90:[90], 150:[150], 220: [220]}
        rho = 1.

        cib_corr_coeffs = None

    if include_planck:

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

    return specs_dic, corr_noise_bands, rho, corr_noise, cib_corr_coeffs
