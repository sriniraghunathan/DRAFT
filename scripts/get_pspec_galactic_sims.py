################################################################################################################

def get_lat_based_masks(nside, max_lat_mask = 50., delta_lat = 10., rotate_to_radec = False):
    """
    #get lat-based masks in ra/dec coordinates
    """
    lat_arr = np.arange( 0., max_lat_mask+1., delta_lat)

    hmask_dic = {}
    npix = H.nside2npix( nside )
    phi_deg, theta_deg = H.pix2ang( nside, np.arange(npix), lonlat = 1 )
    
    for lval in lat_arr:
        curr_mask = np.ones( npix )
        unmask_pixels = np.where( (theta_deg>=-lval) & (theta_deg<=lval) )[0]
        curr_mask[unmask_pixels] = 0.

        if rotate_to_radec:#rotate to ra/dec
            curr_mask = healpix_rotate_coords(curr_mask, coord = ['G', 'C'])

        '''
        H.mollview( curr_mask, coord = ['C', 'G'], sub = (2,2,1) ); H.graticule();
        H.mollview( curr_mask, sub = (2,2,2) ); H.graticule(); show(); #sys.exit()
        '''
        hmask_dic[lval] = curr_mask

    return hmask_dic #these are lat masks in ra/dec coordinates

def split_mask_into_smaller_masks(large_hmask, which_field = 'low_gal', hit_map = None):
    nside = 512 #this is okay to get a rough fsky number
    npix = H.nside2npix(nside) #total number of pixels

    if which_field == 'low_gal': #winter field
        ra1, ra2 = -50., 75.
        dec1, dec2 = -75., -10.

    delra, deldec = 0.1, 0.1
    raarr = np.arange(ra1, ra2, delra)
    decarr = np.arange(dec1, dec2, deldec)

    RA, DEC = np.meshgrid(raarr,decarr)
    RA = RA.flatten()
    DEC = DEC.flatten()

    hmask = np.zeros(npix)
    P=H.pixelfunc.ang2pix(nside,np.radians(90.-DEC),np.radians(RA))
    hmask[P] = 1.

    hmask = H.smoothing(hmask, fwhm = np.radians(2.))
    hmask[hmask<1e-5] = 0.
    hmask[hmask!=0] = 1.    

    nside_out = H.get_nside(large_hmask)
    hmask = H.ud_grade(hmask, nside_out)
    if hit_map is not None:
        hmask *= hit_map

    large_hmask_mod = large_hmask - hmask
    large_hmask_mod[large_hmask_mod<0.] = 0.
    if (0):
        H.mollview(large_hmask, sub = (1,3,1), cmap = cm.jet)
        H.mollview(hmask, sub = (1,3,2), cmap = cm.jet)
        H.mollview(large_hmask_mod, sub = (1,3,3), cmap = cm.jet)    
    return large_hmask_mod, hmask

def healpix_rotate_coords(hmap, coord, beam_to_use_for_smoothing = None, threshold = 0.001):
    """
    coord = ['C', 'G'] to convert a map in RADEC to Gal.    
    """

    #get map pixel
    pixel = np.arange(len(hmap))

    #get angles in this map first
    nside = H.get_nside(hmap)
    angles = H.pix2ang(nside, pixel)

    #roate the angles to the desired new coordinate
    rotated_angles = H.Rotator(coord=coord)(*angles)

    #get the rotated pixel values
    rotated_pixel = H.ang2pix(nside, *rotated_angles)

    #initialise new map
    rot_hmap = np.zeros(len(pixel))

    #push the original map pixel to the new map (in the rotated pixel positions)
    rot_hmap[rotated_pixel] = hmap[pixel]
    if beam_to_use_for_smoothing is not None:
        rot_hmap = H.smoothing(rot_hmap, fwhm = np.radians(beam_to_use_for_smoothing))
        rot_hmap[rot_hmap>threshold] = 1.
        rot_hmap[rot_hmap<threshold] = 0.

    return rot_hmap

def get_spt3g_mask(which_field, nside_out = 4096):

    nside = 512 #this is okay to get a rough fsky number
    npix = H.nside2npix(nside) #total number of pixels

    if which_field == 'winter_field': #winter field
        ra1, ra2 = -50., 50.
        dec1, dec2 = -70., -42.
        #dec1, dec2 = -64., -42.

    elif which_field == 'summer_field_el1c_el2c': #summer field
        #20220225
        #ra1, ra2 = 155., 205.
        #dec1, dec2 = -42., -28.
        #https://pole.uchicago.edu/spt3g/index.php/Observing_Cadence_Summer_2020_2021#Field_Definition
        ra1, ra2 = 155., 225.
        dec1, dec2 = -42., -28.

    elif which_field == 'summer_field_el1b_el2b': #summer field
        ra1, ra2 = 0., 50.
        dec1, dec2 = -42., -28.

    elif which_field == 'summer_field_el1_el5': #summer field
        ra1, ra2 = 50., 100.
        dec1, dec2 = -63., -28.

    #20220225 - https://pole.uchicago.edu/spt3g/index.php/Observing_Cadence_Summer_2020_2021#Field_Definition
    elif which_field == 'summer_field_el1_el2': #summer field
        ra1, ra2 = 50., 100.
        dec1, dec2 = -42., -28.

    delra, deldec = 0.1, 0.1
    raarr = np.arange(ra1, ra2, delra)
    decarr = np.arange(dec1, dec2, deldec)

    RA, DEC = np.meshgrid(raarr,decarr)
    RA = RA.flatten()
    DEC = DEC.flatten()

    HMASK = np.zeros(npix)
    P=H.pixelfunc.ang2pix(nside,np.radians(90.-DEC),np.radians(RA))
    HMASK[P] = 1.
    HMASK = H.smoothing(HMASK, fwhm = np.radians(2.))
    HMASK[HMASK<1e-5] = 0.

    HMASK = H.ud_grade(HMASK, nside_out)

    return HMASK

def set_telescope_observer(lat = -89.991066, lon = -44.65, elevation = 2835.0, astropy = 1):#, location = 'southpole'):
    '''
    if location == 'southpole':
        lat, lon, elevation = -22.57, -67.47, 5200. #Chile
    elif location == 'chile':
        lat, lon, elevation = -22.57, -67.47, 5200. #Chile
    '''
    splat = coord.EarthLocation(lat=lat*u.deg,lon=lon*u.deg, height=elevation*u.meter)

    return splat

def get_footprint_for_min_obs_el(nside, min_obs_el, max_obs_el = 85., timearr = None, ra_dec_arr = None, yyyy = 2029, delta_mm = 1, delta_dd = 30., delta_hh = 24., location = 'southpole', beam_to_use_for_smoothing = 0.5, threshold = 0.001):#, epoch = 'J2000'):
    if ra_dec_arr is None:
        raarr = np.arange(-180, 180., 0.1)
        decarr = np.arange(-90., 90., 0.1)

        ra_dec_arr = []
        for r in raarr:
            for d in decarr:
                ra_dec_arr.append( [r, d])
        ra_dec_arr = np.asarray( ra_dec_arr )

    if timearr is None: #does not matter for pole.
        timearr = []
        for mm in np.arange(1, 13, delta_mm):
            for dd in np.arange(1, 30, delta_dd):
                for hh in np.arange(1, 24, delta_hh):
                    tval = '%04d-%02d-%02dT%02d:00:00' %(yyyy, mm, dd, hh)
                    try:
                        t = Time(tval, format='isot', scale='utc')
                        timearr.append( tval )
                    except:
                        pass
        timearr = Time(timearr, format='isot', scale='utc')

    astropy_coord_radec = SkyCoord(ra_dec_arr[:,0], ra_dec_arr[:,1], unit='deg', frame='icrs')
    splat = set_telescope_observer()
    npix = H.nside2npix(nside)
    hmask = np.zeros(npix)

    for t in timearr:
        astropy_coord_altaz = astropy_coord_radec.transform_to(coord.AltAz(obstime=t,location=splat))
        az, el = astropy_coord_altaz.az.value, astropy_coord_altaz.alt.value
        good_inds = np.where( (el>=min_obs_el) & (el<=max_obs_el) )[0]
        good_ra_dec = ra_dec_arr[good_inds]

        good_pixels = H.ang2pix(nside, np.radians(90.-good_ra_dec[:,1]), np.radians(good_ra_dec[:,0]))

        hmask[good_pixels] = 1.
        if location == 'southpole': break

    hmask = H.smoothing(hmask, fwhm = np.radians(beam_to_use_for_smoothing))    
    hmask[hmask<threshold] = 0.
    hmask[hmask>threshold] = 1.

    return hmask
############################################################
############################################################

import healpy as H, numpy as np, glob, sys, os, argparse
from astropy.time import Time
from astropy import units as u
from astropy import coordinates as coord
from astropy.coordinates import SkyCoord, EarthLocation

local = 1
if str(os.getcwd()).find('sri')>-1: local = 0


#dust_or_sync = sys.argv[1] ##'sync' ##'dust'

parser = argparse.ArgumentParser(description='')
parser.add_argument('-dust_or_sync', dest='dust_or_sync', action='store', help='dust_or_sync', type=str, required=True)
parser.add_argument('-lmax', dest='lmax', action='store', help='lmax', type=int, default=7000)
parser.add_argument('-nuarr', dest='nuarr', action='store', type=int, nargs='+', default= [27, 39, 93, 145, 225, 278], help='nuarr')
#parser.add_argument('-nuarr', dest='nuarr', action='store', type=int, nargs='+', default= [93, 145, 225, 278], help='nuarr')
parser.add_argument('-t_only', dest='t_only', action='store', help='t_only', type=int, default=0)
parser.add_argument('-nside', dest='nside', action='store', help='nside', type=int, default=2048)#4096)
parser.add_argument('-verbose', dest='verbose', action='store', help='verbose', type=int, default=0)
parser.add_argument('-zonca_sims', dest='zonca_sims', action='store', help='zonca_sims', type=int, default=1) #S4 sims by Andrea Zonca
parser.add_argument('-pySM_yomori', dest='pySM_yomori', action='store', help='pySM_yomori', type=int, default=0) #pySM sims by Yuuki Omori

###parser.add_argument('-use_planck_mask', dest='use_planck_mask', action='store', help='use_planck_mask', type=int, default=0) #use Planck galactic mask
parser.add_argument('-which_mask', dest='which_mask', action='store', help='which_mask', type=int, default=-1)

parser.add_argument('-use_lat_step_mask', dest='use_lat_step_mask', action='store', help='use_lat_step_mask', type=int, default=0) #mask based on latitude

parser.add_argument('-use_s4like_mask', dest='use_s4like_mask', action='store', help='use_s4like_mask', type=int, default=0) #rough S4 mask
parser.add_argument('-use_s4like_mask_v2', dest='use_s4like_mask_v2', action='store', help='use_s4like_mask_v2', type=int, default=0) #rough S4 mask v2: split footrpint into clean and unclean region
parser.add_argument('-use_s4like_mask_v3', dest='use_s4like_mask_v3', action='store', help='use_s4like_mask_v3', type=int, default=0) #rough S4 mask v3: split v2 clean into clean and cleaner regions
parser.add_argument('-use_s4delensing_mask', dest='use_s4delensing_mask', action='store', help='use_s4delensing_mask', type=int, default=0) #no mask * delensing footprint (or hit map actually).
parser.add_argument('-cos_el', dest='cos_el', action='store', help='cos_el', type=int, default=40) #rough S4 mask v2: split footrpint into clean and unclean region
parser.add_argument('-use_spt3g_mask', dest='use_spt3g_mask', action='store', help='use_spt3g_mask', type=int, default=0) #SPT-3G summer/winter field
parser.add_argument('-use_planck_mask', dest='use_planck_mask', action='store', help='use_planck_mask', type=int, default=0) #use Planck galactic mask

#20210111 - CMB-S4 SP-LAT - different min obs el and gal cuts
parser.add_argument('-min_obs_el', dest='min_obs_el', action='store', help='min_obs_el', type=float, default=30.)
parser.add_argument('-use_splat_minobsel_galcuts', dest='use_splat_minobsel_galcuts', action='store', help='use_splat_minobsel_galcuts', type=int, default=0)

#20230330 - SPT proposal 2023
parser.add_argument('-which_spt_field', dest='which_spt_field', action='store', help='which_spt_field', type=str, default='nominal')


parser.add_argument('-testing', dest='testing', action='store', help='testing', type=int, default=0)

#nuarr = [20, 27, 39, 93, 145, 225, 278]
#nuarr = [27, 39, 93, 145, 225, 278]


args = parser.parse_args()
args_keys = args.__dict__
for kargs in args_keys:
    param_value = args_keys[kargs]
    if isinstance(param_value, str):
        cmd = '%s = "%s"' %(kargs, param_value)
    else:
        cmd = '%s = %s' %(kargs, param_value)
    exec(cmd)

name_dic = {}
if zonca_sims:
    name_dic[20] = 'ULFL1'
    name_dic[27] = 'LFL1'
    name_dic[39] = 'LFL2'
    name_dic[93] = 'MFL1'
    name_dic[145] = 'MFL2'
    name_dic[225] = 'HFL1'
    name_dic[278] = 'HFL2'

if use_spt3g_mask:
    ##nuarr = [93, 145]##, 225]
    nuarr = [93, 145, 225]

if pySM_yomori:
    #nuarr = [90, 100, 143, 150, 217, 220, 353]
    nuarr = [143, 150]

if testing or local:
    lmax = 2000
    nside = 512
    ##nuarr = [ 145 ]#, 145]
    nuarr = [ 145 ]
    if pySM_yomori:
        nuarr = [ 143 ]
    t_only = 0 #1

'''
if not local:
    mask_folder = '/u/home/s/srinirag/project-nwhiteho/cmbs4/masks/planck/'
    cmbs4_footprint_folder = '/u/home/s/srinirag/project-nwhiteho/cmbs4/footprints/'
else:
    mask_folder = '/Volumes/data_PHD_WD_babbloo/s4/cmbs4/masks/planck/'
    cmbs4_footprint_folder = '/Volumes/data_PHD_WD_babbloo/s4/cmbs4/footprints/'
'''

mask_folder = 'cmbs4/masks/planck/'
cmbs4_footprint_folder = 'cmbs4/footprints/'

if local:
    mask_folder = '/Volumes/data_PHD_WD_babbloo/s4/cmbs4/masks/planck/'
    cmbs4_footprint_folder = '/Volumes/data_PHD_WD_babbloo/s4/cmbs4/footprints/'

if use_s4delensing_mask:
    cmbs4_footprint_folder = '%s/LAT_pole/LAT-MFPL1_pole/' %(cmbs4_footprint_folder)

if zonca_sims or pySM_yomori:
    '''
    if local:
        sim_folder = '/Volumes/data_PHD_WD_babbloo/s4/cmbs4/map_based_simulations/202002_foregrounds_extragalactic_cmb_tophat/4096/'
    else:
        sim_folder = '/u/home/s/srinirag/project-nwhiteho/cmbs4/map_based_simulations/202002_foregrounds_extragalactic_cmb_tophat/4096/'
    '''
    if zonca_sims:
        sim_folder = 'cmbs4/map_based_simulations/202002_foregrounds_extragalactic_cmb_tophat/4096/'
        if local:
            sim_folder = '/Volumes/data_PHD_WD_babbloo/s4/cmbs4/map_based_simulations/202002_foregrounds_extragalactic_cmb_tophat/4096/'

        if (1): #20210323
            sim_folder = sim_folder.replace('202002_foregrounds_extragalactic_cmb_tophat', '202102_design_tool_input')

        if dust_or_sync == 'dust':
            sim_folder = '%s/dust/' %(sim_folder)
        elif dust_or_sync == 'sync':
            sim_folder = '%s/synchrotron/' %(sim_folder)
        elif dust_or_sync == 'freefree':
            sim_folder = '%s/freefree/' %(sim_folder)

        sim_folder = '%s/0000/' %(sim_folder)
    else:
        sim_folder = 'pySM/yomori/galforegrounds/'

    if use_planck_mask:
        ##opfname = '%s/cls_galactic_sims_%s_maskplanck_nside%s_lmax%s.npy' %(sim_folder, dust_or_sync, nside, lmax)
        opfname = '%s/planck_mask/cls_galactic_sims_%s_nside%s_lmax%s.npy' %(sim_folder, dust_or_sync, nside, lmax)
    elif use_lat_step_mask:
        opfname = '%s/lat_steps/cls_galactic_sims_%s_nside%s_lmax%s.npy' %(sim_folder, dust_or_sync, nside, lmax)
    elif use_s4like_mask:
        opfname = '%s/s4like_mask/cls_galactic_sims_%s_nside%s_lmax%s.npy' %(sim_folder, dust_or_sync, nside, lmax)
    elif use_s4like_mask_v2:
        opfname = '%s/s4like_mask_v2/cls_galactic_sims_%s_nside%s_lmax%s.npy' %(sim_folder, dust_or_sync, nside, lmax)
    elif use_s4like_mask_v3:
        opfname = '%s/s4like_mask_v3/cls_galactic_sims_%s_nside%s_lmax%s.npy' %(sim_folder, dust_or_sync, nside, lmax)
    elif use_s4delensing_mask:
        opfname = '%s/s4delensing_mask/cls_galactic_sims_%s_nside%s_lmax%s.npy' %(sim_folder, dust_or_sync, nside, lmax)
    elif use_spt3g_mask:
        if which_spt_field == 'nominal':
            opfname = '%s/spt3g/cls_galactic_sims_%s_nside%s_lmax%s.npy' %(sim_folder, dust_or_sync, nside, lmax)
        else:
            opfname = '%s/spt3g/%s/cls_galactic_sims_%s_nside%s_lmax%s.npy' %(sim_folder, which_spt_field, dust_or_sync, nside, lmax)
    elif use_splat_minobsel_galcuts:
        opfname = '%s/s4_use_splat_minobsel%s_galcuts/cls_galactic_sims_%s_nside%s_lmax%s.npy' %(sim_folder, min_obs_el, dust_or_sync, nside, lmax)
else:
    if local:
        sim_folder = '/Users/sraghunathan/Research/SPTPol/analysis/git/ilc/galactic/CUmilta/ampmod_maps/'
    else:
        sim_folder = 'S4_march_2020/sims_from_others/CUmilta/ampmod_maps/'

    if use_planck_mask:
        #opfname = '%s/cls_gal_%s_nside%s_lmax%s.npy' %(sim_folder, dust_or_sync, nside, lmax)
        ##opfname = '%s/cls_galactic_sims_%s_CUmilta_20200319_maskplanck_nside%s_lmax%s.npy' %(sim_folder, dust_or_sync, nside, lmax)
        opfname = '%s/planck_mask/cls_galactic_sims_%s_CUmilta_20200319_nside%s_lmax%s.npy' %(sim_folder, dust_or_sync, nside, lmax)
    elif use_lat_step_mask:
        opfname = '%s/lat_steps/cls_galactic_sims_%s_CUmilta_20200319_nside%s_lmax%s.npy' %(sim_folder, dust_or_sync, nside, lmax)
    elif use_s4like_mask:
        opfname = '%s/s4like_mask/cls_galactic_sims_%s_CUmilta_20200319_nside%s_lmax%s.npy' %(sim_folder, dust_or_sync, nside, lmax)
    elif use_s4like_mask_v2:
        opfname = '%s/s4like_mask_v2/cls_galactic_sims_%s_CUmilta_20200319_nside%s_lmax%s.npy' %(sim_folder, dust_or_sync, nside, lmax)
    elif use_s4like_mask_v3:
        opfname = '%s/s4like_mask_v3/cls_galactic_sims_%s_CUmilta_20200319_nside%s_lmax%s.npy' %(sim_folder, dust_or_sync, nside, lmax)
    elif use_s4delensing_mask:
        opfname = '%s/s4delensing_mask/cls_galactic_sims_%s_CUmilta_20200319_nside%s_lmax%s.npy' %(sim_folder, dust_or_sync, nside, lmax)

if use_planck_mask: os.system('mkdir %s/planck_mask/' %(sim_folder))
if use_lat_step_mask: os.system('mkdir %s/lat_steps/' %(sim_folder))
if use_s4like_mask: os.system('mkdir %s/s4like_mask/' %(sim_folder))
if use_s4like_mask_v2: os.system('mkdir %s/s4like_mask_v2/' %(sim_folder))
if use_s4like_mask_v3: os.system('mkdir %s/s4like_mask_v3/' %(sim_folder))
if use_s4delensing_mask: os.system('mkdir %s/s4delensing_mask/' %(sim_folder))
if use_spt3g_mask: os.system('mkdir %s/spt3g/' %(sim_folder))
if use_splat_minobsel_galcuts: os.system('mkdir %s/s4_use_splat_minobsel%s_galcuts/' %(sim_folder, min_obs_el))

if not os.path.exists('tmp/'): os.system('mkdir tmp/')
log_file = 'tmp/pspec_%s.txt' %(dust_or_sync)

if t_only:
    opfname = opfname.replace('.npy', '_TTonly.npy')

if which_mask != -1:
    log_file = log_file.replace('.txt', '_mask%s.txt' %(which_mask))     
    opfname = opfname.replace('.npy', '_mask%s.npy' %(which_mask))

if not use_spt3g_mask:##cos_el != 30:
    if use_s4delensing_mask:
        log_file = log_file.replace('.txt', '_delensing.txt') 
        opfname = opfname.replace('.npy', '_delensing.npy')
    elif use_splat_minobsel_galcuts:
        log_file = log_file.replace('.txt', '_minobsel%s.txt' %(min_obs_el))
        opfname = opfname.replace('.npy', '_minobsel%s.npy' %(min_obs_el))
    else:
        log_file = log_file.replace('.txt', '_cos_el_%s.txt' %(cos_el))     
        opfname = opfname.replace('.npy', '_cos_el_%s.npy' %(cos_el))

lf = open(log_file, 'w'); lf.close()
if (1):#testing or not local:

    totthreads = 2
    os.putenv('OMP_NUM_THREADS',str(totthreads))

    #get filename prefix
    if zonca_sims:
        fname_pref = 'cmbs4_%s_uKCMB_LAT-xxx_nside4096_0000.fits' %(dust_or_sync)
        fname_pref = fname_pref.replace('sync', 'synchrotron')
    elif pySM_yomori:
        #fname_pref = 'map4096_yyy_xxxGHz_d1s1.fits'
        fname_pref = 'map4096_yyy_xxxGHz_%s0.fits' %(dust_or_sync[0]) #just pick the first letter: d/s for dust/sync
    else:
        fname_pref = 'Ampmod_map_%s_lmax_5200_freq_xxx_rel0000.fits' %(dust_or_sync)

    logline = '\n'
    lf = open(log_file,'a'); lf.writelines('%s\n' %(logline));lf.close()
    print(logline)

    map_dic = {}
    ###print(nuarr); sys.exit()
    for nucntr, nu in enumerate( nuarr ):
        fname = '%s/%s' %(sim_folder, fname_pref)
        if zonca_sims:
            fname = fname.replace('xxx', name_dic[nu])
        elif pySM_yomori:
            if nu in [90, 150, 220]:
                plancksptstr = 'spt3g'
            elif  nu in [100, 143, 217, 353]:
                plancksptstr = 'planck'
            fname = fname.replace('yyy', plancksptstr).replace('xxx', '%03d' %(nu))

        else:
            fname = fname.replace('xxx', '%03d' %(nu))

        logline = '\t%s\n' %fname
        lf = open(log_file,'a'); lf.writelines('%s\n' %(logline));lf.close()
        print(logline)

        if t_only:
            currmap = H.read_map(fname, verbose = 0)
        else:
            currmap = H.read_map(fname, verbose = 0, field = (0,1,2))
        if H.get_nside(currmap) != nside:
            currmap = H.ud_grade(currmap, nside_out = nside)

        ##from IPython import embed; embed()

        map_dic[nu] = currmap

    if not use_spt3g_mask: #get cmbs4 footprint if not SPT
        if use_splat_minobsel_galcuts:
            cmbs4_hit_map = get_footprint_for_min_obs_el(nside, min_obs_el = min_obs_el)
        else:
            if not use_s4delensing_mask:
                cmbs4_hit_map_fname = '%s/high_cadence_hits_el%s_cosecant_modulation.fits' %(cmbs4_footprint_folder, cos_el)
            else:
                cmbs4_hit_map_fname = '%s/cmbs4_hitmap_LAT-MFPL1_pole_nside4096_1_of_1.fits' %(cmbs4_footprint_folder)
            #cmbs4_hit_map_fname = '%s/high_cadence_hits_el40_cosecant_modulation.fits' %(cmbs4_footprint_folder)
            cmbs4_hit_map = H.read_map(cmbs4_hit_map_fname, verbose = verbose)
            cmbs4_hit_map = cmbs4_hit_map/np.max(cmbs4_hit_map)
            ###20220325 - removing this binary mask
            cmbs4_hit_map[cmbs4_hit_map!=0] = 1.
            if H.get_nside(cmbs4_hit_map) != nside:
                cmbs4_hit_map = H.ud_grade(cmbs4_hit_map, nside_out = nside)

    #now get masks
    if use_planck_mask:
        use_apod_planck = False
        if not use_apod_planck:
            planck_mask_fname = '%s/HFI_Mask_GalPlane-apo0_2048_R2.00.fits' %(mask_folder)
        else:
            planck_mask_fname = '%s/HFI_Mask_GalPlane-apo2_2048_R2.00.fits' %(mask_folder)
        ###planck_mask = H.read_map(planck_mask_fname, verbose = verbose, field = (1,2,3))
        """
        #20230302 - this is true for unapodised masks. Apodised ones have slightly lower fsky.
        TTYPE1  = 'GAL020  '           / 20% sky coverage                               
        TTYPE2  = 'GAL040  '           / 40% sky coverage                               
        TTYPE3  = 'GAL060  '           / 60% sky coverage                               
        TTYPE4  = 'GAL070  '           / 70% sky coverage                               
        TTYPE5  = 'GAL080  '           / 80% sky coverage                               
        TTYPE6  = 'GAL090  '           / 90% sky coverage                               
        TTYPE7  = 'GAL097  '           / 97% sky coverage                               
        TTYPE8  = 'GAL099  '           / 99% sky coverage                               
        """
        planck_mask = H.read_map(planck_mask_fname, verbose = verbose, field = (3, 4, 5)) #GAL070, GAL080, GAL090.
        if (0):

            import healpy as H
            use_apod = False #True
            if use_apod:
                planck_mask_fname = '/Users/sraghunathan/Downloads/HFI_Mask_GalPlane-apo2_2048_R2.00.fits'
                plname = '/Users/sraghunathan/Downloads/planck_masks_20230302.pdf'
            else:
                planck_mask_fname = '/Users/sraghunathan/Downloads/HFI_Mask_GalPlane-apo0_2048_R2.00.fits'
                plname = '/Users/sraghunathan/Downloads/planck_masks_noapod_20230302.pdf'

            planck_mask = H.read_map(planck_mask_fname, field = (3, 4, 5)) #

            clf()
            figure(figsize=(10., 11))
            tit_arr = ['GAL070', 'GAL080', 'GAL090']
            for cntr, (m, t) in enumerate(zip(planck_mask, tit_arr)):
                H.mollview(m, sub = (3, 1, cntr+1), title = r'%s: $f_{\rm sky}$ = %.2f' %(t, np.mean(m)))
            savefig(plname, dpi = 200.)  

        if H.get_nside(planck_mask) != nside:
            planck_mask = H.ud_grade(planck_mask, nside_out = nside)

        '''
        planck_mask = H.smoothing(np.copy(planck_mask), fwhm = np.radians(5.), lmax = lmax, verbose = verbose)
        thresh = 0.4
        for mask_iter in range(len(planck_mask)):
            planck_mask[mask_iter][planck_mask[mask_iter]<thresh] = 0.
            planck_mask[mask_iter][planck_mask[mask_iter]!=0] = 1.
        '''

        tot_masks = len(planck_mask)

    elif use_splat_minobsel_galcuts: #20210111 - CMB-S4 SP-LAT - different min obs el and gal cuts
        lat_mask_dic = get_lat_based_masks(nside)#, use_lat_steps = use_lat_steps)
        splat_minobsel_galcuts_mask_arr = []
        for galcut_bval in lat_mask_dic:
            splat_minobsel_galcuts_mask_arr.append( lat_mask_dic[galcut_bval] )
        tot_masks = len(splat_minobsel_galcuts_mask_arr)

    elif use_lat_step_mask:

        '''
        H.mollview(cmbs4_hit_map, coord = ['C'], sub = (2,2,1)); H.graticule()
        H.mollview(planck_mask, coord = ['G', 'C'], sub = (2,2,2)); H.graticule()

        H.mollview(cmbs4_hit_map, coord = ['C', 'G'], sub = (2,2,3)); H.graticule()
        H.mollview(planck_mask, coord = ['G'], sub = (2,2,4)); H.graticule(); show()
        '''

        if use_lat_step_mask:
            min_lat, max_lat, delta_lat = -60., 30., 15.
        lat_arr = np.arange( min_lat, max_lat + 1., delta_lat )

        npix = H.nside2npix( nside )
        phi_deg, theta_deg = H.pix2ang( nside, np.arange(npix), lonlat = 1 )
        
        lat_mask_arr = []
        for l1 in lat_arr[:-1]:
            l2 = l1 + delta_lat
            curr_mask = np.zeros( npix )
            unmask_pixels = np.where( (theta_deg>=l1) & (theta_deg<l2) )[0]
            curr_mask[unmask_pixels] = 1.

            '''
            print(np.min( theta_deg[unmask_pixels] ), np.max( theta_deg[unmask_pixels] ))
            H.mollview( curr_mask, sub = (2,2,1)); H.graticule(); 
            H.mollview( cmbs4_hit_map, sub = (2,2,2)); H.graticule(); show()
            '''

            lat_mask_arr.append(curr_mask)

        #finally add all masks together
        curr_mask = np.sum(lat_mask_arr, axis = 0)
        lat_mask_arr.append( curr_mask )

        tot_masks = len(lat_mask_arr)

    elif use_s4like_mask or use_s4like_mask_v2 or use_s4like_mask_v3:

        npix = H.nside2npix( nside )
        phi_deg, theta_deg = H.pix2ang( nside, np.arange(npix), lonlat = 1 )

        s4_mask_dic = {2: 10., 1: 15., 0: 20.}

        if use_s4like_mask_v3: #20210204: only 1 clean / dirty mas
            s4_mask_dic = {0: 10.}
        
        s4like_mask_arr = []
        for mask_iter in sorted( s4_mask_dic ):
            min_lat, max_lat = -s4_mask_dic[mask_iter], s4_mask_dic[mask_iter]
            curr_mask = np.ones( npix )
            unmask_pixels = np.where( (theta_deg>=min_lat) & (theta_deg<max_lat) )[0]
            curr_mask[unmask_pixels] = 0.

            s4like_mask_arr.append( curr_mask )

        tot_masks = len(s4like_mask_arr)

    elif use_s4delensing_mask:

        npix = H.nside2npix( nside )
        s4delensing_mask_arr = [ np.ones(npix) ]
        tot_masks = len(s4delensing_mask_arr)       

    elif use_spt3g_mask: #get SPT-3G summer/winter field mask

        if which_spt_field == 'nominal':
            spt_mask_winter = get_spt3g_mask('winter_field', nside_out = nside)
            spt_mask_summer_el1c_el2c = get_spt3g_mask('summer_field_el1c_el2c', nside_out = nside)
            spt_mask_summer_el1b_el2b = get_spt3g_mask('summer_field_el1b_el2b', nside_out = nside)
            spt_mask_summer_el1_el5 = get_spt3g_mask('summer_field_el1_el5', nside_out = nside)

            spt_mask_arr = [spt_mask_winter, spt_mask_summer_el1c_el2c, spt_mask_summer_el1b_el2b, spt_mask_summer_el1_el5]

            tot_masks = len(spt_mask_arr)

        elif which_spt_field == 'spt_proposal_2023_13k_sqdeg_field_mask':

            mask_fd = 'masks/spt_proposal_2023_13k_sqdeg_field/'
            fname_arr = ['%s/spt_proposal_2023_13k_sqdeg_field_planckGAL070.fits' %(mask_fd), 
                        '%s/spt_proposal_2023_13k_sqdeg_field_planckGAL080.fits' %(mask_fd),
                        '%s/spt_proposal_2023_13k_sqdeg_field_planckGAL090.fits' %(mask_fd)]

            spt_mask_arr = []
            for tmpfname in fname_arr:
                spt_mask_arr.append( H.read_map( tmpfname ) )

        tot_masks = len(spt_mask_arr)

    logline = '\tget masks now\n'
    lf = open(log_file,'a'); lf.writelines('%s\n' %(logline));lf.close()
    print(logline)

    mask_arr = []

    for mask_iter in range(tot_masks):

        if use_planck_mask:

            mask = np.copy(planck_mask[mask_iter])

        elif use_splat_minobsel_galcuts:

            mask = splat_minobsel_galcuts_mask_arr[mask_iter]

        elif use_lat_step_mask:

            mask = lat_mask_arr[mask_iter]

        elif use_s4like_mask or use_s4like_mask_v2 or use_s4like_mask_v3:

            mask = s4like_mask_arr[mask_iter]

        elif use_s4delensing_mask:

            mask = s4delensing_mask_arr[mask_iter]

        elif use_spt3g_mask:

            mask = spt_mask_arr[mask_iter]

        rotated = 0
        if not use_spt3g_mask:
            if use_planck_mask:
                mask = H.ud_grade(mask, 128)
            else:
                mask = H.ud_grade(mask, 256)

            mask = healpix_rotate_coords(mask, coord = ['G', 'C'])
            if use_planck_mask and not use_apod_planck:
                mask = H.smoothing(mask, fwhm = np.radians(2.))
                mask[mask<1e-5] = 0.

            if use_planck_mask:
                import matplotlib
                matplotlib.use('Agg')
                from pylab import *
                clf()
                H.mollview(mask, title = r'%s: $f_{\rm sky}$ = %.2f' %(mask_iter, np.mean(mask)), cbar = False)
                pl_folder = 'plots/planck_masks'
                if not os.path.exists( pl_folder ): os.system('mkdir -p %s' %(pl_folder))
                plname = '%s/planck_mask_%s.png' %(pl_folder, mask_iter)
                savefig(plname); close()

            rotated = 1

        if pySM_yomori: #all masks must be in galactic coordinates for original pySM sims
            mask = H.ud_grade(mask, 128)
            mask = healpix_rotate_coords(mask, coord = ['C', 'G'])
            rotated = 1

        if rotated:
            #simple rotation from gal to celestial
            mask = H.smoothing(np.copy(mask), fwhm = np.radians(10.), verbose = verbose)#, lmax = lmax)
            mask = H.ud_grade(mask, nside_out = nside)
            thresh = 0.4
            mask[mask<thresh] = 0.
            mask[mask!=0] = 1.

        mask_arr.append( mask )

    if use_planck_mask or use_s4like_mask:

        #first full sky
        npix = H.nside2npix( nside )    
        no_mask = np.ones( npix )
        #mask_arr = np.concatenate(([no_mask], mask_arr))
        mask_arr.append( no_mask )

        mask_arr = mask_arr * cmbs4_hit_map
        tot_masks = len(mask_arr)

    #if use_s4like_mask_v2 or use_s4like_mask_v3: #20200627
    if use_s4like_mask_v2 or use_s4like_mask_v3 or use_s4delensing_mask: #20210415

        mask_arr = mask_arr * cmbs4_hit_map
        #add additional masks that cover the uncleaned CMB-S4 region.
        unclean_mask_arr = []
        for m in range(tot_masks):
            curr_clean_mask = mask_arr[m]
            curr_unclean_mask = cmbs4_hit_map - curr_clean_mask
            unclean_mask_arr.append( curr_unclean_mask )
            '''
            tit1 = np.mean(curr_clean_mask);tit2 = np.mean(curr_unclean_mask);print(tit1 + tit2)
            H.mollview(curr_clean_mask, sub = (1,2,1), title = '%.2f' %tit1);H.mollview(curr_unclean_mask, sub = (1,2,2), title = '%.2f' %tit2); show()
            '''
        mask_arr = np.concatenate( (mask_arr, unclean_mask_arr) )
        tot_masks = len(mask_arr)

    if use_splat_minobsel_galcuts: #20220111 - CMB-S4 SP-LAT
        mask_arr = mask_arr * cmbs4_hit_map
        if (0):
            for i in range(len(mask_arr)):
                H.mollview(mask_arr[i]); show()

    if use_s4like_mask_v3: #20210204: split clean regiosn into 2 sub-fields
        tots4masks = len(s4like_mask_arr)
        clean_masks = mask_arr[:tots4masks]
        dirty_masks = mask_arr[tots4masks:]
        mod_clean_masks_arr = []
        print(tots4masks)
        for m in range(tots4masks):
            curr_clean_mask = clean_masks[m]
            curr_clean_mask_large, curr_clean_mask_small = split_mask_into_smaller_masks(curr_clean_mask, hit_map = cmbs4_hit_map)        
            mod_clean_masks_arr.append(curr_clean_mask_small)
            mod_clean_masks_arr.append(curr_clean_mask_large)
            mod_clean_masks_arr.append(dirty_masks[m])
        
        mask_arr = np.copy(mod_clean_masks_arr)


    #print(len(mask_arr))
    if which_mask != -1:
        mask_arr = [mask_arr[which_mask]]

    mask_arr = np.asarray(mask_arr)
    tot_masks = len(mask_arr)

    fsky_arr = np.mean(mask_arr, axis = 1)

    logline = '\t\t all masks obtained\n'
    lf = open(log_file,'a'); lf.writelines('%s\n' %(logline));lf.close()
    print(logline)

    if testing:
        from IPython import embed; embed()
        from pylab import *

        from matplotlib import rc;rc('text', usetex=True);rc('font', weight='bold');matplotlib.rcParams['text.latex.preamble'] = r'\boldmath'
        rcParams['figure.dpi'] = 150
        rcParams["figure.facecolor"] = 'white'
        rcParams['font.family'] = 'serif'

        rc('text.latex',preamble=r'\usepackage{/Users/sraghunathan/.configs/apjfonts}')

        if use_lat_step_mask:
            tot_masks = tot_masks - 1
        clf()
        for mask_iter in range(tot_masks):
            fsky = np.mean(mask_arr[mask_iter])
            if not use_s4like_mask_v2:
                H.mollview(mask_arr[mask_iter], sub = (1, tot_masks,mask_iter+1), title = r'Mask: %s: f$_{\rm sky} = %.2f$' %(mask_iter, fsky), cbar = 0, title_fontsize = 10); 
            else:
                if mask_iter<3:
                    maskno = '%s(a)' %(mask_iter)
                else:
                    maskno = '%s(b)' %(mask_iter-3)
                H.mollview(mask_arr[mask_iter], sub = (2, tot_masks/2,mask_iter+1), title = r'Mask: %s: f$_{\rm sky} = %.2f$' %(maskno, fsky), cbar = 0, title_fontsize = 10); 

        show(); sys.exit()
        plname = '/Users/sraghunathan/Research/SPTPol/analysis/git/DRAFT/reports/galactic_sims/maps_masks/masks.pdf'
        if use_lat_step_mask:
            plname = '/Users/sraghunathan/Research/SPTPol/analysis/git/DRAFT/reports/galactic_sims/maps_masks/masks_S4wide_cluster_search.pdf'
        if use_s4like_mask:
            plname = '/Users/sraghunathan/Research/SPTPol/analysis/git/DRAFT/reports/galactic_sims/maps_masks/masks_S4_Neff.pdf'
        if use_s4like_mask_v2:
            plname = '/Users/sraghunathan/Research/SPTPol/analysis/git/DRAFT/reports/galactic_sims/maps_masks/masks_S4_v2_Neff_cos_el_%s.png' %(cos_el)
        plfolder = '/'.join( plname.split('/')[:-1] )
        os.system('mkdir -p %s' %(plfolder))
        savefig(plname, dpi = 150)
        show()

        #mask_arr = [mask_arr[2], mask_arr[5]]
        #tot_masks = len(mask_arr)
        totiter = 1
        for iter in range(totiter):
            clf()
            if iter == 0:
                if not zonca_sims:
                    vmin, vmax = -80., 80. #None, None
                else:
                    if dust_or_sync == 'dust':
                        vmin, vmax = 0., 100.
                    else:
                        vmin, vmax = None, None
            else:
                vmin, vmax = None, None ##-100., 100. #None, None

            if t_only:
                currmap = [currmap]

            mask_str_arr = []
            if use_s4like_mask_v2:
                mask_str_arr = ['Mask 1: S4-Clean: fsky = ', 'Mask 2: S4-Dirty: fsky = ', 'Mask 3: {\it Planck} rest: fsky = ']

                planck_mask_fname = '%s/HFI_Mask_GalPlane-apo0_2048_R2.00.fits' %(mask_folder)
                planck_mask = H.read_map(planck_mask_fname, verbose = 1, field = (2))#,2,3))
                if (1):
                    planck_mask = H.ud_grade(planck_mask, 128)
                    #simple rotation from gal to celestial
                    planck_mask = healpix_rotate_coords(planck_mask, coord = ['G', 'C'])
                    planck_mask = H.smoothing(np.copy(planck_mask), fwhm = np.radians(10.), verbose = verbose)#, lmax = lmax)
                    planck_mask = H.ud_grade(planck_mask, nside_out = nside)
                    thresh = 0.4
                    planck_mask[planck_mask<thresh] = 0.
                    planck_mask[planck_mask!=0] = 1.                
                planck_mask = H.ud_grade(planck_mask, nside)
                planck_rest_of_sky = (1 - cmbs4_hit_map) * planck_mask
                mask_arr.append( planck_rest_of_sky )
                tot_masks = len(mask_arr)

            if use_spt3g_mask:
                mask_str_arr = ['Winter', 'Summer: el1c/el2c', 'Summer: el1b/el2b', 'Summer: el1-el5']

            vmin, vmax = 0., 60.
            for mask_iter in range(tot_masks):
                fsky = np.mean(mask_arr[mask_iter])
                tit = r'%s @ %s GHz + Mask %s: f$_{\rm sky} = %.2f$' %(dust_or_sync, nuarr[0], mask_iter, fsky)
                if len(mask_str_arr)>0:
                    tit = r'%s @ %s GHz: %s: f$_{\rm sky}$=%.2f' %(dust_or_sync, nuarr[0], mask_str_arr[mask_iter], fsky)
                #H.mollview(currmap[0] * mask_arr[mask_iter], sub = (1,tot_masks,mask_iter+1), title_fontsize = 8, title = tit, min = vmin, max = vmax, cbar=True,); #unit = r'$\mu K$')
                H.mollview(currmap[0] * mask_arr[mask_iter], sub = (2,2,mask_iter+1), title_fontsize = 8, title = tit, min = vmin, max = vmax, cbar=True, unit = r'$\mu K$')
            
            if iter == 0:
                if zonca_sims:
                    plname = '/Users/sraghunathan/Research/SPTPol/analysis/git/DRAFT/scripts/reports/galactic_sims/maps_masks/%s_%s.pdf' %(dust_or_sync, nu)
                else:
                    plname = '/Users/sraghunathan/Research/SPTPol/analysis/git/DRAFT/scripts/reports/galactic_sims/maps_masks/%s_%s_fixedcolourscale.pdf' %(dust_or_sync, nu) 
            else:
                plname = '/Users/sraghunathan/Research/SPTPol/analysis/git/DRAFT/scripts/reports/galactic_sims/maps_masks/%s_%s_freecolourscale.pdf' %(dust_or_sync, nu)
            #savefig(plname)
            if use_lat_step_mask:
                    plname = '/Users/sraghunathan/Research/SPTPol/analysis/git/DRAFT/scripts/reports/galactic_sims/maps_masks/%s_%s_S4wide_cluster_search.pdf' %(dust_or_sync, nu)
            if use_s4like_mask_v2:
                #plname = '/Users/sraghunathan/Research/SPTPol/analysis/git/DRAFT/reports/galactic_sims/maps_masks/%s_%s_masks_S4_v2_Neff_cos_el_%s.png' %(dust_or_sync, nu, cos_el)
                plname = '/Users/sraghunathan/Research/SPTPol/analysis/git/DRAFT/reports/galactic_sims/maps_masks/%s_%s_masks_S4_v2_Neff_cos_el_%s_v2.png' %(dust_or_sync, nu, cos_el)
            savefig(plname, dpi = 150)
            show(); #sys.exit()

        clf()
        cmbs4_hit_map_flist = glob.glob('%s/high_cadence_hits_*_cosecant_modulation.fits' %(cmbs4_footprint_folder))
        for cntr, cmbs4_hit_map_fname in enumerate( sorted( cmbs4_hit_map_flist ) ):
            fname_str = cmbs4_hit_map_fname.split('/')[-1].replace('.fits', '').replace('_', '\_')
            cmbs4_hit_map = H.read_map(cmbs4_hit_map_fname, verbose = verbose)
            cmbs4_hit_map_dummy = np.copy(cmbs4_hit_map)
            cmbs4_hit_map_dummy[cmbs4_hit_map_dummy!=0] = 1.
            fsky = np.mean(cmbs4_hit_map_dummy)
            H.mollview(cmbs4_hit_map, sub = (1,3,cntr+1), title = r'%s: f$_{\rm sky} = %.2f$' %(fname_str, fsky), title_fontsize = 6); 
        plname = '/Users/sraghunathan/Research/SPTPol/analysis/git/DRAFT/reports/galactic_sims/maps_masks/S4_hitmaps.png'
        savefig(plname, dpi = 150)
        #show()
        sys.exit()


    logline = '\tget power spectra now\n'
    lf = open(log_file,'a'); lf.writelines('%s\n' %(logline));lf.close()
    print(logline)
    
    if not os.path.exists(opfname):
        resdic = {}
        resdic['fsky_arr'] = fsky_arr
        resdic['lmax'] = lmax
        resdic['cl_dic'] = {}
    else:
        resdic = np.load(opfname, allow_pickle = 1).item()

    for mask_iter in range(tot_masks):
        if mask_iter not in resdic['cl_dic']:
            resdic['cl_dic'][mask_iter] = {}
        for nu1 in nuarr:
            for nu2 in nuarr:

                print(nu1, nu2)

                logline = '\tMask = %s: (%s,%s)\n' %(mask_iter, nu1, nu2)
                lf = open(log_file,'a'); lf.writelines('%s\n' %(logline));lf.close()
                print(logline)

                if nu1 != 39 and nu2 != 39:
                    if (nu2, nu1) in resdic['cl_dic'][mask_iter] or (nu1, nu2) in resdic['cl_dic'][mask_iter]: 
                        logline = '\t\talready complete\n'
                        lf = open(log_file,'a'); lf.writelines('%s\n' %(logline));lf.close()
                        print(logline)
                        continue

                map1, map2 = map_dic[nu1], map_dic[nu2]

                curr_mask = mask_arr[mask_iter]
                fsky = np.mean(curr_mask)
                map1 = map1 * curr_mask
                map2 = map2 * curr_mask

                if testing:
                    H.mollview(map1, sub = (1,2,1)); H.mollview(map1, sub = (1,2,2)); show()


                curr_cl = H.anafast(map1, map2, lmax = lmax)
                curr_cl /= fsky
                resdic['cl_dic'][mask_iter][(nu1, nu2)] = curr_cl

                if not testing:
                    np.save(opfname, resdic)
    if testing:
        from IPython import embed; embed()
    sys.exit()



#plot
#opfname = 'cl_websky.npy'
from pylab import *
from matplotlib import rc;rc('text', usetex=True);rc('font', weight='bold');matplotlib.rcParams['text.latex.preamble'] = [r'\boldmath']
#mpl.rcParams['hatch.linewidth'] = 0.5
rcParams['font.family'] = 'serif'
rc('text.latex',preamble=r'\usepackage{/Users/sraghunathan/.configs/apjfonts}')

resdic = np.load(opfname, allow_pickle = 1).item()
cl_arr = resdic['cl_arr']
fname_arr = resdic['fname_arr']
lmax = resdic['lmax']

spt_el = 3000.
spt_points_arr = [0, (9.6, 3.46), (66.1, 50.6), 2.9]
spt_points_err_arr = [(0., 0.), (0.36, 0.54), (3.1, 4.4), (1.3, 1.3)]

cl_mod_arr = np.zeros( (len(fname_arr)-1, lmax + 1) )
cl_mod_arr[:3] = np.copy(cl_arr)[:3]
cl_mod_arr[3] = np.sum(cl_arr[4:], axis = 0)


el = np.arange( lmax + 1 )
elnorm = 3000
elind = np.where(el == elnorm)
dl = cl_mod_arr[1] * el * ( el + 1)/ 2 / np.pi    
fac = ( np.sum(spt_points_arr[1]) / dl[elind] )**0.5
print(dl[elind], fac, dl[elind] * fac)

clf()
colorarr = ['navy', 'orangered', 'darkred', 'black']
ax = subplot(111, yscale = 'log')#, xscale = 'log')
for cntr, cl in enumerate( cl_mod_arr[:-1] ):
    if cntr == 0: continue
    lab = fname_arr[cntr].split('/')[-1].replace('_','-').replace('.fits', '').replace('0','').replace('nu','').upper()
    colorval = colorarr[cntr]
    if lab.lower().find('ksz')==-1:
        #lab = r'Unmasked: %s' %(lab) 
        lab = r'Masked + scaled: %s' %(lab) 

    ##cl = cl  * div_arr[cntr]**2.
    cl = ( cl ) * fac**2.
    el = np.arange( len(cl) )
    dl = cl * el * ( el + 1)/ 2 / np.pi    
    plot(el, dl, lw = 1., color = colorval, label = lab)
    spt_data_point = np.sum(spt_points_arr[cntr])
    spt_data_point_err = np.sqrt( spt_points_err_arr[cntr][0]**2. + spt_points_err_arr[cntr][1]**2. )
    print(spt_data_point, spt_data_point_err)
    errorbar(spt_el, spt_data_point, yerr = spt_data_point_err, color = colorval, marker = 'o', capsize = 2.)

    if cntr == 1 or cntr == 2:
        if cntr == 1:
            freq1 = 150
        elif cntr == 2:
            freq1 = 220

        el, cl_dg_po, cl_dg_clus = fg.get_cl_dust(freq1, freq1, freq0 = 150)
        dl_dg_po = cl_dg_po * el * ( el + 1)/ 2 / np.pi    
        dl_dg_clus = cl_dg_clus * el * ( el + 1)/ 2 / np.pi    
        dl_spt = dl_dg_po + dl_dg_clus
        spt_lwval = 2.
        plot(el, dl_spt, lw = spt_lwval, color = colorval,alpha = 0.5)#, label = r'SPT')
        plot(el, dl_dg_po, lw = spt_lwval, ls = '--', color = colorval, alpha = 0.5)#, label = r'SPT: Po')
        plot(el, dl_dg_clus, lw = spt_lwval, ls = ':' , color = colorval, alpha = 0.5)#, label = r'SPT: Clus')

if (0): #add sehgal
    sehgal_dic = np.load('cross_power_CIB_sehgal_SO_sims.fits.npy', allow_pickle = 1).item()
    sehgal_cl_dic = sehgal_dic['cl_dic']
    sehgal_lmax = sehgal_dic['lmax']
    sehgal_cl = sehgal_cl_dic[(217, 217)]

    sehgal_el = np.arange( len(sehgal_cl) )
    sehgal_dl = sehgal_cl * sehgal_el * ( sehgal_el + 1)/ 2 / np.pi    
    plot(sehgal_el, sehgal_dl, lw = 1., color = 'darkred', label = r'Sehgal 217 $\times$ 217', alpha = 0.5)

if (1):
    plot([], [], lw = spt_lwval, color = 'k',alpha = 0.5, label = r'SPT')
    plot([], [], lw = spt_lwval, ls = '--', color = 'k',alpha = 0.5, label = r'SPT: Po')
    plot([], [], lw = spt_lwval, ls = ':' , color = 'k',alpha = 0.5, label = r'SPT: Clus')
    #plot([], [], marker = '.', label = r'Masked: SPT (G15)', ls = 'None')

ylim(0.5, 2e3)
#plot([], [], marker = 'None', label = r'Unmasked Websky')
legend(loc = 2, fancybox = 1, fontsize = 12)
xlabel(r'Multipole $\ell$'); ylabel(r'D$_{\ell}$')
xlim(100., 4000)
show();


