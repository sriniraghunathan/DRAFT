import pysm3
import pysm3.units as u
import healpy as H
import numpy as np

import os
##os.environ['KMP_WARNINGS'] = 'off'

import matplotlib
matplotlib.use('Agg')
from pylab import *

import warnings
warnings.filterwarnings("ignore")

band=145
nside = 2048
dust_arr = ['d1', 'd2', 'd3', 'd4']
op_fd = 'sims/pysm/'
if not os.path.exists( op_fd ): os.system('mkdir -p %s' %(op_fd))
for dcntr, d in enumerate( dust_arr ):
    fname = '%s/%s.fits' %(op_fd, d)
    if not os.path.exists( fname ):
        sky = pysm3.Sky(nside=nside, preset_strings=[d])

        print(sky.components)

        map_145GHz_arr = sky.get_emission(band * u.GHz)
        map_145GHz_arr = map_145GHz_arr.to(u.uK_CMB, equivalencies=u.cmb_equivalencies(band*u.GHz))
        print(len(map_145GHz_arr))

        H.write_map( fname, map_145GHz_arr[0])
    else:

        hmap = H.read_map( fname )

        clf()
        H.mollview(hmap, title = r'Dust at %s GHz: pySM3: Model = %s' %(band, d), min = 0., max = 150.)
        plname = 'plots/pysm_gal_dust_%s.png' %(d)
        savefig(plname)
        close()