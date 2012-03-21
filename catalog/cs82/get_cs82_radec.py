from __future__ import division
from ast import literal_eval as lit
import sys
import os
import errno
import numpy as np
import math
import cmath
from scipy import optimize
import collections
import fit_mag as fmg
import glob
import pyfits
import matplotlib.pyplot as pl
import plot_templates as pltemp
from time import strftime


sys.path.append("/home/brunomor/Documents/Trabalho/Repositorio") ### Append Bruno
from sltools.catalog import fits_data as fd ### Import Bruno
from sltools.catalog import ascii_data as ascd ### Import Bruno
from sltools.catalog import table_matching as tbm
from sltools.coordinate import wcs_conversion as wcscnv

########################################
### SETTING FOLDERS, FILES, ETC   ######
########################################

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST:
            pass
        else: raise


input_cs82 = sys.argv[1]

folder_path = os.path.dirname(input_cs82).replace('/AFTER_PSFEX','')



# =============================== Getting CS82 Catalog Information

  
cs82_tbhdu = pyfits.open(input_cs82,ignore_missing_end=True,memmap=True)[2]

cs82_cut_hdu = fd.sample_entries(cs82_tbhdu,CLASS_STAR=(0,0.95), MAG_AUTO=(0,21), MAGERR_AUTO=(0,0.2172))

cs82_data = cs82_cut_hdu.data


cs82_ra = cs82_data.field('alpha_J2000')
cs82_dec = cs82_data.field('delta_J2000')

results = np.column_stack((cs82_ra,cs82_dec))

np.savetxt(folder_path + '/sdss_radec_W4.txt', results, fmt='%.3f    %.5f')

