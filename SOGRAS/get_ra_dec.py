""" Module to deal with tables matching """

import pyfits
import glob
import os

from sltools.coordinate.wcs_conversion import xy2radec

def get_image_limits(image_path):
    """ Function get the image limits in RA and DEC
        
    FIXME
    
    Input:
     - image_path : a path to the image
    
    Output:
     - [(index_BtoA0,distance_BtoA0), ]

    ---
    """

    image_hdulist = pyfits.open(image_path)
    image_header = image_hdulist[0].header

    n_axis1 = image_header['NAXIS1']
    n_axis2 = image_header['NAXIS2']

    radec_1 = xy2radec(image_path, 0, 0)
    radec_2 = xy2radec(image_path, n_axis1, n_axis2)

    ra_min = min( [radec_1[0], radec_2[0]] )
    ra_max = max( [radec_1[0], radec_2[0]] )
    dec_min = min( [radec_1[1], radec_2[1]] )
    dec_max = max( [radec_1[1], radec_2[1]] )

    #print ra_min, dec_min
    #print ra_max, dec_max




    return ra_min, dec_min, ra_max, dec_max
#Main Program

IM_PATH = '/home/gbcaminha/SOGRAS/2008b_data/2008b_reduced/'





FOLDERS = os.listdir(IM_PATH)
FOLDERS.sort()

IMAGE_FILES = []

for i in FOLDERS:
    IMAGE_FILES.append( glob.glob( os.path.join(IM_PATH+i, '*r.fits'))[0] )


RADEC_LIMITS = []

for i in range(len(IMAGE_FILES)):
    RADEC_LIMITS.append(get_image_limits(IMAGE_FILES[i]))
    STR_COND1 = 'WHERE g.ra>' + str(RADEC_LIMITS[i][0])
    STR_COND2 = 'AND g.ra<' + str(RADEC_LIMITS[i][2])
    STR_COND3 = 'AND g.dec>' + str(RADEC_LIMITS[i][1])
    STR_COND4 = 'AND g.dec<' + str(RADEC_LIMITS[i][3])
    STR_COND5 = 'AND g.objID=p.objID'

    SELECT1 = 'SELECT g.objID, g.ra, g.raErr, g.dec, g.decErr, p.z, p.zErr, '
    SELECT2 = 'g.modelMag_u, g.modelMag_g, g.modelMag_r, g.modelMag_i, '
    SELECT3 = 'g.modelMag_z, g.modelMagErr_u, g.modelMagErr_g,'
    SELECT4 = 'g.modelMagErr_r, g.modelMagErr_i, g.modelMagErr_z'

    INTO1 = 'into mydb.SOGRAS_' + str(FOLDERS[i])

    FROM1 = 'FROM Galaxy g, Photoz p'

    print SELECT1, SELECT2, SELECT3, SELECT4, INTO1
    print FROM1
    print STR_COND1, STR_COND2, STR_COND3, STR_COND4, STR_COND5
    print
    print '**********************************************************'
    print 

#RADEC_LIMITS[i][0], RADEC_LIMITS[i][1], RADEC_LIMITS[i][2], RADEC_LIMITS[i][3]







