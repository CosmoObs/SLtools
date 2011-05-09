""" Module to deal with tables matching """

import pyfits
import glob
import os

from sltools.coordinate.wcs_conversion import xy2radec

def get_image_limits(image_path):
    """ Function get the image limits in RA and DEC
        
    
    
    Input:
     - image_path : a path to the image
    
    Output:
     - ra_min, dec_min, ra_max, dec_max

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

def get_subdirectories(directory_path):
    """ Function to get all subdirectories in a directory
        
      
    Input:
     - directory_pat: a path to the directory
    
    Output:
     - [directory1, directory2 ...]

    ---
    """
    folders = os.listdir(directory_path)
    folders.sort()
    return folders


def make_sdss_query():
    """ Function to criate a query for the SDSS CasJob system.
        
    
    
    Input:
     - none
    
    Output:
     - none

    ---
    """

    im_path = '/home/gbcaminha/SOGRAS/2008b_data/2008b_reduced/'

    folders = os.listdir(im_path)
    folders.sort()

    image_files = []

    for i in folders:
        image_files.append( glob.glob( os.path.join(im_path+i, '*r.fits'))[0] )


    radec_limits = []

    for i in range(len(image_files)):
        radec_limits.append(get_image_limits(image_files[i]))
        str_cond1 = 'WHERE g.ra>' + str(radec_limits[i][0])
        str_cond2 = 'AND g.ra<' + str(radec_limits[i][2])
        str_cond3 = 'AND g.dec>' + str(radec_limits[i][1])
        str_cond4 = 'AND g.dec<' + str(radec_limits[i][3])
        str_cond5 = 'AND g.objID=p.objID'

        select1 = 'SELECT g.objID, g.ra, g.raErr, g.dec, g.decErr, p.z, p.zErr,'
        select2 = ' g.modelMag_u, g.modelMag_g, g.modelMag_r, g.modelMag_i,'
        select3 = ' g.modelMag_z, g.modelMagErr_u, g.modelMagErr_g,'
        select4 = ' g.modelMagErr_r, g.modelMagErr_i, g.modelMagErr_z'

        into1 = 'into mydb.SOGRAS_' + str(folders[i])

        #from1 = 'from Galaxy g, Photoz p'

        print select1, select2, select3, select4, into1
        print 'from Galaxy g, Photoz p' #from1
        print str_cond1, str_cond2, str_cond3, str_cond4, str_cond5
        print
        print '**********************************************************'
        print 


#Main Program


