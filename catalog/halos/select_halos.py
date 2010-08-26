import sys;
import numpy as np;
from sltools.coordinate import *;
from sltools.image  import *;
from sltools.catalog  import *;

def select_halos( hdulist, minimum_halo_mass=0, ra=(), dec=(), image_file='' ): 
    """Select halos with mass greater than a given bottom value.

    select_halos( hdulist, minimum_halo_mass, ra=(), dec=(), image_file='' )

    Input:
     - hdulist : pyfits catalog object
     - minimum_halo_mass : bottom value for Halos selection
     - ra - (RA_min,RA_max) : Right Ascension limiting values for halos selection
     - dec - (DEC_min,DEC_max) : Declination limiting values for halos selection
     - image_file : Sky image filename

    Output:
     - list with selected Halo IDs

    """


    # First, let us verify the passed arguments..
    #
    if ( len(ra)==1  or  len(dec)==1 ):
        print >> sys.stderr,"ra/dec should be a tuple with minimum and maximum limit values";
        return (False);

    if ( image_file==''  and  ra==() and dec==() ):
        print >> sys.stderr,"No sky region information were given. Try to pass 'image_file' or ra and dec";
        return (False);

    elif ( image_file=='' ):
        ra,dec = readout_radec( ra, dec );
        ra_min, ra_max = ra; 
        dec_min, dec_max = dec;

    else: 
        ra, dec = get_image_limits( image_file );
        ra_min, ra_max = ra; 
        dec_min, dec_max = dec; 


    # Read out some columns of interest from halo catalogue..
    #
    dic = get_fits_data( hdulist, 'HALOID', 'RA', 'DEC', 'M200' );
    if not (dic):
        print >> sys.stderr, "Error: empty output from get_fits_data().";
        return (False);


    haloids = dic['HALOID'];
    ra = dic['RA'];
    dec = dic['DEC'];
    m200 = dic['M200'];


    select = ( ra > ra_min) & ( ra < ra_max ) & ( dec > dec_min) & (dec < dec_max) & (m200 > minimum_halo_mass) 

    index = np.where(select)[0]

    return  haloids[np.where(select)[0]] 

