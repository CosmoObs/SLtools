##@package select_halos
# Select halos with mass greater than a given bottom value.
# 
#@param HDUlist : pyfits catalog object
#@param minimum_halo_mass : bottom value for Halos mass in Msun units
#@param ra : a tuple (ra_min,ra_max) with Right Ascension limit values for sky region
#@param dec : a tuple (dec_min,dec_max) with Declination limit values for sky region
#@param image_file : sky image/tile fits filename
#@return HaloIDs : a list containing all Halo IDs above the mass criterium and within given sky region 

import sys;

# AddArcs/SLtools..
from tools.catalog.halos import get_catalog_data;
from tools.image import get_image_limits;

# AddArcs/Pipeline..
from functions import readout_radec;


#####################################
#
def select_halos( hdulist, minimum_halo_mass, ra=(), dec=(), image_file='' ): 
    """Select halos with mass greater than a given bottom value.

    select_halos( hdulist, minimum_halo_mass, ra=(), dec=(), image_file='' )
    select_halos( ) ---> list with selected Halo IDs.

    Input:
    hdulist : pyfits catalog object
    minimum_halo_mass : bottom value for Halos selection
    ra - (RA_min,RA_max) : Right Ascension limiting values for halos selection
    dec - (DEC_min,DEC_max) : Declination limiting values for halos selection
    image_file : Sky image filename

    Output:
    List with selected Halo IDs
    """


    # First, let us verify the passed arguments..
    #
    if ( len(ra)==1  or  len(dec)==1 ):
        print >> sys.stderr,"ra/dec should be a tuple with minimum and maximum limiting values";
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

    del ra,dec;

    # Read out some columns of interest from halo catalogue..
    #
    dic = get_catalog_data( hdulist, 'HALOID', 'RA', 'DEC', 'M200' );

    haloid = dic['HALOID'];
    _ra = dic['RA'];
    _dec = dic['DEC'];
    m200 = dic['M200'];

    ids = [];
    controle = [];
    for i in range( 0, len( haloid )):

        # Check whether the current HaloID was already used..
        if ( haloid[i] in controle or haloid[i] == 0 ):
            continue;

        # Select Halos by Sky region..
        if ( ra_min <= _ra[i] < ra_max  and  dec_min <= _dec[i] < dec_max ):
            if ( minimum_halo_mass <= m200[i] ):
                ids.append( haloid[i] );

        controle.append( haloid[i] );

    # Return the big ones
    return (ids);


# ----- 


