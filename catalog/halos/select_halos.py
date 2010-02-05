

##@package select_halos [formerly SelectHalos]
# Select halos with mass greater than minimum_halo_mass. 

#
# Example:
# After reading a halo catalog, get the mass of halo id 234

#

# First, read the halo catalog

#
# selected_halos=select_halo(halo_catalog, minimum_halo_mass, sky_region)
#
# Second, find the index of halo id 234
#
# index=selected_halos['id'].index(234) 
# 
# then return the mass

#

# m200=selected_halos[m200][index]

# 
#@param HDUlist : pyfits catalog
#@param minimum_halo_mass
#@param ra : a tuple (ra_min,ra_max) with Right Ascension limit values for sky region
#@param dec : a tuple (dec_min,dec_max) with Declination limit values for sky region
#@param image_file
#@return HaloIDs : a list containing all Halo IDs above the mass criterium and within given sky region 
def select_halos( hdulist, minimum_halo_mass, ra=(), dec=(), image_file='' ): 
    import sys;
    import tools;
    import tools.catalog;
    from tools.catalog.halos.get_catalog_data import get_catalog_data;
    from tools.image.get_image_limits import get_image_limits;

    if ( len(ra)==1  or  len(dec)==1 ):
        return (False);

    if ( image_file==''  and  ra==() and dec==() ): 
        return (False);

    elif ( image_file=='' ): 
        ra_min, ra_max = ra; 
        dec_min, dec_max = dec;

    else: 
        ra, dec = get_image_limits( image_file );
        ra_min, ra_max = ra; 
        dec_min, dec_max = dec; 

    dic = get_catalog_data( hdulist, 'HALOID', 'RA', 'DEC', 'M200' );

    haloid = dic['HALOID'];
    ra = dic['RA'];
    dec = dic['DEC'];
    m200 = dic['M200'];

    ids = [];
    controle = [];
    for i in range( 0, len( haloid )):

        # Check whether the current HaloID was already used..
        if ( haloid[i] in controle or haloid[i] == 0 ):
            continue;

        # Select Halos by Sky region..
        if ( ra_min <= ra[i] < ra_max  and  dec_min <= dec[i] < dec_max ):
            if ( minimum_halo_mass <= m200[i] ):
                ids.append( haloid[i] );

        controle.append( haloid[i] );

    # Return the big ones
    return (ids);


# ----- 


