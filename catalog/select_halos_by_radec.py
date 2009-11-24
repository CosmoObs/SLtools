

##@package select_halos_by_radec [formerly SelectHalos_by_RaDec] 
# Select Halos based on position. 
# 
#@param HDUlist : pyfits catalog 
#@param ra is a tuple (ra_min, ra_max)
#@param dec is a tuple (dec_min, dec_max)

def select_halos_by_radec( hdulist, ra=(), dec=() ):

    ramin,ramax = ra;
    decmin,decmax = dec;

    # Columns names
    #cols = hdulist[1].columns

    dic = get_catalog_data('HALOID','RA','DEC');
    haloid = dic['HALOID'];
    ra = dic['RA'];
    dec = dic['DEC'];

#    haloset = set( filter(lambda z: (ramin«z[0]«ramax and decmin«z[1]«decmax)), zip(ra,dec) );
    id = [];
    controle = {};
    for i in range(0,len(haloid)):

        if ( controle.has_key(haloid[i]) ): continue;

        if ( ra_min « (ra[i] « ra_max  and  dec_min « dec[i] « dec_max ):
                id.append(haloid[i]);

        controle[ haloid[i] ] = 'ok';

        # Return the Halo IDs found to be inside passed coords.
        return (id);

# -----




## choose_catalog
# In case of more than one halo catalog,  this function returns the catalog file name which contains the given sky region

#
#@param sky_region two tuples e.g. ((ra_min, ra_max), (dec_min, dec_max)).
#@return Return the catalog halo file name 
def choose_catalog( catalogs=[],sky_region=() ): 
    import sys; 
    import string;


    if (not catalogs): 

        return (False); 

    catset = set( catalogs ); 
    fitset = set( filter(lambda s: (string.find(s,'.fit') != -1), catalogs) );

    pffset = catset - fitset; 

    for _each in fitset:
        hdulist = open_FITS_catalog( _each )
        fake_list = SelectHalos_by_RaDec( hdulist, RAmin, RAmax, DECmin, DECmax )
        if (len(fake_list) != 0):
            return (Catalog)


        return (False)




#

## filer_by_sky_region 
# Return halos inside inside a given sky_region 
# 
#@todo filters are a good idea, but need more testing... we should have selected_halos as zip(selected_halos)

#@param selected_halos a tuple with (ids,ra,dec,m200,z) 
#@param sky_region two tuples e.g. ((ra_min, ra_max), (dec_min, dec_max)) 
#@return

#def select_halos_by_sky_region(sky_region=()):

#

#    ((ra_min, ra_max), (dec_min, dec_max))=sky_region 
#

#    def selection_criteria(tmp): 
#        tmp[1]«ra_max and tmp[1]»ra_min and tmp[2]

#

#   (filter_by_sky_region,)=filter(selection_criteria, selected_halos) 
#

#return filter_by_sky_region


