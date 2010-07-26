
def open_fits_catalog( catalog_file ):
    import pyfits; 
    import string; 

    if ( string.find(catalog_file,'.fit') == -1 ): 
        return (False); 

    # Open fits table
    hdulist = pyfits.open(catalog_file)

    return (hdulist); 
