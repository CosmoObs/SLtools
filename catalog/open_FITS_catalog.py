

##@package open_FITS_catalog [formerly OpenFITSCatalog]

# Open FITS catalog to HDUlist.

#@param catalog_file name of FITS file 
#@return hdulist list of FITS file slices

def open_FITS_catalog( catalog_file ):
    import pyfits; 
    import string; 

    if ( string.find(catalog_file,'.fit') == -1 ): 
        return False; 

    # Open fits table
    hdulist = pyfits.open(catalog_filename)

    return (hdulist); 

# ----- 
 

