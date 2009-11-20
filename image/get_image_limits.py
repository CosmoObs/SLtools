


##@package get_image_limits [formerly RaDecLimits]

# This function reads RA and DEC values correponding to the spatial limits of the sky image. 
# 
#@param image_file, the name of image(tile) fits file.
#@return ra,dec, two tuples with image sky region limits.

def get_image_limits( image_file ):

        # Image Header..
        hdr = pyfits.getheader( ImageFitsFile )

        # Some image coordinates and measurements..
        naxis1 = hdr['NAXIS1']
        naxis2 = hdr['NAXIS2']
        ra_centre = hdr['CRVAL1']
        dec_centre = hdr['CRVAL2']

        # Pixel (0,0) to sky (degrees) coordinates..
        _out = string.split(_out)
        ra_min = string.atof(_out[0])
        dec_min = string.atof(_out[1])

        # Now, pixel (N,N) to sky coordinates..
        _out = string.split(_out)
        ra_max = string.atof(_out[0])
        dec_max = string.atof(_out[1])

        # Function returns sky coordinates in degrees..
        return ((ra_min,ra_max),(dec_min,dec_max)); 

# -----


