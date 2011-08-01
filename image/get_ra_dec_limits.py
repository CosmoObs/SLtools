
##@package get_image_limits
#
# This function reads RA and DEC values correponding to the spatial limits of the sky image. 
# 
#@param image_file, the name of image(tile) fits file.
#@return ra,dec, two tuples with image sky region limits.

# ---
import pyfits;
import pywcs
import numpy

def get_image_limits(image):

    hdr = pyfits.getheader( image)
    
    wcs = pywcs.WCS(hdr) 

    naxis1 = hdr['NAXIS1']
    naxis2 = hdr['NAXIS2']
    ra_centre = hdr['CRVAL1']
    dec_centre = hdr['CRVAL2']
    pixcrd = numpy.array([[1,1],[naxis1,naxis2]])

    skycrd = wcs.wcs_pix2sky(pixcrd,1)

    ra_min = skycrd[0][0]
    dec_min = skycrd[0][1]
    ra_max = skycrd[1][0]
    dec_max = skycrd[1][1]


    return ((ra_min,ra_max),(dec_min,dec_max))



