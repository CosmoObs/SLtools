
##@package get_image_limits [formerly RaDecLimits]

# This function reads RA and DEC values correponding to the spatial limits of the sky image. 
# 
#@param image_file, the name of image(tile) fits file.
#@return ra,dec, two tuples with image sky region limits.

def get_image_limits( image_file ):
	import pyfits;
	import commands;

        # Image Header..
        hdr = pyfits.getheader( image_file )

        # Some image coordinates and measurements..
        naxis1 = hdr['NAXIS1']
        naxis2 = hdr['NAXIS2']
        ra_centre = hdr['CRVAL1']
        dec_centre = hdr['CRVAL2']

        # Pixels (0,0) and (10000,10000) to sky (degrees) coordinates..
	out = commands.getoutput( "xy2sky -d %s %s %s %s %s" % ( image_file, 0, 0, 10000, 10000 ));
        out_min,out_max = out.split('\n');
	out_min = out_min.split(None);
	out_max = out_max.split(None);
        ra_min = float(out_min[0]);
        dec_min = float(out_min[1]);
	ra_max = float(out_max[0]);
	dec_max = float(out_max[1]);

        # Function returns sky coordinates in degrees..
        return ((ra_min,ra_max),(dec_min,dec_max)); 

# -----


