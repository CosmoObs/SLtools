# -*- coding: utf-8 -*-
#!/usr/bin/env python

"""Package to deal with (FITS) image copy/cuts"""

##@package imcp
#
#
# This package contains functions to cutout smaller regions
# from a original image, paste the image into a bigger array,
# as well as any other function for image copy manipulation.
#
# Executable package: YES
#
# To see the package help message, type:
#
# > python imcp.py --help
#
# 'imcp' has one mandatory argument: the (input) image FITS file.
# If none of the options are given, imcp outputs a file named
# 'output.fits' with 1/4th of the original image area, with the 
# same central point/pixel.
#
# > python imcp.py input_file.fits
#
# The options available make possible to define the central point
# where the cut will be centered, as well as the side lengths for
# output image. These values can be passed in 'pixel' or 'degrees'
# units.


import sys;

import os
import pyfits
import string
import commands
import numpy as np;


#============================================================================================
def _hdr_dimpix( hdr ):
	"""Read header params and output image size and degrees-to-pixel convertion factor"""


#	x_img_size = NAXIS1 = int(hdr['NAXIS1']);
#	y_img_size = NAXIS2 = int(hdr['NAXIS2']);

	CUNIT1 = str(hdr['CUNIT1']);
	CUNIT2 = str(hdr['CUNIT2']);

	CD1_1 = float(hdr['CD1_1']);
	CD2_2 = float(hdr['CD2_2']);

	# Convertion from degrees to arcsec units..
	#
	if ( CUNIT1 == 'deg' and CUNIT2 == 'deg' ):
		conv_factor = 3600.0;
	else:
		print >> sys.stderr, "Error: Don't know how to deal with coordinates units of given image.";
		print >> sys.stderr, "       This was written to deal with degrees-to-arcsec transformation.";
		print >> sys.stderr, "       Take a look at 'CUNIT*' parameters on header files and email us.";
		return (False,False);

# Scale is given, using CD1_1, CD1_2, CD2_1, CD2_2 header keys:
#   scale = 3600. * sqrt ((cd11**2+cd21**2+cd12**2+cd22**2)/2.) 

	x_arcsec_to_pixel = abs(CD1_1) * conv_factor;
	y_arcsec_to_pixel = abs(CD2_2) * conv_factor;


	return (y_arcsec_to_pixel, x_arcsec_to_pixel);

# ---

#====================================================================
def _hdr_update(hdr, x_ini, y_ini):
	"""Update header information about world coordinate system"""


	LTV1 = None;
	LTV2 = None;

	NAXIS1 = int(hdr['NAXIS1']);
	NAXIS2 = int(hdr['NAXIS2']);

	CRPIX1 = float(hdr['CRPIX1']);
	CRPIX2 = float(hdr['CRPIX2']);

	CD1_1 = float(hdr['CD1_1']);
	CD2_2 = float(hdr['CD2_2']);

	# Check the edges and update header keyworks for new image..
	#
#	if ( x_ini >= 0  and  y_ini >= 0 ):
#		LTV1 = -1*x_ini;
#		LTV2 = -1*y_ini;
#		CRPIX1 = CRPIX1 + LTV1;
#		CRPIX2 = CRPIX2 + LTV2;
#	if ( x_ini < 0  and  y_ini >= 0 ):
#		LTV2 = -1*y_ini;
#		CRPIX2 = CRPIX2 + LTV2;
#	if ( y_ini < 0  and  x_ini >= 0 ):
#		LTV1 = -1*x_ini;
#		CRPIX1 = CRPIX1 + LTV1;
	LTV1 = -1*x_ini;
	LTV2 = -1*y_ini;
	CRPIX1 = CRPIX1 + LTV1;
	CRPIX2 = CRPIX2 + LTV2;


	# Add some keys to header..
	#
	WCSDIM = 2
	CDELT1  =  CD1_1; # -7.5000002980000E-5
	CDELT2  =  CD2_2; # 7.50000029800001E-5
	LTM1_1 = 1.0
	LTM2_2 = 1.0
	WAT0_001= 'system=image'
	WAT1_001= 'wtype=tan axtype=ra'
	WAT2_001= 'wtype=tan axtype=dec'


	# Header update..
	#
	if (LTV1 != None) :
#		hdr.update('LTV1',LTV1)
		hdr.update('CRPIX1',CRPIX1)
	if (LTV2 != None) :
#		hdr.update('LTV2',LTV2)
		hdr.update('CRPIX2',CRPIX2)
	hdr.update('WCSDIM',WCSDIM)
	hdr.update('CDELT1',CDELT1)
	hdr.update('CDELT2',CDELT2)
	hdr.update('LTM1_1',LTM1_1)
	hdr.update('LTM2_2',LTM2_2)
	hdr.update('WAT0_001',WAT0_001)
	hdr.update('WAT1_001',WAT1_001)
	hdr.update('WAT2_001',WAT2_001)


	return (hdr);

# ---

#===============================================================================================================
def cutout( image, header, coord_unit='pixel', xo=None, yo=None, size_unit='pixel', x_size=None, y_size=None ):
	"""
	Do a snapshot from given fits image.

	cutout( image.array, header [,...] ) -> (ndarray, header)

	The function can deal with 'pixel' and 'degrees' cordinate units. As input, besides the image FITS 
	filename (input_file), it receives a central position (xo,yo) and side lengths (x_size,y_size) for
	output image position and dimensioning. Function returns two values: a numpy.array with resultant 
	image pixel intensities and respective header object.

	If 'xo=None' and 'yo=None', given image (input_file) central pixel will be chosen as (xo,yo).
	If 'x_size=None' and 'y_size=None', half length of each side will be used for output dimensions.

	Input:
	 - image      : Image (FITS) ndarray
	 - header     : Image (FITS) header object
	 - coord_unit : 'pixel' or 'degrees' for position (xo,yo) values
	 - xo         : Horizontal central position for output (cut) image
	 - yo         : Vertical central position for output (cut) image
	 - size_unit  : 'pixel' or 'degrees' for size (x_size,y_size) values
	 - x_size     : Horizontal size (in pixels) of output image
	 - y_size     : Vertical size (in pixels) of output image

	Output:
	 - (ndarray, header) : image array and header object

	"""

	# Initialize some variables..
	#
	x_diff = 0;   x_edge = 0;
	y_diff = 0;   y_edge = 0;
	x_fin = 0;   x_ini = 0;
	y_fin = 0;   y_ini = 0;


	# Read fits file..
	#
#	imagem, hdr = pyfits.getdata( input_file, header=True );
	imagem = image;
	hdr = header;

	if (hdr != None):
		x_arcsec_to_pixel, y_arcsec_to_pixel = _hdr_dimpix( hdr );

	y_img_size, x_img_size = imagem.shape;


	# Verify the type of 'size_unit' param..
	#
	if ( size_unit == 'pixel' ):
		if ( x_size != None ):
			x_cut_size = int(x_size);
		if ( y_size != None ):
			y_cut_size = int(y_size);

	elif ( size_unit == 'degrees' ):
		if ( x_size != None ):
			x_cut_size = int(float(x_size)/x_arcsec_to_pixel);
		if ( y_size != None ):
			y_cut_size = int(float(y_size)/y_arcsec_to_pixel);

	else:
		print >> sys.stderr, "Error: Side sizes type (degrees/pixel) not correctly defined."
		return (False,False);


	# And if no side size was given, define a default value correspondig to half of original image..
	#
	if ( x_size == None ):
		x_cut_size = int(x_img_size/2);
		print >> sys.stdout, "Warning: 'x_size' not given. Using half of image x side(%d)" % (x_cut_size);

	if ( y_size == None ):
		y_cut_size = int(y_img_size/2);
		print >> sys.stdout, "Warning: 'y_size' not given. Using half of image y side(%d)" % (x_cut_size);


	# Check requested output.vs.input image sizes..
	#
	if ( x_cut_size == x_img_size and y_cut_size == y_img_size ):
		print >> sys.stdout, "Warning: Requested output sizes are the same as input image. Returning original image and header.";
		return (imagem,hdr);


	imagemnova = np.zeros( (y_cut_size,x_cut_size), dtype=float );


	# Verify central coordinates values..
	#
	if ( coord_unit == 'pixel' and (xo != None and yo != None) ):
		x_halo = int(xo)-1;
		y_halo = int(yo)-1;

	elif ( coord_unit == 'degrees' and (xo != None and yo != None) ):
		_out = commands.getoutput( 'sky2xy %s %s %s' % (input_file, xo, yo) );
		_out = string.split(_out);
		x_halo = int(_out[-2])-1;
		y_halo = int(_out[-1])-1;

	elif ( xo == None and yo == None ):
		x_halo = int(x_img_size/2);
		y_halo = int(y_img_size/2);
		print >> sys.stdout, "Warning: No central coordinates were given for snapshot.";
		print >> sys.stdout, "         Using image central pixel as the cut center.";
	else:
		print >> sys.stderr, "Error: Central positioning is out of valid values. Exiting.";
		return (False,False);


	# Define the images (in/out) slices to be copied..
	#
	x_ini = x_halo - int(x_cut_size/2) -1;   x_fin = x_ini + x_cut_size;
	y_ini = y_halo - int(y_cut_size/2) -1;   y_fin = y_ini + y_cut_size;

	x_ini_old = max( 0, x_ini );   x_fin_old = min( x_img_size, x_fin );
	y_ini_old = max( 0, y_ini );   y_fin_old = min( y_img_size, y_fin );

	x_ini_new = abs( min( 0, x_ini ));   x_fin_new = x_cut_size - (x_fin - x_fin_old);
	y_ini_new = abs( min( 0, y_ini ));   y_fin_new = y_cut_size - (y_fin - y_fin_old);

	if (hdr != None):
		hdr = _hdr_update(hdr, x_ini, y_ini);

	imagemnova[ y_ini_new:y_fin_new, x_ini_new:x_fin_new ] = imagem[ y_ini_old:y_fin_old, x_ini_old:x_fin_old ];


	return (imagemnova, hdr);

# ---

#=======================================================================================
def poststamp( image, header=None, output_size=None, increase=0, relative_increase=False ):
    """
    Create poststamp from given image array

    poststamp( image.array [,...] ) -> image_out.array

    Function resize the given 'image' array according to image features. The region of interest
    is defined as the largest rectangular region containing non-null pixels, this region will be 
    called "object". Resizing will be done based on the object dimensions.
    Resizing can be done in three different ways:
     i)   Fixing the returned image size by 'output_size' argument.### A single value (for square output)
          or a tuple (Num-raws,Num-columns) can be passed for fixed output size:
          Use:
              output_size != None

     ii)  Fixing the increase in number of pixels for output image border
          Use:
              increase != 0
              relative_increase = False

     iii) Defining the border size by a relative factor based on the object dimensions
          Use:
              increase != 0
              relative_increase = True

    If just 'image' array is passed, function will return an image array with just the region of interest
    pixels on it. No increase is done.

    Input:
     - image             : 2D ndarray
     - header            : If there is a FITS header object related to 'image', give it here to be updated
     - output_size       : A single value (for square output) or a tuple (Num_raws, Num_columns) can be passed
                           If given, none of following parameters will be used
     - increase          : Absolute or relative factor. 'relative_increase' defines if it is a fixed number 
                           of pixels or relative to image object size
     - relative_increase : If True, 'increase' will be used as a multiplicative factor of object size
                           If False, 'increase' will be the border size over object pixels

    Output:
     - (ndarray,header)  : Resized image array and header instance

    """

    # If output image is required to have a fixed size, do that and finish the funtion..
    #
    if ( output_size != None ):
        try:
            y_size, x_size = output_size;
        except:
            y_size = x_size = output_size;

        image_out, header = cutout( image, header, x_size=x_size, y_size=y_size );

        return (image_out,header);


    # Otherwise, if the output has variable size, according to feature dimensions..
    #
    y_inds, x_inds = np.where(image != 0);

    y_ind_max = max(y_inds);
    y_ind_min = min(y_inds);

    x_ind_max = max(x_inds);
    x_ind_min = min(x_inds);
    
    if ( increase != 0 ):

        x_size = (x_ind_max - x_ind_min)/2.;
        y_size = (y_ind_max - y_ind_min)/2.;
        xo = x_size + x_ind_min;
        yo = y_size + y_ind_min;

        if (relative_increase == True):

            x_size = x_size*increase;
            y_size = y_size*increase;

        else:

            x_size = x_size + 2*increase;
            y_size = y_size + 2*increase;


        image_out, header = cutout( image, header, xo=xo, yo=yo, x_size=x_size, y_size=y_size );

    else:

        image_out = image;


    return (image_out, header);

# ---

#======================================================================
def _write_fitsimg( outfile, image, header ):
	"""Just write down the image array and header to a FITS file"""
	os.system('[ -f %s ] && rm %s' % (outfile,outfile));
	pyfits.writeto( outfile, image, header );
	return (True);
# ---


# @cond
#==========================
if __name__ == "__main__" :

	from optparse import OptionParser;

	usage="Usage: %prog [options] <image_fits_file>"

	parser = OptionParser(usage=usage);

#	parser.add_option('-f',
#			  dest='input_file', default=None,
#			  help='Input FITS image');
	parser.add_option('-o',
			  dest='output_file', default='output.fits',
			  help='Output FITS image [output.fits]');
	parser.add_option('--coord',
			  dest='coord_unit', default='pixel',
			  help='Coordinate units, for xo & yo [pixel]');
	parser.add_option('--xo',
			  dest='xo', default=None,
			  help='X centre');
	parser.add_option('--yo',
			  dest='yo', default=None,
			  help='Y centre');
	parser.add_option('--size',
			  dest='size_unit', default='pixel',
			  help='Side length untis, for x_size & y_size [pixel]');
	parser.add_option('--x_size',
			  dest='x_size', default=None,
			  help='Size of x side on output');
	parser.add_option('--y_size',
			  dest='y_size', default=None,
			  help='Size of y side on output');
	parser.add_option('--increase',
			  dest='incr', default=0,
			  help='Amount of border/enlarge of image. It is an additive or multiplicative factor. See "relative_increase" flag');
	parser.add_option('--relative_increase', action='store_true'
			  dest='rel_incr', default=False,
			  help='Turn on the multiplicative factor for increase image.')

	(opts,args) = parser.parse_args();

	outfits = opts.output_file;
	coord = opts.coord_unit;
	x = opts.xo;
	y = opts.yo;
	size = opts.size_unit;
	dx = opts.x_size;
	dy = opts.y_size;
	incr = opts.incr;
	rel_incr = opts.rel_incr;

	if ( len(args) != 1 ):
		parser.error("Wrong number of arguments. Try option '--help' for usage and help messages.")

	infits = args[0];

	image, hdr = pyfits.getdata(infits,header=True);

	if (incr == 0):
		imgcut, hdr = cutout( image, header=hdr, coord_unit=coord, xo=x, yo=y, size_unit=size, x_size=dx, y_size=dy );
	else:
		imgcut, hdr = poststamp( image, header=hdr, increase=incr, relative_increase=rel_incr )

	_write_fitsimg( outfits, imgcut, hdr );

	sys.exit(0);

# ------------------
# @endcond
