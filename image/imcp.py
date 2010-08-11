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


# \cond
# INTERNAL FUNCTIONS #
#----------------------
def _hdr_dimpix( hdr ):
	"""Read header params and output image size and degrees-to-pixel convertion factor"""


#	x_img_size = NAXIS1 = int(hdr['NAXIS1']);
#	y_img_size = NAXIS2 = int(hdr['NAXIS2']);

	try:
		CUNIT1 = str(hdr['CUNIT1']);
		CUNIT2 = str(hdr['CUNIT2']);
	except:
		print >> sys.stderr, "Warning: Unable to find 'CUNIT[1|2]' parameters in header instance for image coordinates unit.";
		print >> sys.stderr, "         Degrees ('deg') is being used as default value.";
		print >> sys.stderr, "         Take a look inside image header parameters for coordinate units and email us so we can improve the code.";
		CUNIT1 = 'deg';
		CUNIT2 = 'deg';

	CD1_1 = float(hdr['CD1_1']);
	CD2_2 = float(hdr['CD2_2']);

	# Convertion from degrees to arcsec units..
	#
	if ( CUNIT1 == 'deg' and CUNIT2 == 'deg' ):
		conv_factor = 3600.0;

# Scale is given, using CD1_1, CD1_2, CD2_1, CD2_2 header keys:
#   scale = 3600. * sqrt ((cd11**2+cd21**2+cd12**2+cd22**2)/2.) 

	x_arcsec_to_pixel = abs(CD1_1) * conv_factor;
	y_arcsec_to_pixel = abs(CD2_2) * conv_factor;


	return (y_arcsec_to_pixel, x_arcsec_to_pixel);

# ---

#----------------------------------
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
		hdr.update('LTV1',LTV1)
		hdr.update('CRPIX1',CRPIX1)
	if (LTV2 != None) :
		hdr.update('LTV2',LTV2)
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
# \endcond


# =============================================================================================================================
def cutout( image, header=None, coord_unit='pixel', xo=None, yo=None, size_unit='pixel', x_size=None, y_size=None, mask=None ):
	"""
	Do a snapshot from given fits image.

	cutout( image [,...] ) -> (ndarray, header)

	The function can deal with 'pixel' and 'degrees' cordinate units. As input, besides the image FITS 
	filename (input_file), it receives a central position (xo,yo) and side lengths (x_size,y_size) for
	output image position and dimensioning. Function returns two values: a numpy.array with resultant 
	image pixel intensities and respective header object.

	In addition, an image mask information can be given to use as interest region selection. 'mask' is
	expected to be a numpy.where like structure (i.e, mask=numpy.where()). If given, returned image will
	have all pixels null except the ones listed in 'mask' parameter.

	If 'xo=None' and 'yo=None', given image (input_file) central pixel will be chosen as (xo,yo).
	If 'x_size=None' and 'y_size=None', half length of each side will be used for output dimensions.

	Input:
	 - image      : Image (numpy ndarray)
	 - header     : Image header object
	 - coord_unit : 'pixel' or 'degrees' for position (xo,yo) values
	 - xo         : Horizontal central position for output (cut) image
	 - yo         : Vertical central position for output (cut) image
	 - size_unit  : 'pixel' or 'degrees' for size (x_size,y_size) values
	 - x_size     : Horizontal size (in pixels) of output image
	 - y_size     : Vertical size (in pixels) of output image
	 - mask       : tuple with arrays (y,x) with image/array indexes of interest

	Output:
	 - (ndarray, header) : resultant image array and (updated) header instance

	"""


	imagem = image;
	hdr = header;

	# Initialize some variables..
	#
	x_diff = 0;   x_edge = 0;
	y_diff = 0;   y_edge = 0;
	x_fin = 0;   x_ini = 0;
	y_fin = 0;   y_ini = 0;


	if ( hdr ):
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
		print >> sys.stderr, "Error: Side sizes type (unit = degrees|pixel) not correctly defined."
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
		print >> sys.stdout, "Warning: Requested output sizes are the same as input image. Returning image and header unchanged.";
		return (imagem,hdr);


	# Verify central coordinates values..
	#
	if ( coord_unit == 'pixel' and (xo != None and yo != None) ):
		x_halo = int(xo);
		y_halo = int(yo);

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
	x_ini = x_halo - int(x_cut_size/2) #-1;
	x_fin = x_ini + x_cut_size;
	y_ini = y_halo - int(y_cut_size/2) #-1;
	y_fin = y_ini + y_cut_size;

	x_ini_old = max( 0, x_ini );   x_fin_old = min( x_img_size, x_fin );
	y_ini_old = max( 0, y_ini );   y_fin_old = min( y_img_size, y_fin );

	x_ini_new = abs( min( 0, x_ini ));   x_fin_new = x_cut_size - (x_fin - x_fin_old);
	y_ini_new = abs( min( 0, y_ini ));   y_fin_new = y_cut_size - (y_fin - y_fin_old);

	# If header, update pixel<->sky information..
	#
	if ( hdr ):
		hdr = _hdr_update(hdr, x_ini, y_ini);

	# Initialize new image, and take all index list..
	#
	imagemnova = np.zeros( (y_cut_size,x_cut_size), dtype=imagem.dtype );
	ind_z = np.where(imagemnova == 0);

	# Copy requested image slice..
	#
	imagemnova[ y_ini_new:y_fin_new, x_ini_new:x_fin_new ] = imagem[ y_ini_old:y_fin_old, x_ini_old:x_fin_old ];

	# If 'mask', maintain just "central" object on it..
	#
	if ( mask ):
		msk = ( mask[0]-y_ini, mask[1]-x_ini )

		zip_m = zip( msk[0], msk[1] );
		zip_z = zip( ind_z[0], ind_z[1] );

		L = list(set(zip_z) - set(zip_m));

		try:
			ind_0, ind_1 = zip(*L);
			indx = ( np.array(ind_0), np.array(ind_1) );
			imagemnova[ indx ] = 0;
		except:
			pass;

	return (imagemnova, hdr);

# ---

#==========================================================================================
def sextamp(seg_img, obj_img, header=None, increase=0, relative_increase=False, objIDs=[]):
    """
    Identify objects on given images by their IDs and return object images

    sextamp( seg_img.ndarray, obj_img.ndarray [,...] )

    By default, if 'objIDs' is not given, postamp will scan segmentation image ('seg_img') for
    the list of object ID numbers. If 'objIDs' is given, those IDs will be used for object
    poststamps creation.

    'increase' and 'relative_increase' define whether the poststamps will have a size different
    from object-only dimensions or just the object (pixels) itself, and, also, if this increasement
    value is a multiplicative factor (True) or an additive one (False).

    Since a list with object IDs is given, a list with arrays, with each IDentified
    object, is returned.

    Input:
     - seg_img           : image with segmented objects (e.g, SEx's segmentation image)
     - obj_img           : image with objects (observed pixel values)
     - header            : FITS header to be updated and passed for each poststamp
     - increase          : float value for poststamp resizing (>0)
     - relative_increase : Is 'increase' a additive value (default) or multiplicative one(?)
     - objIDs            : List with objects ID in 'seg_img'. Default is to scan 'seg_img' for IDs

    Output:
     - (ndarrays,headers)  : lists with image (poststamps) arrays and corresponding headers

    """


    # Initialize outputs:
    objs_list = [];
    hdrs_list = [];

    # If object IDs are not given, scan segmentation image for those ids..
    #
    if ( objIDs == [] ):
	    objIDs = list( set( seg_img.flatten() ) - set([0]) );

    # For each object (id) scan the respective indices for image mask and 'cutout' params
    #
    for id in objIDs:

        ind = np.where(seg_img == id);

	y_min = min( ind[0] );
	x_min = min( ind[1] );

	y_idx = ind[0] - y_min;
	x_idx = ind[1] - x_min;

	y_size = max( y_idx ) + 1;
	x_size = max( x_idx ) + 1;

	# Central pixel on original image:
	yo = y_size/2 + y_min;
	xo = x_size/2 + x_min;

	if ( increase != 0 ):

		if (relative_increase == True):

			x_size = x_size*increase;
			y_size = y_size*increase;

		else:

			x_size = x_size + 2*increase;
			y_size = y_size + 2*increase;


	try:
		hdr = header.copy();
	except:
		hdr = None;


	# Get the final image from 'cutout' output..
	#
        image_out, hdr = cutout( obj_img, header=hdr, xo=int(xo), yo=int(yo), x_size=int(x_size), y_size=int(y_size), mask=ind );

        objs_list.append( image_out );
	hdrs_list.append( hdr );


    return ( objs_list,hdrs_list );

# ---


# \cond
#==========================
if __name__ == "__main__" :

	from optparse import OptionParser;

	usage="Usage: %prog [options] <image.fits>"
	parser = OptionParser(usage=usage);

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
#	parser.add_option('--increase',
#			  dest='incr', default=0,
#			  help='Amount of border/enlarge of image. It is an additive or multiplicative factor. See "relative_increase" flag');
#	parser.add_option('--relative_increase', action='store_true',
#			  dest='rel_incr', default=False,
#			  help='Turn on the multiplicative factor for increase image.');
	parser.add_option('--header', action='store_true',
			  dest='use_header', default=False,
			  help="Use image header for WCS adjustment(?)");


	(opts,args) = parser.parse_args();

	outfits = opts.output_file;
	coord = opts.coord_unit;
	x = opts.xo;
	y = opts.yo;
	size = opts.size_unit;
	dx = opts.x_size;
	dy = opts.y_size;
#	incr = opts.incr;
#	rel_incr = opts.rel_incr;
	use_header = opts.use_header;

	if ( len(args) != 1 ):
		parser.error("Wrong number of arguments. Try option '--help' for usage and help messages.")

	infits = args[0];

	if (use_header):
		image, hdr = pyfits.getdata(infits,header=True);
	else:
		image = pyfits.getdata(infits, header=False);
		hdr = None;

#	if (incr == 0):
#		imgcut, hdr = cutout( image, header=hdr, coord_unit=coord, xo=x, yo=y, size_unit=size, x_size=dx, y_size=dy );
#	else:
#		imgcut, hdr = poststamp( image, header=hdr, increase=incr, relative_increase=rel_incr )

	imgcut, hdr = cutout( image, header=hdr, coord_unit=coord, xo=x, yo=y, size_unit=size, x_size=dx, y_size=dy );

	os.system('[ -f %s ] && rm %s' % (outfits,outfits));
	pyfits.writeto( outfits, imgcut, hdr );


	sys.exit(0);

# ------------------
# \endcond
