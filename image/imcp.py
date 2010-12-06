#!/usr/bin/env python

"""Package to deal with image copy/cuts"""

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

import sltools;
from sltools.string import regexp;
from sltools.image import segobjs;

# \cond
#----------------------
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
		CUNIT1 = 'deg';
		CUNIT2 = 'deg';

	try:
		CD1_1 = float(hdr['CD1_1']);
		CD1_2 = float(hdr['CD1_2']);
		CD2_1 = float(hdr['CD2_1']);
		CD2_2 = float(hdr['CD2_2']);
	except:
		CD1_1 = CD1_2 = CD2_1 = CD2_2 = 1;

	# Convertion from degrees to arcsec units..
	#
	if ( CUNIT1 == 'deg' and CUNIT2 == 'deg' ):
		conv_factor = 3600.0;

    # Scale is given, using CD1_1, CD1_2, CD2_1, CD2_2 header keys:
    #scale = 3600. * sqrt ((cd11**2+cd21**2+cd12**2+cd22**2)/2.) 

	x_arcsec_to_pixel = abs(CD1_1) * conv_factor;
	y_arcsec_to_pixel = abs(CD2_2) * conv_factor;


	return (y_arcsec_to_pixel, x_arcsec_to_pixel);



def _hdr_update(hdr, x_ini, y_ini):
	"""Update header information about world coordinate system"""


	LTV1 = None;
	LTV2 = None;

	NAXIS1 = int(hdr['NAXIS1']);
	NAXIS2 = int(hdr['NAXIS2']);

	try:
		CRPIX1 = float(hdr['CRPIX1']);
		CRPIX2 = float(hdr['CRPIX2']);

		CD1_1 = float(hdr['CD1_1']);
		CD1_2 = float(hdr['CD1_2']);
		CD2_1 = float(hdr['CD2_1']);
		CD2_2 = float(hdr['CD2_2']);
	except:
		CRPIX1 = CRPIX2 = 1.0;
		CD1_1 = CD1_2 = CD2_1 = CD2_2 = 1.0;
		
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

#----------------------
# \endcond

#=============================================================================
def copy_objimgs(seg_img, obj_img, objIDs):
    """
    Identify objects on given images by their IDs and return object images

    objects_IDfy( objIDs.list, seg_img.ndarray, obj_img.ndarray )

    Function identifies the pixels of each ID(int) in 'objIDs' on segmentation
    image(array) - seg_img. These identified pixels (their values) are copied 
    from object image, obj_img, to a new blank array with same shape as given images.
    (BTW, it is clear that seg_img.shape == obj_img.shape, right...)

    Since a list with object IDs is given, a list with arrays, with each IDentified
    object, is returned.

    Input:
     - seg_img : image with segmented objects, with pixels values in objIDs
     - obj_img : image with objects (observed pixel values)
     - objIDs  : a list with (int) numbers that identify the objects in seg_img

    Output:
     - objs_list : a list with arrays containing each identified object

    """

    objs_list = [];

    if (objIDs == []):
        objIDs = list(set(seg_img.flatten()) - set([0]))
        
    # In case the user just give a number (one object ID..), we deal with that
    # Lets try first the list, if not, use the single ID and do the job..
    #
    try:
        for objID in objIDs:
            merged_img = np.zeros( seg_img.shape, dtype = obj_img.dtype );
            _ind = np.where(seg_img == objID);
            merged_img[_ind] = obj_img[_ind];

            objs_list.append( merged_img );

    except:
        objID = objIDs;
        # Identify and copy the pixels in obj_img to a new blank image
        #
        merged_img = np.zeros( seg_img.shape, dtype = obj_img.dtype );
        _ind = np.where(seg_img == objID);
        merged_img[_ind] = obj_img[_ind];

        objs_list.append( merged_img );


    return (objs_list);

# ---

# =================================================================================================================
def cutout( image, header=None, coord_unit='pixel', xo=0, yo=0, size_unit='pixel', x_size=0, y_size=0, mask=None ):
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

	If 'xo=0' and 'yo=0', given image (input_file) central pixel will be chosen as (xo,yo).
	If 'x_size=0' and 'y_size=0', half length of each side will be used for output dimensions.

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


	# Get the side sides (at least, 1!) and transform for the size_unit if necessary..
	#
	x_cut_size = max( 1, int(float(x_size)) );
	y_cut_size = max( 1, int(float(y_size)) );

	if ( size_unit == 'degrees' ):
		x_cut_size = int(float(x_cut_size)/x_arcsec_to_pixel);
		y_cut_size = int(float(y_cut_size)/y_arcsec_to_pixel);


	# And if no side size was given, define a default value correspondig to half of original image..
	#
	if not ( x_size ):
		x_cut_size = int(x_img_size/2);
		print >> sys.stdout, "Warning: 'x_size' not given. Using half of image x side(%d)" % (x_cut_size);

	if not ( y_size ):
		y_cut_size = int(y_img_size/2);
		print >> sys.stdout, "Warning: 'y_size' not given. Using half of image y side(%d)" % (x_cut_size);


	# Check requested output.vs.input image sizes..
	#
	if ( x_cut_size == x_img_size and y_cut_size == y_img_size ):
		print >> sys.stdout, "Warning: Requested output sizes are the same as input image. Returning image and header unchanged.";
		return (imagem,hdr);


	# Verify central coordinates values..
	#
	if ( coord_unit == 'pixel' and (xo != 0 and yo != 0) ):
		x_halo = int(xo);
		y_halo = int(yo);

	elif ( coord_unit == 'degrees' and (xo != 0 and yo != 0) ):
		_out = commands.getoutput( 'sky2xy %s %s %s' % (input_file, xo, yo) );
		_out = string.split(_out);
		x_halo = int(_out[-2])-1;
		y_halo = int(_out[-1])-1;

	elif ( xo == 0 and yo == 0 ):
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
		hdr.update('NAXIS1',x_cut_size);
		hdr.update('NAXIS2',y_cut_size);

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

# ==========================================================================================
def segstamp(obj_img, seg_img, header=None, increase=0, relative_increase=False, objIDs=[]):
    """
    Identify objects on given images by their IDs and return object images

    segstamp( seg_img.ndarray, obj_img.ndarray [,...] )

    By default, if 'objIDs' is not given, postamp will scan segmentation image ('seg_img') for
    the list of object ID numbers. If 'objIDs' is given, those IDs will be used for object
    poststamps creation.

    'increase' and 'relative_increase' define whether the poststamps will have a size different
    from object-only dimensions or just the object (pixels) itself, and, also, if this increasement
    value is a multiplicative factor (True) or an additive one (False).

    Since a list with object IDs is given, a list with arrays, with each IDentified
    object, is returned.

    Input:
     - obj_img           : image with objects (observed pixel values)
     - seg_img           : image with segmented objects (e.g, SEx's segmentation image)
     - header            : FITS header to be updated and passed for each poststamp
     - increase          : float value for poststamp resizing (>0)
     - relative_increase : Is 'increase' a additive value (default) or multiplicative one(?)
     - objIDs            : List with objects ID(integers) in 'seg_img'. Default is to scan 'seg_img' for IDs

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

    if type(objIDs) is int :
        objIDs = [objIDs];
    
    for _id in objIDs:
    
        if (_id == 0):
            objs_list.append(None);
            hdrs_list.append(None);
            continue;
            
        ind = segobjs.create_IDmask(seg_img, _id)[0];

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
# To keep compatibility with the old "sextamp" function:
def sextamp(seg_img, obj_img, header=None, increase=0, relative_increase=False, objIDs=[]):
	return segstamp(obj_img, seg_img, header, increase, relative_increase, objIDs)
# ---

# \cond
# ============================================================================
if __name__ == "__main__" :

	from optparse import OptionParser;

	usage="Usage: %prog [options] <image.fits>"
	parser = OptionParser(usage=usage);

	parser.add_option('-o',
			  dest='output_file', default='pstamp_',
			  help='Output (image) file rootname ["pstamp_"]');
	parser.add_option('--coord',
			  dest='coord_unit', default='pixel',
			  help='Coordinate units, for xo & yo [pixel]');
	parser.add_option('--xo',
			  dest='xo', default=0,
			  help='X centre');
	parser.add_option('--yo',
			  dest='yo', default=0,
			  help='Y centre');
	parser.add_option('--size_unit',
			  dest='size_unit', default='pixel',
			  help='Side length unit, for x_size & y_size [pixel]');
	parser.add_option('--x_size',
			  dest='x_size', default=0,
			  help='Size of x side on output');
	parser.add_option('--y_size',
			  dest='y_size', default=0,
			  help='Size of y side on output');
	parser.add_option('--no-header', action='store_false',
			  dest='use_header', default=True,
			  help="Do not use existing image header.");
	parser.add_option('--outdir',
			  dest='dir_name', default='out_imcp',
			  help="Directory name to store output images(poststamps)");
	parser.add_option('-s', '--segimg',
			  dest='segimg', default=None,
			  help="Segmentation image (e.g, Sextractor's SEGMENTATION output)");
	parser.add_option('--objIDs',
			  dest='IDs',default='',
			  help="Comma-separated list of object IDs found on 'segimg' to extract");
	parser.add_option('--increase',
			  dest='incr', default=0,
			  help="Amount of border/enlarge of object poststamp. It is an additive or multiplicative factor. See 'relative_increase' flag");
	parser.add_option('--relative_increase', action='store_true',
			  dest='rel_incr', default=False,
			  help="Turn on the multiplicative factor for increase image, otherwise, 'increase' will be treated as an absolute value.");

	(opts,args) = parser.parse_args();

	outfits = opts.output_file;
	coord = opts.coord_unit;
	x = opts.xo;
	y = opts.yo;
	size = opts.size_unit;
	dx = opts.x_size;
	dy = opts.y_size;
	segimg = opts.segimg;
	IDs = opts.IDs;
	incr = opts.incr;
	rel_incr = opts.rel_incr;
	use_header = opts.use_header;
	dir_name = opts.dir_name;

	if ( len(args) != 1 ):
		parser.error("Wrong number of arguments. Try option '--help' for usage and help messages.")

	infits = args[0];

	if (use_header):
		image, hdr = pyfits.getdata(infits,header=True);
	else:
		image = pyfits.getdata(infits, header=False);
		hdr = None;

	try:
		os.mkdir( dir_name );
	except OSError:
		pass;

	owd = os.getcwd()+"/";
	os.chdir( dir_name );

	if (segimg):
		IDs = regexp.str2lst(IDs);
		imglst, hdrlst = segstamp(image, segimg, header=hdr, increase=incr, relative_increase=rel_incr, objIDs=objIDs)
	else:
		imgcut, hdr = cutout( image, header=hdr, coord_unit=coord, xo=x, yo=y, size_unit=size, x_size=dx, y_size=dy );
		IDs = ['0'];
		imglst = [imgcut];
		hdrlst = [hdr];

	for i in range(len(IDs)):
		outfits = outfits+"%s.fits" % (IDs[i]);
		os.system('[ -f %s ] && rm -f %s' % (outfits,outfits));
		pyfits.writeto( outfits, imglst[i], hdrlst[i] );

	print >> sys.stdout, "Done.";
	os.chdir( owd );
	sys.exit(0);

# ------------------
# \endcond
