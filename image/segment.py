#!/usr/bin/env python

"""Package to deal with (Sextractor) segmentation procedures"""

##@ segment
#
#
# This package contains functions to deal with image objects segmentation.
# Functions are designed primarily to work with Sextractor
#
# Executable package : No
#


import sys;

import os;
import re;
import string;
import pyfits;
import numpy as np;


#=======================================================================
def run_sex(fits_image, params=["NUMBER","X_IMAGE","Y_IMAGE"], args={}):
    """
    Run Sextractor on given (fits) image, with optional parameters

    run_sex( 'image.fits' [ [...],{...} ] )

    Function runs Sextractor on a given FITS image (fits_image) to generate 
    sex's "OBJECTS" and "SEGMENTATION" outputs. Output parameters (catalogs) 
    can be passed in a list of strings. SE command-line arguments can
    be passed as a dictionary {'key' = 'value'}.
    * Output parameter "NUMBER" is always placed as the first column on catalog

    Two catalogue files are created during sex's runs: 'catalog_obj.cat' and 
    'catalog_seg.cat', with "objects" and "segmentation" features output 
    values, respectively.
    Two FITS files - images - are also created during the runs, one for each
    run type: "objects" and "segmentation". File names are just the given 
    filename with extension ".fits" changed by "_obj.fits" and "_seg.fits", 
    respectively.

    Input:
     - fits_image : FITS image filename
     - params     : List of parameters to output on catalogues
     - args       : Dictionary with sextractor line arguments. If not given,
                    'sex -dd' will be used as default options.

    Output:
     * Two ascii files and two fits images. See above paragraph.

    """

    # Sextractor config file..
    #
#    os.system( 'sex -dd > default.sex' );

    # Object features to output..
    #
    fparam = open('default.param','w');
    fparam.write( "NUMBER\n" );               # I've fixed the first param for compatibility reasons.
    for _param in params:
            fparam.write( "%s\n" % (_param) );
    fparam.close();


    # Build sextractor command line..
    #
    args = {};
    cmd_line = ' ';
    for key in args:
        cmd_line = cmd_line + ' -' + key + ' ' + kwargs[key];


    rootname = string.replace( fits_image,".fits","" );

    # OBJECTS output
    outobjimage = rootname + "_obj.fits";
    status= os.system( 'sex %s -CHECKIMAGE_TYPE OBJECTS -CATALOG_TYPE ASCII -CATALOG_NAME catalog_obj.cat -CHECKIMAGE_NAME %s %s' % (fits_image,outobjimage,cmd_line));
    if ( status != 0 ): print >> sys.stderr, "Warning: Sextractor raised an error code '%s' during 'OBJECTS' run. Continuing." % (status);

    # SEGMENTATION output
    outsegimage = rootname + "_seg.fits";
    status= os.system( 'sex %s -CHECKIMAGE_TYPE SEGMENTATION -CATALOG_TYPE ASCII -CATALOG_NAME catalog_seg.cat -CHECKIMAGE_NAME %s %s' % (fits_image,outsegimage,cmd_line));
    if ( status != 0 ): print >> sys.stderr, "Warning: Sextractor raised an error code '%s' during 'SEGMENTATION' run. Continuing." % (status);


    if (status):
        return (False);

    return (True);

# ---

#=============================================================================
def objects_IDfy(objIDs, seg_img, obj_img):
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
     - objIDs  : a list with (int) numbers that identify the objects in seg_img
     - seg_img : image with segmented objects, with pixels values in objIDs
     - obj_img : image with objects (observed pixel values)

    Output:
     - objs_list : a list with arrays containing each identified object

    """

    objs_list = [];

    # In case the user just give a number (one object ID..), we deal with that
    # Lets try first the list, if not, use the single ID and do the job..
    #
    try:
        for objID in objIDs:
            merged_img = np.zeros( seg_img.shape, dtype=int );
            _ind = np.where(seg_img == objID);
            merged_img[_ind] = obj_img[_ind];

            objs_list.append( merged_img );

    except:
        objID = objIDs;
        # Identify and copy the pixels in obj_img to a new blank image
        #
        merged_img = np.zeros( seg_img.shape, dtype=int );
        _ind = np.where(seg_img == objID);
        merged_img[_ind] = obj_img[_ind];

        objs_list.append( merged_img );


    return (objs_list);

# ---

#==========================================================================================================================
def image_IDfy(catalog, seg_image, obj_image, out_rootname="IDfy_out_"):
    """
    Identify objects on given images by their IDs on catalog and write object images

    image_IDfy(catalog.txt, seg_image.fits, obj_image.fits, out_rootname.string, relative_size.bool, force_relative.bool)

    Identifies the pixels of each ID(int) in catalog (1st column) on segmentation
    image - seg_image. These identified pixels (their intensities) are copied 
    from object image, obj_image, to a new blank image with same shape as given ones.
    (BTW, it is clear that seg_image and obj_image have the same dimensions)

    For each object in catalog an image with corresponding object pixels will be created.
    The name of new (output) images will be composed by "out_rootname"+"ID"+".fits".

    Input:
     - catalog        : ASCII catalog. The first column should the ID of objects in seg_image
     - seg_image      : FITS image with segmented objects, with pixels values in catalog
     - obj_image      : FITS image with objects (observed pixel values)
     - out_rootname   : Basename for output images naming, plus the object ID (FITS file)

    Output:
    * Depending on the size of catalog, fits images will be created

    """


    # Read input images and the ascii catalog (header read from obj image)..
    #
    obj, hdr = pyfits.getdata(obj_image, header=True);
    seg = pyfits.getdata(seg_image);
    cat = np.loadtxt(catalog, comments="#");


    # Object IDs
    objIDs = list(cat[:,0]);   # Notice that I'm supposing that the object IDs
                               # are placed in the first column of given catalog

    imgs_list = objects_IDfy(objIDs, seg, obj);

    # Write down each returned fits image..
    #
    for img in imgs_list:

        pyfits.writeto( out_rootname+'%s.fits', img, hdr % (ID));


    return (True);

# ---
# I'll leave it here as a good record of basic necessary reg-exp use.
# It was used before the (only!) numpy line above.
#
#    fp_cat = open(catalog,'r');
#    objIDs = [];
#    for line in fp_cat.readlines():
#        line.rstrip("\n");
#        if ( re.search("^$",line) ): continue;
#        if ( re.search("^\s*#",line) ): continue;
#        line = re.sub("^\s*","",line);
#        line = string.split(line);
#        objIDs.append( int(line[0]) );
#    fp_cat.close();
# ---

# @cond
#==========================
if __name__ == "__main__" :

	from optparse import OptionParser;

	usage="Usage: %prog [options] <image_fits_file>"

	parser = OptionParser(usage=usage);

	parser.add_option('--run_sex', action='store_true'
			  dest='run_sex', default=True,
			  help='Run sextractor on given FITS image');
        parser.add_option('--params_file',
                          dest='params', default='default.param',
                          help='Sextractor params file [default.param]');
        parser.add_option('--config_file',
                          dest='config', default='default.sex',
                          help='Sextractor config file [default.sex]');

	(opts,args) = parser.parse_args();

        run_sex = opts.run_sex;
        params_file = opts.params;
        config_file = opts.config;

	if ( len(args) != 1 ):
		parser.error("Wrong number of arguments. Try option '--help' for usage and help messages.")

	infits = args[0];

	if ( run_sex ):

            fp_params = open(params_file);
            params = [ line.rstrip("\n")   for line in fp_params.readlines() ]
            config = config_parser(config_file);


	sys.exit(0);

# ------------------
# @endcond
