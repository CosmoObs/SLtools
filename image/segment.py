#!/usr/bin/env python

"""Module to deal with (Sextractor) segmentation procedures"""

##@ segment
#
#
# This package contains functions to deal with image objects segmentation.
# Functions are designed primarily to work with Sextractor
#
# Executable package : No


import sys;

import os;
import re;
import string;
import pyfits;
import numpy as np;

import sextractor as se;


#=============================================================================
# Compatibility for old 'segment' module:
run_sex = se.run;

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

#==========================================================================================================================
def images_IDfy(catalog, seg_image, obj_image, out_rootname="IDfy_out_"):
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
    cat = np.loadtxt(catalog, comments="'");


    # Object IDs
    objIDs = list(cat[:,0]);   # Notice that I'm supposing that the object IDs
                               # are placed in the first column of given catalog

    objIDs = [ int(id) for id in objIDs ];
    imgs_list = objects_IDfy(objIDs, seg, obj);

    # Write down each returned fits image..
    #
    for i in range(len(imgs_list)):

        _id = objIDs[i];
        pyfits.writeto( out_rootname+str(_id)+'.fits', imgs_list[i], hdr );


    return (True);

