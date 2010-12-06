#!/usr/bin/env python

"""Module to deal with objects identification in segmented images"""

##@ segobjs
#
#
# This package contains functions to deal with image objects segmentation.
# Functions are designed primarily to work with Sextractor
#
# Executable package : No


import sys;

import pyfits;
import numpy as np;


# =======================================================================
def centroid2ID(segimg, centroids=[]):
    """Select objects based on given positions

    Input:
     - segimg : ndarray
       ndarray with a segmented image (int), e.g, SE's SEGMENTATION 
     - centroids : [(int,int),] :
       List of tuples with positions to search for objects [(x0,y0),(x1,y1),...]

    Output:
     - object IDs : [int,]
       List of IDs (int) for each given (x,y) points
       
    """

    objIDs = [];
    
    # Is there an identified object on given (xo,yo) point? If yes, store its ID..
    #
    try:
        for xo,yo in centroids:
            objid = segimg[int(yo),int(xo)];
            objIDs.append(objid);
            if not (objid):
                print >> sys.stdout, "No objects were identified for (x=%s,y=%s) position." % (xo,yo);
                print >> sys.stdout, "Be careful, \"ID\" 0 was appended to the output to fill the position.";
    except TypeError:
        print >> sys.stderr, "Error: Points were not given as expected. Try something like: centroids=[(x1,y1), (,), ...]"
        del objIDs;
        return (False);
        
    return (objIDs);


#=============================================================================
def create_IDmask(segimg, objIDs=[]):
    """ Generate mask for each object in given image.
    
    Each ID in 'objIDs' is used to create a mask for found objects in 'segimg'.
    
    Input:
     - segimg : ndarray(ndim=2,dtype=int)
        Image array with int numbers as object identifiers
     - objIDs : [int,]
        List with IDs for objects inside 'segimg'
        
    Output:
     -> [(ndarray,ndarray),]
        List of tuples with index arrays in it. The output is a list of "numpy.where()" outs.
        
    """

    indx = [];
    
    # For each object (id) scan the respective indices for image mask and 'cutout' params
    #
    try:
        for id in objIDs:
            # 'ind' will be used (also) as the 'mask' information on "cutout" function
            ind = np.where(segimg == int(id));
            indx.append(ind);
    except:
        id = objIDs;
        ind = np.where(segimg == int(id));
        indx.append(ind);

    return indx;
    
"""
#=============================================================================

    FUNCTION MOVED TO imcp.py/copy_objimgs().

def copy_IDfied(seg_img, obj_img, objIDs):
    \"""
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

    \"""

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

#===============================================================================================
if __name__ == "__main__" :

    import sys;
    from optparse import OptionParser;

    usage_msg = "%prog  [options]  <segmentation_image.fits> <objects_image.fits>";
    epilog_msg = "
    Identify objects on given images by their IDs or centroid positions.
    
    IDs are meant to be the objects IDentification used inside \"segmentation\" image('seg').
    IDs are expected to be found in the first column of the ASCII table passed through
    'cat' argument. If no catalog is given, a list of all IDs (non-zero values) will be
    created from 'seg' image.
    \"--XY\" option is used, in place of IDs, catalog's first and second columns should contain
    Xo Yo point positions.
    
    * The catalog is expected to be a data block table, column separator is any whitespace,
    centroids (x,y) should be placed as such: 'X' 'Y'.
    I.e.,
      if --XY:  X0 Y0 ...
                X1 Y1 ...
                ...
            or
                X0 Y0 ...
                X1 Y1 ...
                ...
    
    
    Identified object pixels (their intensities) are copied from \"object\" image('obj')
    to a new blank image with same shape as given ones.
    (BTW, it should the clear that 'seg' and 'obj' have the same dimensions)
    
    For each object in catalog an image with corresponding object pixels will be created.
    The name of new (output) images will be composed by \"outname + ID + .fits\".
    "
    
#    oParser = OptionParser(usage=usage_msg, epilog=epilog_msg);
    oParser = OptionParser(usage=usage_msg);

    oParser.add_option('-c', '--cat',
                        dest='table', default=None,
                        help="ASCII catalog with object IDs. The IDs should be in catalog's first column.")
    oParser.add_option('--XY', action='store_true',
                        dest='is_xy', default=False,
                        help="Use image positions (x,y) instead of IDs in given catalog/table (@ 1st column)");
    oParser.add_option('-o','--outname',
                        dest='rootname', default='outobj_',
                        help="Output file rootname to use: 'outname + ID + .fits'");
						
    (opts,args) = oParser.parse_args();
    if not (len(args) == 2):
        oParser.print_usage();
        sys.exit(1);
        
    cat = opts.table;
    out_rootname = opts.rootname;
    use_centroid = opts.is_xy;
    
    # Read input images and the ascii catalog (header read from obj image)..
    #
    obj, hdr = pyfits.getdata(args[1], header=True);
    seg = pyfits.getdata(args[0]);

    try:
        cat = np.loadtxt(cat, comments="#", dtype=int);
        if (use_centroid):
            objIDs = centroids2IDs(seg, centroids=zip(cat[:,0],cat[:,1]));
        else:
            objIDs = list(cat[:,0]);
    except:
        objIDs = [];

    imgs_list = objects_IDfy(seg, obj, objIDs);

    # Write down each returned fits image..
    #
    for i in range(len(imgs_list)):

        _id = objIDs[i];
        pyfits.writeto( out_rootname+str(_id)+'.fits', imgs_list[i], hdr );


    sys.exit(0);
"""
