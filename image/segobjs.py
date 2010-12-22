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
     - segimg : ndarray(ndim=2,dtype=int)
       ndarray with a segmented image (int), e.g, SE's SEGMENTATION 
     - centroids : [(int,int),]
       List of tuples with positions to search for objects [(x0,y0),(x1,y1),...]

    Output:
     - object IDs : [int,]
       List of IDs (int) for each given (x,y) points
       
       
    Example:
    >>> x = [1267,768,963]
    >>> y = [851,564,733]
    >>> centroids = zip(x,y)
    >>> objIDs = centroid2ID(SegImg_array,centroids)
    >>> objIDs
    [2,4,3]
    
    """

    objIDs = [];
    
    # Is there an identified object on given (xo,yo) point? If yes, store its ID..
    #
    try:
        for xo,yo in centroids:
            objid = segimg[int(float(yo)),int(float(xo))];
            objIDs.append(objid);
            if not (objid):
                print >> sys.stdout, "No objects were identified for (x=%s,y=%s) position." % (xo,yo);
                print >> sys.stdout, "\"ID\" 0 means \"no object identified\".";
    except TypeError:
        print >> sys.stderr, "Error: Points were not given as expected. centroids=[(x1,y1), ...]"
        del objIDs;
        return (False);
        
    return (objIDs);


#=============================================================================
def create_IDmask(segimg, objIDs=[]):
    """ Generate mask for each object in given image.
    
    ID in 'objIDs' is used to create a mask for each object (ID) in 'segimg'.
    
    Input:
     - segimg : ndarray(ndim=2,dtype=int)
        Image array with int numbers as object identifiers
     - objIDs : [int,]
        List with IDs for objects inside 'segimg'
        
    Output:
     -> [(ind_array,ind_array),]
        List of tuples with index arrays in it. The output is a list of "numpy.where()" arrays.
        
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
    
