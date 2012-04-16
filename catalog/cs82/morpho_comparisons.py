#!/usr/bin/env python
# ====================================================
# Authors:
# Bruno Moraes - bruno.a.l.moraes@gmail.com - 01/Nov/2011
# ====================================================


"""Module to ..."""

##@package morpho_comparisons
#
# This package does
# 
# The code 
#
# Executable package: YES
#
# To see the package help message, type:
#
# > python morpho_comparisons.py --help
#
# To run the code:
#
# > python morpho_comparisons.py 'bla' 'bla'



from __future__ import division
import sys
import os
import errno
import numpy as np
import pyfits
import matplotlib.pyplot as pl
from time import strftime

from sltools.catalog import data_subsets as dts
from sltools.catalog import fits_data as fd
from sltools.catalog import ascii_data as ascd
from sltools.catalog import table_matching as tbm
from sltools.coordinate import wcs_conversion as wcscnv
from sltools.io import io_addons as ioadd
from sltools.plot import plot_templates as pltemp




#################################################
####           MATCHING FUNCTIONS            ####
#################################################


def create_matching_arrays(cat_A, cat_B, coord_A_names, coord_B_names, radius):

    """
    Creates a distance array and a matching array for cat_A and cat_B objects.

    This function uses a few Complex matrix operations to build a (dim A) x
    (dim B) Real matrix where the entry (i,j) gives the squared distance 
    between the ith object in cat_A and the jth object in cat_B. From this 
    matrix, a simple numpy logical test returns a matrix of same dimensions, 
    with True entries for dist_matrix[i][j] < radius^2 and False entries
    otherwise.


    Input:
     - cat_A             hdudata : FITS Table data catalog
     - cat_B             hdudata : FITS Table data catalog
     - coord_A_names   [str,str] : ra_A and dec_A column names
     - coord_B_names   [str,str] : ra_B and dec_B column names
     - radius              float : Radius to perform matching
    
    Output:
     - dist_matrix   numpy.ndarray : [i,j] distance array (ndim=2,dtype=float) 
     - match_matrix  numpy.ndarray : [i,j] matching array (ndim=2,dtype=bool)
    
    ---
    """

    # Creating complex coordinate vectors

    cat_A_ra = cat_A.field(coord_A_names[0])
    cat_A_dec = cat_A.field(coord_A_names[1])
    cat_A_radec = cat_A_ra + 1j*cat_A_dec

    cat_B_ra = cat_B.field(coord_B_names[0]).real
    cat_B_dec = cat_B.field(coord_B_names[1]).real
    cat_B_radec = cat_B_ra + 1j*cat_B_dec
    
    
    # Matrix operation: C = (pA).1_(1xn) - 1_(mx1).(pB)    

    coord_matrix = np.outer(cat_A_radec,np.ones(len(cat_B_radec)))-np.outer(np.ones(len(cat_A_radec)),cat_B_radec)


    # Element-wise complex conjugate multiplication and square-root

    dist_matrix = np.real(coord_matrix*np.conjugate(coord_matrix))
     

    # Logical matching matrix

    match_matrix = (dist_matrix < radius*radius)

    return dist_matrix, match_matrix



def exclude_unmatched_objects(cat_A, cat_B, dist_matrix, match_matrix, output_unmatched = False):

    """
    Finds unmatched objects in catalogs A and B, excludes them from catalogs.

    This function works with the outputs of the function create_matching_arrays
    (cf. above). It finds lines and columns in the matching matrix that contain
    only False (i.e., lines/columns where no matching was found) and excludes
    them from the catalogs and from the distance matrix and matching matrix.
    As an option, it can output catalogs of the unmatched objects.


    Input:
     - cat_A                hdudata : FITS Table data catalog
     - cat_B                hdudata : FITS Table data catalog
     - dist_matrix    numpy.ndarray : [i,j] distance array (ndim=2,dtype=float)
     - match_matrix   numpy.ndarray : [i,j] matching array (ndim=2,dtype=bool)
     - output_unmatched        bool : Outputs unmatched objects if True
    
    Output:
     - cat_A_matched        hdudata : FITS Table data catalog of matched objs
     - cat_B_matched        hdudata : FITS Table data catalog of matched objs
     - dist_matrix    numpy.ndarray : [i,j] distance array (ndim=2,dtype=float)
     - match_matrix   numpy.ndarray : [i,j] matching array (ndim=2,dtype=bool)

    Optional Outputs:
     - cat_A_unmatched    hdudata : FITS Table data catalog of unmatched objs
     - cat_B_unmatched    hdudata : FITS Table data catalog of unmatched objs    
    
    ---
    """


    matched_A = np.any(match_matrix, axis=1)

    print "There are " + str((~matched_A).sum()) + " unmatched objects in catalog A"

    matched_B = np.any(match_matrix, axis=0)

    print "There are " + str((~matched_B).sum()) + " unmatched objects in catalog B"

    cat_A_matched = cat_A[matched_A]
    cat_B_matched = cat_B[matched_B]

    dist_matrix = dist_matrix[matched_A]
    dist_matrix = dist_matrix[:,matched_B]

    match_matrix = match_matrix[matched_A]
    match_matrix = match_matrix[:,matched_B]

    if output_unmatched:

        cat_A_unmatched = cat_A[~matched_A]

        cat_B_unmatched = cat_B[~matched_B]

        return cat_A_matched, cat_B_matched, dist_matrix, match_matrix, cat_A_unmatched, cat_B_unmatched

    else:

        return cat_A_matched, cat_B_matched, dist_matrix, match_matrix


def match_neighbors(cat_A_matched, cat_B_matched, dist_matrix, match_matrix, exclude_multiple=True, output_multiple=False):

    """
    Gets reordered catalog B with nearest neighbors to each catalog A entry

    This function uses the outputs of the function exclude_unmatched_objects
    (cf. above). It takes the two matched catalogs and the matching matrices, 
    and finds the nearest neighbor in catalog B to each of catalog A objects.
    It then returns cat_A and cat_B reordered so that the nearest neighbors
    are in the same line.

    There are two possible options concerning cases where there's multiple
    matching. The first (default is True) allows the user to find the nearest
    neighbor only in the cases where the matching is unique in both senses 
    (i.e. a bijection). Without this option, the nearest neighbor will still
    be found, but the secondary matchings will be ignored. Some of the 
    comparisons resulting may be physically unsound. The second option 
    (default is False), *not yet implemented*, creates .reg masks of the 
    multiple matching cases, allowing the user to cross-check these cases
    visually.


    Input:
     - cat_A_matched        hdudata : FITS Table data catalog of matched objs
     - cat_B_matched        hdudata : FITS Table data catalog of matched objs
     - dist_matrix    numpy.ndarray : [i,j] distance array (ndim=2,dtype=float)
     - match_matrix   numpy.ndarray : [i,j] matching array (ndim=2,dtype=bool)
     - exclude_multiple        bool : Gets only neighbors in the case of
                                      bijective matching if True. If False,
                                      gets nearest of multiple neighbors
     - output_multiple         bool : Outputs .reg masks of multiple matchings
                                      to a file if True
    
    Output:
     - cat_A_neighbor        hdudata : FITS Table data catalog of objs
     - cat_B_neighbor        hdudata : FITS Table data catalog of nearest
                                       neighbors of objects in cat_A_neighbor
     
    ---
    """

   
    if np.any(match_matrix.sum(1)==0) or np.any(match_matrix.sum(0)==0):

        print "These catalogs still contain unmatched objects. Run the function exclude_unmatched_objects on them before using this function."
        sys.exit(1)


    filter_A = (match_matrix.sum(1) == 1)

    print "There are " + str((~filter_A).sum()) + " multiply-matched objects in catalog A"

    filter_B = (match_matrix.sum(0) == 1)
        
    print "There are " + str((~filter_B).sum()) + " multiply-matched objects in catalog B"


    if exclude_multiple:

        print "Keeping only one-to-one matching..."

        cat_A_matched = cat_A_matched[filter_A]
        cat_B_matched = cat_B_matched[filter_B]

        dist_matrix = dist_matrix[filter_A]
        dist_matrix = dist_matrix[:,filter_B]

        match_matrix = match_matrix[filter_A]
        match_matrix = match_matrix[:,filter_B]

    else:

        print "Keeping nearest neighbor in the case of multiple matching..."


    if output_multiple:

        print "Creating reg masks of multiple-matching objects..."
    
        multiple_A = ~filter_A
        multiple_B = ~filter_B

        match_matrix_multiple_A = match_matrix[multiple_A]
        match_idxs_tuples_A = match_matrix_multiple_A.nonzero()
        match_idxs = np.transpose(np.array([match_idxs_tuples_A[0],match_idxs_tuples_A[1]]))

        match_matrix_multiple_B = match_matrix[:,multiple_B]
        match_idxs_tuples_B = match_matrix_multiple_B.nonzero()
        match_idxs = np.transpose(np.array([match_idxs_tuples_B[0],match_idxs_tuples_B[1]]))
        
        # Writing this output to reg masks (including vectors to connect 
        # neighbors) will be highly non-trivial. This will have to wait.
        # For this process, the idea would be to use the indexes above 
        # to select the coordinates of each of the objects in the catalogs.
        # This will give both the center coordinates of each circle but also
        # the coordinates and direction of a .reg vector, which will allow 
        # a graphical indication of who is multiply-matched to whom.

    
    matching_mask = np.argmin(dist_matrix,axis=1)

    cat_A_neighbor = cat_A_matched
    cat_B_neighbor = cat_B_matched[matching_mask]

    return cat_A_neighbor, cat_B_neighbor



def matching_pipeline(cat_A, cat_B, coord_A_names, coord_B_names, radius, exclude_multiple=True, output_multiple=False, output_unmatched=False):

    """
    Performs nearest neighbor matching between two catalogs.

    This pipeline matches two catalogs in coordinates specified as input,
    excludes unmatched objects and finds the nearest neighbor in catalog B to
    each of catalog A objects. It then returns cat_A and cat_B reordered so 
    that the nearest neighbors are in the same line. Catalogs of the unmatched
    A and B objects may be given as optional outputs (default is False).

    There are two additional options concerning cases where there's multiple
    matching. The first (default is True) excludes multiple matching cases 
    (i.e. cases where the matching is not unique in both senses). The second
    option (default is False), *not yet implemented*, creates .reg masks of 
    the multiple matching cases, allowing the user to cross-check these cases
    visually.

    For more information on the details of the procedure, refer to the
    documentation of the functions create_matching_arrays, 
    exclude_unmatched_objects and match_neighbors.


    Input:
     - cat_A_matched        hdudata : FITS Table data catalog of matched objs
     - cat_B_matched        hdudata : FITS Table data catalog of matched objs
     - coord_A_names      [str,str] : ra_A and dec_A column names
     - coord_B_names      [str,str] : ra_B and dec_B column names
     - radius                 float : Radius to perform matching
     - output_unmatched        bool : Outputs unmatched objects if True
     - exclude_multiple        bool : Gets only neighbors in the case of
                                      bijective matching if True. If False,
                                      gets nearest of multiple neighbors
     - output_multiple         bool : Outputs .reg masks of multiple matchings
                                      to a file if True
    
    Output:
     - cat_A_matched        hdudata : FITS Table data catalog of matched objs
     - cat_B_matched        hdudata : FITS Table data catalog of matched objs

    Optional Outputs:
     - cat_A_unmatched    hdudata : FITS Table data catalog of unmatched objs
     - cat_B_unmatched    hdudata : FITS Table data catalog of unmatched objs 
    
    ---
    """


    dist_matrix, match_matrix = create_matching_arrays(cat_A, cat_B, coord_A_names, coord_B_names, radius)


    if output_unmatched:

        cat_A_matched, cat_B_matched, dist_matrix, match_matrix, cat_A_unmatched, cat_B_unmatched = exclude_unmatched_objects(cat_A, cat_B, dist_matrix, match_matrix, output_unmatched)

        cat_A_neighbor, cat_B_neighbor = match_neighbors(cat_A_matched, cat_B_matched, dist_matrix, match_matrix, exclude_multiple, output_multiple)

        return cat_A_neighbor, cat_B_neighbor, cat_A_unmatched, cat_B_unmatched


    else:

        cat_A_matched, cat_B_matched, dist_matrix, match_matrix = exclude_unmatched_objects(cat_A, cat_B, dist_matrix, match_matrix, output_unmatched)

        cat_A_neighbor, cat_B_neighbor = match_neighbors(cat_A_matched, cat_B_matched, dist_matrix, match_matrix, exclude_multiple, output_multiple)

        return cat_A_neighbor, cat_B_neighbor


