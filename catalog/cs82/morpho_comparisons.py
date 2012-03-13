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


sys.path.append("/home/brunomor/Documents/Trabalho/Repositorio") ### Append Bruno
from sltools.catalog import fits_data as fd ### Import Bruno
from sltools.catalog import ascii_data as ascd ### Import Bruno
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

    match_matrix = (dist_matrix < radius**2)

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





#################################################
####          DATA SUBSETS FUNCTIONS         ####
#################################################

def create_data_subsets(data,col_name,cuts):

    """
    Creates subsets from data, cutting in numerical value ranges of a given field

    This function takes different sets of input value ranges for a given field in
    the data and uses them to create numpy masks that are then applied to the 
    original data in order to create subsets of it. These subsets and the 
    corresponding masks are returned.

    An arbitrary number of cuts can be performed, but they must be of one of the
    following three formats:

        - ['.',b] --> performs the cut data.field(col_name) < b
        - [a,b] ----> performs the cut a <= data.field(col_name) < b
        - [a,'.'] --> performs the cut a <= data.field(col_name)

    There's currently no support for cuts in string fields or to change from
    semi-open to open or closed intervals.


    Input:
     - data                        hdr : Input data
     - col_name                    str : Name of the column chosen to perform
                                         the cuts
     - cuts  [[a_1,b_1],...,[a_n,b_n]] : Cuts for different subsets

    Output:
     - cut_data         [pyfits_bla_1,...,pyfits_bla_n] : List of data subsets
                                            (each data subset: ndim=2, dtype=)
     - masks    [ndarray_1,...,ndarray_n] : List of masks used to create the 
                                            data subsets 
                                            (each mask: ndim=1,dtype=bool) 
    """

    # Try/Except if the field doesn't exist in the catalog?    

    cut_data = []
    masks = []

    for i in range(len(cuts)):

        if cuts[i][0] == '.':

            masks.append(data.field(col_name) < cuts[i][1])

        elif cuts[i][1] == '.':
    
            masks.append(cuts[i][0] <= data.field(col_name))

        else:

            masks.append((cuts[i][0] <= data.field(col_name)) & (data.field(col_name) < cuts[i][1]))


        cut_data.append(data[masks[i]])

    return cut_data, masks


def create_matched_subsets(data_A,data_B,col_name,cuts):

    """
    Creates subsets of data A and B, cutting in numerical value ranges of a 
    given field in data A

    This function creates subsets of matched data A and B according to cuts
    performed in a given field of data A, returning two lists, containing the 
    subsets of A and B. Input data A and B must have the same number of rows,
    ideally coming from a matching procedure.

    For details of the cutting procedure, check the documentation of the
    function create_data_subsets.


    Input:
     - data_A                      hdr : Input data set A
     - data_B                      hdr : Input data set B (same size as A)
     - col_name                    str : Name of the column chosen to perform
                                         the cuts
     - cuts  [[a_1,b_1],...,[a_n,b_n]] : Cuts for different subsets

    Output:
     - cut_data_A     [pyfits_bla_1,...,pyfits_bla_n] : List of data subsets
                                            (each data subset: ndim=2, dtype=)
     - cut_data_B     [pyfits_bla_1,...,pyfits_bla_n] : List of data subsets
                                            (each data subset: ndim=2, dtype=)
    """

    cut_data_A, masks_A = create_data_subsets(data_A,col_name,cuts)

    cut_data_B = []

    for i in range(len(cut_data_A)):
        
        cut_data_B.append(data_B[masks_A[i]])

    return cut_data_A, cut_data_B


###############################################################################
###############################################################################

########################################
### SETTING FOLDERS, FILES, ETC   ######
########################################

input_cs82 = sys.argv[1]

sdss_csv = sys.argv[2]

folder_path = os.path.dirname(input_cs82).replace('/AFTER_PSFEX','')

filename = os.path.basename(input_cs82)

ioadd.create_folder(folder_path + '/plots_sdss/')

# Choose plot fields:

if filename.find('deV') > 0:

    fieldnames_sdss = ['deVRad_i','deVAB_i','deVPhi_i']
    fieldnames_cs82 = ['SPHEROID_REFF_WORLD','SPHEROID_ASPECT_WORLD','THETAMODEL_IMAGE']

elif filename.find('exp') > 0:

    fieldnames_sdss = ['expRad_i','expAB_i','expPhi_i']
    fieldnames_cs82 = ['DISK_SCALE_WORLD','DISK_ASPECT_WORLD','DISK_THETA_IMAGE']

else:

    print "There's something wrong with the input file"
    sys.exit(1)


######################################################################################
### IMPORTING THE DATA AND DOING THE NECESSARY MATCHING, BASIC CUTS AND REORDERING ###
######################################################################################

# ===================== Transformation from the SDSS skyserver CSV to a FITS catalog


sdss_csv = sys.argv[2]

fieldnames=['objID','ra','dec','psffwhm_i','deVRad_i','deVRadErr_i','deVAB_i','deVABErr_i','deVPhi_i','deVMag_i','deVMagErr_i','deVFlux_i','deVFluxIvar_i',
'expRad_i','expRadErr_i','expAB_i','expABErr_i','expPhi_i','expMag_i','expMagErr_i','expFlux_i','expFluxIvar_i','cModelMag_i','cModelMagErr_i',
'cModelFlux_i','cModelFluxIvar_i','modelMagErr_i','modelFlux_i','modelFluxIvar_i','run','rerun','camcol','field','type','modelMag_u','modelMag_g',
'modelMag_r','modelMag_i','modelMag_z']

sdss_dict = ascd.dict_from_csv(sdss_csv,fieldnames, header_lines=1, delimiter=',',dialect='excel')

sdss_tbhdu = fd.dict_to_tbHDU(sdss_dict)

sdss_data = sdss_tbhdu.data


# =============================== Getting CS82 Catalog Information

  
cs82_tbhdu = pyfits.open(input_cs82,ignore_missing_end=True,memmap=True)[2]

cs82_cut_hdu = fd.sample_entries(cs82_tbhdu,CLASS_STAR=(0,0.95), MAG_AUTO=(18.5,20.5), MAGERR_AUTO=(0,0.2172))

cs82_data = cs82_cut_hdu.data

print len(cs82_data)



# ============================= Doing a Nearest Neighbour Match and Getting Full Reordered cs82_data.

radius = 0.5*0.000277 # degrees ~ 1 arc-second

cs82_full, sdss_data = matching_pipeline(cs82_data, sdss_data, ['alpha_J2000','delta_J2000'], ['ra','dec'], radius)
'''
sdss_ra = sdss_data.field('ra').real
sdss_dec = sdss_data.field('dec').real
sdss_radec = zip(sdss_ra,sdss_dec)

cs82_ra = cs82_data.field('alpha_J2000')
cs82_dec = cs82_data.field('delta_J2000')
cs82_radec = zip(cs82_ra, cs82_dec)


near_neigh = np.array(tbm.nearest_neighbour(sdss_radec,cs82_radec))

true_match = np.array([ float(_d) < float(radius) for _d in near_neigh[:,1]])

sdss_data = sdss_data[true_match]
near_neigh = near_neigh[true_match]

mask = near_neigh[:,0].astype(int)

cs82_full = cs82_data[mask]
'''
# I haven't made sure that the match was unique from the cs82 side, 
# i.e. two different sdss objects may be matching the same cs82 one.
# But in view of the cut-off radius chosen, this is a completely secondary 
# effect.

'''
sdss_ra = sdss_data.field('ra').real
sdss_dec = sdss_data.field('dec').real
sdss_radec = zip(sdss_ra,sdss_dec)

cs82_ra = cs82_full.field('alpha_J2000')
cs82_dec = cs82_full.field('delta_J2000')
cs82_radec = zip(cs82_ra, cs82_dec)

near_neigh3 = np.array(tbm.nearest_neighbour(cs82_radec,sdss_radec))


print near_neigh3[:,1].max()
'''
print len(cs82_full)

# Keeping only FLAGS = 0 for the analysis

mask_flags = (cs82_full.field('FLAGS') == 0)

cs82_full = cs82_full[mask_flags]
sdss_data = sdss_data[mask_flags]

print len(sdss_data)

'''
# Another ds9 mask

image_name = '/home/brunomor/Documents/Trabalho/S82p28m_y.V2.7A.swarp.cut.fits'

hdu_image = pyfits.open(image_name)

image_header = hdu_image[0].header

sdss_ra = sdss_data.field('ra').real
sdss_dec = sdss_data.field('dec').real

sdss_ra_list = sdss_ra.tolist()
sdss_dec_list = sdss_dec.tolist()

pixcrd = []
xlist = []
ylist = []

for i in range(len(sdss_ra)):
    pixcrd.append(wcscnv.radec2xy(image_header,sdss_ra_list[i],sdss_dec_list[i]))
    xlist.append(pixcrd[i][0])
    ylist.append(pixcrd[i][1])

ascd.write_ds9cat(xlist,ylist,30,'circle','yellow', outputfile='boss_new_matching.reg')


# Another ds9 mask


cs82_ra = cs82_full.field('alpha_J2000').real
cs82_dec = cs82_full.field('delta_J2000').real

cs82_ra_list = cs82_ra.tolist()
cs82_dec_list = cs82_dec.tolist()

pixcrd = []
xlist = []
ylist = []

for i in range(len(sdss_ra)):
    pixcrd.append(wcscnv.radec2xy(image_header,cs82_ra_list[i],cs82_dec_list[i]))
    xlist.append(pixcrd[i][0])
    ylist.append(pixcrd[i][1])

ascd.write_ds9cat(xlist,ylist,40,'circle','green',outputfile='cs82_new_matching.reg')
'''

# Selecting the fields to be used for the comparison

# bla bla

################################################
##########    SEVERAL CUTS & PLOTS    ##########
################################################

# SECTION CUTS - IMPROVED Feb/28/2012

# Magnitude cuts 

cs82_mag, sdss_mag = create_matched_subsets(cs82_full,sdss_data,'MAG_AUTO',[['.',19],[19,20],[20,21]])


# Size cuts with CS82 FWHM

s_cut = 8

cs82_size, sdss_size = create_matched_subsets(cs82_mag[0],sdss_mag[0],'FWHM_IMAGE',[['.',s_cut],[s_cut,'.']])


# Size cuts with SDSS PSF FWHM

psf_cut = 0.75

sdss_psf, cs82_psf = create_matched_subsets(sdss_data,cs82_full,'psffwhm_i',[['.',psf_cut],[psf_cut,'.']])

# WARNING: There was an issue with the ingestion of SDSS data here, the field 
# was complex and needed to be transformed back to real.

###########################################################################################################

# SECTION PLOTS - CAN I FORMULATE A SUFFICIENTLY GENERAL FUNCTION TO SIMPLIFY THIS PART?

plot_lims = [[[0,5],[0,1],[-90,90]],[[0,5],[0,1],[0,180]]]

plot_titles = ['Comparison 1 - Effective Radius', 'Comparison 2 - Aspect Ratio', 'Comparison 3 - Theta']

units = [[' (arcsec)',' (arcsec)'],['',''],[' (deg)',' (deg)']]


for i in range(len(fieldnames_cs82)):

    fit_x = np.arange(plot_lims[1][i][0],plot_lims[1][i][1],0.01)

    if i == 0:

        if fieldnames_cs82[i] == 'SPHEROID_REFF_WORLD':
            j=0
        else:
            j=1
        

        color_mag_data = [[cs82_mag[2].field(fieldnames_cs82[i])*3600*(1.68**j),sdss_mag[2].field(fieldnames_sdss[i]).real],[cs82_mag[1].field(fieldnames_cs82[i])*3600*(1.68**j),sdss_mag[1].field(fieldnames_sdss[i]).real],[cs82_mag[0].field(fieldnames_cs82[i])*3600*(1.68**j),sdss_mag[0].field(fieldnames_sdss[i]).real],[fit_x,fit_x]]


        cs82_size_data = [[cs82_size[1].field(fieldnames_cs82[i])*3600*(1.68**j),sdss_size[1].field(fieldnames_sdss[i]).real],[cs82_size[0].field(fieldnames_cs82[i])*3600*(1.68**j),sdss_size[0].field(fieldnames_sdss[i]).real],[fit_x,fit_x]]


        sdss_size_data = [[cs82_psf[1].field(fieldnames_cs82[i])*3600*(1.68**j),sdss_psf[1].field(fieldnames_sdss[i]).real],[cs82_psf[0].field(fieldnames_cs82[i])*3600*(1.68**j),sdss_psf[0].field(fieldnames_sdss[i]).real],[fit_x,fit_x]]

        print fieldnames_cs82[i]


    elif i == 2:


        color_mag_data = [[cs82_mag[2].field(fieldnames_cs82[i]),sdss_mag[2].field(fieldnames_sdss[i]).real],[cs82_mag[1].field(fieldnames_cs82[i]),sdss_mag[1].field(fieldnames_sdss[i]).real],[cs82_mag[0].field(fieldnames_cs82[i]),sdss_mag[0].field(fieldnames_sdss[i]).real],[fit_x-90,fit_x]]


        cs82_size_data = [[cs82_size[1].field(fieldnames_cs82[i]),sdss_size[1].field(fieldnames_sdss[i]).real],[cs82_size[0].field(fieldnames_cs82[i]),sdss_size[0].field(fieldnames_sdss[i]).real],[fit_x-90,fit_x]]


        sdss_size_data = [[cs82_psf[1].field(fieldnames_cs82[i]),sdss_psf[1].field(fieldnames_sdss[i]).real],[cs82_psf[0].field(fieldnames_cs82[i]),sdss_psf[0].field(fieldnames_sdss[i]).real],[fit_x-90,fit_x]]

        print fieldnames_cs82[i]


    else:

        color_mag_data = [[cs82_mag[2].field(fieldnames_cs82[i]),sdss_mag[2].field(fieldnames_sdss[i]).real],[cs82_mag[1].field(fieldnames_cs82[i]),sdss_mag[1].field(fieldnames_sdss[i]).real],[cs82_mag[0].field(fieldnames_cs82[i]),sdss_mag[0].field(fieldnames_sdss[i]).real],[fit_x,fit_x]]


        cs82_size_data = [[cs82_size[1].field(fieldnames_cs82[i]),sdss_size[1].field(fieldnames_sdss[i]).real],[cs82_size[0].field(fieldnames_cs82[i]),sdss_size[0].field(fieldnames_sdss[i]).real],[fit_x,fit_x]]


        sdss_size_data = [[cs82_psf[1].field(fieldnames_cs82[i]),sdss_psf[1].field(fieldnames_sdss[i]).real],[cs82_psf[0].field(fieldnames_cs82[i]),sdss_psf[0].field(fieldnames_sdss[i]).real],[fit_x,fit_x]]

        print fieldnames_cs82[i]

    
    pltemp.plot_data(color_mag_data, ['.','.','.','-'], ['b','g','r','k'], ['20 < MAG_AUTO < 21', '19 < MAG_AUTO < 20', 'MAG_AUTO < 19',''], 'CS82 - ' + fieldnames_cs82[i] + units[i][0], 'SDSS - ' + fieldnames_sdss[i] + units[i][1], plot_lims[0][i], plot_lims[1][i], plot_titles[i])
    ioadd.create_folder(folder_path + '/plots_sdss/')
    pl.savefig(folder_path + '/plots_sdss/' + fieldnames_sdss[i] + '_mags.png')    
    pl.clf()

    pltemp.plot_data(cs82_size_data, ['.','.','-'], ['b','r','k'], ['CS82 FWHM_IMAGE > '+ str(s_cut*0.186)+'"','CS82 FWHM_IMAGE < '+ str(s_cut*0.186)+'"',''], 'CS82 - ' + fieldnames_cs82[i] + units[i][0], 'SDSS - ' + fieldnames_sdss[i] + units[i][1], plot_lims[0][i], plot_lims[1][i], plot_titles[i])
    pl.savefig(folder_path + '/plots_sdss/' + fieldnames_sdss[i] + '_cs82_fwhm.png')    
    pl.clf()

    pltemp.plot_data(sdss_size_data, ['.','.','-'], ['b','r','k'], ['SDSS psffwhm_i > '+ str(psf_cut)+'"','SDSS psffwhm_i < '+ str(psf_cut)+'"',''], 'CS82 - ' + fieldnames_cs82[i] + units[i][0], 'SDSS - ' + fieldnames_sdss[i] + units[i][1], plot_lims[0][i], plot_lims[1][i], plot_titles[i])
    ioadd.create_folder(folder_path + '/plots_sdss/')
    pl.savefig(folder_path + '/plots_sdss/' + fieldnames_sdss[i] + '_sdss_psf.png')    
    pl.clf()


# Masters' Figure 11 Mag Color:

print j

rel_reff_21 = [cs82_mag[2].field(fieldnames_cs82[0])*3600*1.68**j,(sdss_mag[2].field(fieldnames_sdss[0]).real-cs82_mag[2].field(fieldnames_cs82[0])*3600*1.68**j)/(cs82_mag[2].field(fieldnames_cs82[0])*3600*1.68**j)]

rel_reff_20 =[cs82_mag[1].field(fieldnames_cs82[0])*3600*1.68**j,(sdss_mag[1].field(fieldnames_sdss[0]).real-cs82_mag[1].field(fieldnames_cs82[0])*3600*1.68**j)/(cs82_mag[1].field(fieldnames_cs82[0])*3600*1.68**j)]

rel_reff_19 =[cs82_mag[0].field(fieldnames_cs82[0])*3600*1.68**j,(sdss_mag[0].field(fieldnames_sdss[0]).real-cs82_mag[0].field(fieldnames_cs82[0])*3600*1.68**j)/(cs82_mag[0].field(fieldnames_cs82[0])*3600*1.68**j)]

rel_reff = [rel_reff_21, rel_reff_20, rel_reff_19, [np.arange(0, 5, 0.01),0*np.arange(0, 5, 0.01)]]


pltemp.plot_data(rel_reff, ['.','.','.','--'], ['b','g','r','k'], ['20 < MAG_AUTO < 21 / median: ' + str('%.2f' % np.median(rel_reff_21[1])), '19 < MAG_AUTO < 20 / median: '+ str('%.2f' % np.median(rel_reff_20[1])), 'MAG_AUTO < 19 / median: '+ str('%.2f' % np.median(rel_reff_19[1])),''], 'CS82 - ' + fieldnames_cs82[0] + units[0][0], '(reff_sdss - reff_cs82)/reff_cs82', [0, 4], [-1, 5], 'Relative Effective Radius Difference SDSS x CS82')

pl.savefig(folder_path + '/plots_sdss/fig_11_'+fieldnames_sdss[0][0:3]+'.png')    
pl.clf()


# Masters' Figure 11 Size:

print j

rel_reff_big = [cs82_size[1].field(fieldnames_cs82[0])*3600*1.68**j,(sdss_size[1].field(fieldnames_sdss[0]).real-cs82_size[1].field(fieldnames_cs82[0])*3600*1.68**j)/(cs82_size[1].field(fieldnames_cs82[0])*3600*1.68**j)]

rel_reff_small = [cs82_size[0].field(fieldnames_cs82[0])*3600*1.68**j,(sdss_size[0].field(fieldnames_sdss[0]).real-cs82_size[0].field(fieldnames_cs82[0])*3600*1.68**j)/(cs82_size[0].field(fieldnames_cs82[0])*3600*1.68**j)]


rel_reff_size = [rel_reff_big, rel_reff_small, [np.arange(0, 5, 0.01),0*np.arange(0, 5, 0.01)]]


pltemp.plot_data(rel_reff_size, ['.','.','--'], ['b','r','k'], ['CS82 FWHM_IMAGE > '+ str(s_cut*0.186)+'" / median: ' + str('%.2f' % np.median(rel_reff_big[1])), 'CS82 FWHM_IMAGE < '+ str(s_cut*0.186)+'" / median: '+ str('%.2f' % np.median(rel_reff_small[1])),''], 'CS82 - ' + fieldnames_cs82[0] + units[0][0], '(reff_sdss - reff_cs82)/reff_cs82', [0, 4], [-1, 5], 'Relative Effective Radius Difference SDSS x CS82')

pl.savefig(folder_path + '/plots_sdss/fig_11_size_'+fieldnames_sdss[0][0:3]+'.png')    
pl.clf()



# Half-light Radius Figures CS82:

print j

hlrad_21 = [cs82_mag[2].field(fieldnames_cs82[0])*3600*1.68**j,cs82_mag[2].field('FLUX_RADIUS')*0.186]

hlrad_20 =[cs82_mag[1].field(fieldnames_cs82[0])*3600*1.68**j,cs82_mag[1].field('FLUX_RADIUS')*0.186]

hlrad_19 =[cs82_mag[0].field(fieldnames_cs82[0])*3600*1.68**j,cs82_mag[0].field('FLUX_RADIUS')*0.186]

hlrad = [hlrad_21, hlrad_20, hlrad_19, [np.arange(0, 5, 0.01),np.arange(0, 5, 0.01)]]


pltemp.plot_data(hlrad, ['.','.','.','--'], ['b','g','r','k'], ['20 < MAG_AUTO < 21', '19 < MAG_AUTO < 20', 'MAG_AUTO < 19',''], 'CS82 - ' + fieldnames_cs82[0] + units[0][0], 'CS82 - FLUX_RADIUS [0.5]' + units[0][0], [0, 4], [0, 4], 'CS82 Half-Light Radius x Effective Radius')

pl.savefig(folder_path + '/plots_sdss/hlradius_'+fieldnames_sdss[0][0:3]+'.png')    
pl.clf()


# Half-light Radius Figures SDSS:

print j

hlrad_21 = [sdss_mag[2].field(fieldnames_sdss[0]).real,cs82_mag[2].field('FLUX_RADIUS')*0.186]

hlrad_20 =[sdss_mag[1].field(fieldnames_sdss[0]).real,cs82_mag[1].field('FLUX_RADIUS')*0.186]

hlrad_19 =[sdss_mag[0].field(fieldnames_sdss[0]).real,cs82_mag[0].field('FLUX_RADIUS')*0.186]

hlrad = [hlrad_21, hlrad_20, hlrad_19, [np.arange(0, 5, 0.01),np.arange(0, 5, 0.01)]]


pltemp.plot_data(hlrad, ['.','.','.','--'], ['b','g','r','k'], ['20 < MAG_AUTO < 21', '19 < MAG_AUTO < 20', 'MAG_AUTO < 19',''], 'CS82 - ' + fieldnames_cs82[0] + units[0][0], 'CS82 - FLUX_RADIUS [0.5]' + units[0][0], [0, 4], [0, 4], 'CS82 Half-Light Radius x SDSS Effective Radius')

pl.savefig(folder_path + '/plots_sdss/hlradius_sdss_'+fieldnames_sdss[0][0:3]+'.png')    
pl.clf()

# The following piece of code should be placed elsewhere, decide later where.

'''
# =============================== Creating the masking file for the full cs82 mag 21 S/N 5 CLASS_STAR < 0.95


# Transform from radec to xy:

image_name = '/home/brunomor/Documents/Trabalho/S82p28m_y.V2.7A.swarp.cut.fits'

hdu_image = pyfits.open(image_name)

image_header = hdu_image[0].header

cs82_ra_list = cs82_ra.tolist()
cs82_dec_list = cs82_dec.tolist()

pixcrd = []
xlist = []
ylist = []

for i in range(len(cs82_ra)):
    pixcrd.append(wcscnv.radec2xy(image_header,cs82_ra_list[i],cs82_dec_list[i]))
    xlist.append(pixcrd[i][0])
    ylist.append(pixcrd[i][1])

ascd.write_ds9cat(xlist,ylist,10,'circle','red', outputfile = 'cs82_full_cuts_objects.reg')


# =============================== Creating the Mask File for the full boss/sdss data


# Transform from radec to xy:

image_name = '/home/brunomor/Documents/Trabalho/S82p28m_y.V2.7A.swarp.cut.fits'

hdu_image = pyfits.open(image_name)

image_header = hdu_image[0].header

sdss_ra_list = sdss_ra.tolist()
sdss_dec_list = sdss_dec.tolist()

pixcrd = []
xlist = []
ylist = []

for i in range(len(sdss_ra)):
    pixcrd.append(wcscnv.radec2xy(image_header,sdss_ra_list[i],sdss_dec_list[i]))
    xlist.append(pixcrd[i][0])
    ylist.append(pixcrd[i][1])

ascd.write_ds9cat(xlist,ylist,30,'circle','yellow', outputfile = 'boss_full_objects.reg')
'''
