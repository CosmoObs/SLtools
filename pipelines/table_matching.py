#!/usr/bin/env python

""" Module to deal with tables matching """

import sys;

import numpy as np;

import sltools;
from sltools.catalog import fits_data,ascii_data;

def nearest_neighbour(centroids_A, centroids_truth):
    """ Function to compute and return the set of nearest point
        
    This function computes for each entry of 'centroids_A', the
    nearest point in 'centroids_B'. Then, it does the same for
    'centroids_B', computing the nearest 'centroids_A' to each
    "B" entry.
    It is returned, respective to each entry, a "index" and a
    "distante" measure correspondig to each identified nearest
    point.
    
    Input:
     - centroids_A : [(x_A0,y_A0),(x_A1,y_A1), ]
     - centroids_B : [(x_B0,y_B0),(x_B1,y_B1), ]
    
    Output:
     - [(index_BA0,distance_BA0), ]
     - [(index_AB0,distance_AB0), ]

    """
    
    centroids_B = centroids_truth;
    
    length_A = len(centroids_A);
    length_B = len(centroids_B);
    
    x_A,y_A = zip(*centroids_A);
    x_B,y_B = zip(*centroids_B);
    
    Mat = np.zeros((length_A,length_B));
    
    for i in xrange(length_A):
        Mat[i] = (x_A[i] - x_B[:])**2 + (y_A[i] - y_B[:])**2;
    
    # Distance(minimum) from all points in "_A" to each point in "_B".
    # So, the following vector, "_Ato_B", means "minimum distance of 
    # all points in vector _A to each point in _B, and respective (_A)
    # vector index". Size of both vector should be equal to input "_B"
    # data.('centroids_B')
#    dist_min_AtoB = np.sqrt(np.amin(Mat,axis=0));
#    indx_min_AtoB = np.argmin(Mat,axis=0);
    
#    ind_dist_AtoB = zip(indx_min_AtoB,dist_min_AtoB);

    # Again here for "_B". All the points in "_B" near(est) to each
    # point in "_A". Index of the nearest point in "_B" and the distance
    # to that point. Size iqual to "_A"-input.(len('centroids_A'))
    dist_min_BtoA = np.sqrt(np.amin(Mat,axis=1));
    indx_min_BtoA = np.argmin(Mat,axis=1);
    
    ind_dist_BtoA = zip(indx_min_BtoA,dist_min_BtoA);
    
#    return (ind_dist_BtoA,ind_dist_AtoB);
    return (ind_dist_BtoA);

# ---

def match_positions(tbhdu_A, tbhdu_truth, radius):
    """ Check if positions in tables match within a given radius
    
    "tbhdu"s 'x' and 'y' fields are expected for the processing
    
    Input:
     - tbhdu_A (main)
     - tbhdu_B (reference)
     - radius : float
    
    Output:
     -  
    """
    
    tbhdu_B = tbhdu_truth;
    
    tbname_A = tbhdu_A.name;
    tbname_B = tbhdu_B.name;
    try:
        x_A = tbhdu_A.data.field('X');
        y_A = tbhdu_A.data.field('Y');
    except:
        x_A = tbhdu_A.data.field('X_IMAGE');
        y_A = tbhdu_A.data.field('Y_IMAGE');
    try:
        x_B = tbhdu_B.data.field('X');
        y_B = tbhdu_B.data.field('Y');
    except:
        x_B = tbhdu_B.data.field('X_IMAGE');
        y_B = tbhdu_B.data.field('Y_IMAGE');
    
    centroids_A = zip(x_A,y_A);
    centroids_B = zip(x_B,y_B);
    
    try:
        len_A = len(centroids_A);
        len_B = len(centroids_B);
    except:
        return False;
    
    A_near_neibors = nearest_neighbour(centroids_A,centroids_B);
    
    A_nn_indx,A_nn_dist = zip(*A_near_neibors);
    
    A_match_neibors = [ _d < radius for _d in A_nn_dist ];
    
    A_neibor_BX = [ x_B[i] for i in A_nn_indx ];
    A_neibor_BY = [ y_B[i] for i in A_nn_indx ];
    
    # Alias:
    naa = np.asarray;

    tbname = tbname_B+"_points_near_to_"+tbname_A;
    dict_Aout = {'X':naa(x_A),'Y':naa(y_A),'XY_MATCH':naa(A_match_neibors),'X_NEAR':naa(A_neibor_BX),'Y_NEAR':naa(A_neibor_BY)};
    new_tbhdu_A = fits_data.dict_to_tbHDU(dict_Aout,tbname);
    
#    A_near_neibors, B_near_neibors = nearest_neighbour(centroids_A,centroids_B);
#    B_nn_indx,B_nn_dist = zip(*B_near_neibors);
#    B_match_neibors = [ _d < radius for _d in B_nn_dist ];
#    B_neibor_AX = [ x_A[i] for i in B_nn_indx ];
#    B_neibor_AY = [ y_A[i] for i in B_nn_indx ];
#    tbname = tbname_A+"_points_near_to_"+tbname_B;
#    dict_Bout = {'X':naa(x_B),'Y':naa(y_B),'XY_MATCH':naa(B_match_neibors),'X_NEAR':naa(B_neibor_AX),'Y_NEAR':naa(B_neibor_AY)};
#    new_tbhdu_B = fits_data.dict_to_tbHDU(dict_Bout,tbname);
#    return (new_tbhdu_A,new_tbhdu_B);

    return (new_tbhdu_A);

# ---

#def run(filename_A,filename_B,radius=1,output_A='tableA',output_B='tableB', write_FITS=False):
def run(filename_A,filename_B,radius=1,output_A='tableA', write_FITS=False):
    """
    Input:
     - filename_A : DS9 region file
     - filename_B : DS9 region file
     
    Output:
    
    """
    # Open two tables, for now, ds9 region files:
    #
    DA_in = ascii_data.read_ds9cat(filename_A);
    DB_in = ascii_data.read_ds9cat(filename_B);
    
    # Translate them to table HDU (FITS tables Header Table Unit)
    tbhdu_A = fits_data.dict_to_tbHDU(DA_in,tbname=filename_A);
    tbhdu_B = fits_data.dict_to_tbHDU(DB_in,tbname=filename_B);
    
#    tbhdu_A_wneib,tbhdu_B_wneib = match_positions(tbhdu_A,tbhdu_B,radius);
    tbhdu_A_wneib = match_positions(tbhdu_A,tbhdu_B,radius);
    
    return tbhdu_A_wneib;

# ==========================

if __name__ == '__main__' :

    # Initializing code, command-line options..
    #
    import optparse;
    
    usage="\n  %prog [options] <ds9_CheckTable.reg> <ds9_TruthTable.reg>"
    parser = optparse.OptionParser(usage=usage);

    parser.add_option('-r',
                        dest='radius', default=10,
                        help='Distance for the objects matching in both tables [10]');
    parser.add_option('-o',
                        dest='filename', default='matching_points',
                        help='Output tables (DS9 and FITS) filename [matching_points]');
    parser.add_option('--write_FITS', action='store_true',
                        dest='write_fits', default=False,
                        help='If used, write output FITS files (tables)');
                        
    (opts,args) = parser.parse_args();

    radius = opts.radius;
    outfilename = opts.filename;
    write_FITS = opts.write_fits;
    
    if len(args) < 2 :
        parser.print_help();
        sys.exit(1);
        
    regionfile_A = args[0];
    regionfile_B = args[1];

    tbhdu_A = run(regionfile_A, regionfile_B, radius)
    
    if write_FITS:
        tbhdu_A.writeto(outfilename+'.fit');
    
    # Set the dictionary to output to a ds9 region file
    # the points that matched with the truth table will
    # be written as yellow in ds9 region file.
    
    tbdata = tbhdu_A.data;

    l_len = range( int(tbhdu_A.header['NAXIS2']) );
    
    TF = tbdata.field('XY_MATCH');
    
    x = [ tbdata.field('X_NEAR')[i] for i in l_len if TF[i] ];
    y = [ tbdata.field('Y_NEAR')[i] for i in l_len if TF[i] ];
    marker = [ 'circle' for i in l_len if TF[i] ];
    color = [ 'green' for i in l_len if TF[i] ];
    size = [ 30 for i in l_len if TF[i] ];
    
    ascii_data.write_ds9cat(x, y, size, marker, color, outputfile=outfilename+'.reg')
    
    print >> sys.stdout, "Done.";
    sys.exit(0);
