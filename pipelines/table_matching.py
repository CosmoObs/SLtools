#!/usr/bin/env python

""" Module to deal with tables matching """

import sys;

import numpy as np;

import sltools;
from sltools.catalog import fits_data,ascii_data;

def nearest_neighbour(centroids_A, centroids_B):
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
    dist_min_AtoB = np.sqrt(np.amin(Mat,axis=0));
    indx_min_AtoB = np.argmin(Mat,axis=0);
    
    # Again here for "_B". All the points in "_B" near(est) to each
    # point in "_A". Index of the nearest point in "_B" and the distance
    # to that point. Size iqual to "_A"-input.(len('centroids_A'))
    dist_min_BtoA = np.sqrt(np.amin(Mat,axis=1));
    indx_min_BtoA = np.argmin(Mat,axis=1);
    
    ind_dist_BtoA = zip(indx_min_BtoA,dist_min_BtoA);
    ind_dist_AtoB = zip(indx_min_AtoB,dist_min_AtoB);
    
    return (ind_dist_BtoA,ind_dist_AtoB);

# ---

def match_positions(tbhdu_A, tbhdu_B, radius):
    """ Check if positions in tables match within a given radius
    
    "tbhdu"s 'x' and 'y' fields are expected for the processing
    
    Input:
     - tbhdu_A (main)
     - tbhdu_B (reference)
     - radius : float
    
    Output:
     -  
    """
    
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
    
    A_near_neibors, B_near_neibors = nearest_neighbour(centroids_A,centroids_B);
    
    A_nn_indx,A_nn_dist = zip(*A_near_neibors);
    B_nn_indx,B_nn_dist = zip(*B_near_neibors);

    A_match_neibors = [ _d < radius for _d in A_nn_dist ];
    B_match_neibors = [ _d < radius for _d in B_nn_dist ];
    
    A_neibor_BX = [ x_B[i] for i in A_nn_indx ];
    A_neibor_BY = [ y_B[i] for i in A_nn_indx ];
    
    B_neibor_AX = [ x_A[i] for i in B_nn_indx ];
    B_neibor_AY = [ y_A[i] for i in B_nn_indx ];
    
    naa = np.asarray;

    tbname = tbname_B+"_points_near_to_"+tbname_A;
    dict_Aout = {'X':naa(x_A),'Y':naa(y_A),'XY_MATCH':naa(A_match_neibors),'X_NEAR':naa(A_neibor_BX),'Y_NEAR':naa(A_neibor_BY)};
    new_tbhdu_A = fits_data.dict_to_tbHDU(dict_Aout,tbname);
    
    tbname = tbname_A+"_points_near_to_"+tbname_B;
    dict_Bout = {'X':naa(x_B),'Y':naa(y_B),'XY_MATCH':naa(B_match_neibors),'X_NEAR':naa(B_neibor_AX),'Y_NEAR':naa(B_neibor_AY)};
    new_tbhdu_B = fits_data.dict_to_tbHDU(dict_Bout,tbname);
    
    return (new_tbhdu_A,new_tbhdu_B);

# ---

def run(filename_A,filename_B,radius=1,output_A='tableA',output_B='tableB', write_FITS=False):
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
    
    tbhdu_A_wneib,tbhdu_B_wneib = match_positions(tbhdu_A,tbhdu_B,radius);
    
    DA_out = DA_in;
    TF = tbhdu_A_wneib.data.field('XY_MATCH');
    _len = int(tbhdu_A_wneib.header['NAXIS2']);
    _ctmp = tbhdu_A_wneib.data.field('X_NEAR');
    DA_out['x'].extend( [ _ctmp[i] for i in range(_len) if TF[i] ] );
    _ctmp = tbhdu_A_wneib.data.field('Y_NEAR');
    DA_out['y'].extend( [ _ctmp[i] for i in range(_len) if TF[i] ] );
    DA_out['color'].extend( [ 'green' for i in range(_len) if TF[i] ] );
    DA_out['marker'].extend( [ 'circle' for i in range(_len) if TF[i] ] );
    DA_out['size'].extend( [ 30 for i in range(_len) if TF[i] ] );
    
    DB_out = DB_in;
    TF = tbhdu_B_wneib.data.field('XY_MATCH');
    _len = int(tbhdu_B_wneib.header['NAXIS2']);
    _ctmp = tbhdu_B_wneib.data.field('X_NEAR');
    DB_out['x'].extend( [ _ctmp[i] for i in range(_len) if TF[i] ] );
    _ctmp = tbhdu_B_wneib.data.field('Y_NEAR');
    DB_out['y'].extend( [ _ctmp[i] for i in range(_len) if TF[i] ] );
    DB_out['color'].extend( [ 'green' for i in range(_len) if TF[i] ] );
    DB_out['marker'].extend( [ 'circle' for i in range(_len) if TF[i] ] );
    DB_out['size'].extend( [ 30 for i in range(_len) if TF[i] ] );
    
    ascii_data.write_ds9cat(x=DA_out['x'],y=DA_out['y'],size=DA_out['size'],marker=DA_out['marker'],color=DA_out['color'],outputfile=output_A+'.reg')
    ascii_data.write_ds9cat(x=DB_out['x'],y=DB_out['y'],size=DB_out['size'],marker=DB_out['marker'],color=DB_out['color'],outputfile=output_B+'.reg')
    
    if write_FITS:
        tbhdu_A_wneib.writeto(output_A+'.fit');
        tbhdu_B_wneib.writeto(output_B+'.fit');
    
    return
    

# ==========================

if __name__ == '__main__' :

    # Initializing code, command-line options..
    #
    import optparse;
    print "";
    
    usage="Usage:  %prog [options] <ds9regionfile_A.reg> <ds9regionfile_A.reg>"
    parser = optparse.OptionParser(usage=usage);

    parser.add_option('--radius',
                        dest='radius', default=1,
                        help='Distance for the objects matching in both tables');
    parser.add_option('--outname_A',
                        dest='outname_A', default='table_new_A',
                        help='Output tables (DS9 and FITS) file rootnames["table_new_A"]');
    parser.add_option('--outname_B',
                        dest='outname_B', default='table_new_B',
                        help='Output tables (DS9 and FITS) file rootnames["table_new_B"]');
    parser.add_option('--write_FITS', action='store_true',
                        dest='write_fits', default=False,
                        help='If used, write output FITS files (tables)');
                        
    (opts,args) = parser.parse_args();

    radius = opts.radius;
    outname_A = opts.outname_A;
    outname_B = opts.outname_B;
    write_FITS = opts.write_fits;
    
    if len(args) < 2 :
        parser.print_help();
        sys.exit(1);
        
    regionfile_A = args[0];
    regionfile_B = args[1];

    run(regionfile_A, regionfile_B, radius, outname_A, outname_B, write_FITS)

    print >> sys.stdout, "Done.";
    sys.exit(0);
