#!/usr/bin/env python

""" Module to deal with tables matching """

import sys;
import os;
import logging;

import numpy as np;
import pyfits;

import sltools;
from sltools.catalog import fits_data,ascii_data;
from sltools.io import log;

def nearest_neighbour(centroids_A, centroids_B):
    """ Function to compute and return the set of nearest point
        
    This function computes for each entry of 'centroids_A', the
    nearest point in 'centroids_B'.
    
    Is returned, respective to each entry, the index and the
    distance measured correspondig to each identified nearest
    point.
    
    Input:
     - centroids_A : [(x_A0,y_A0),(x_A1,y_A1), ]
     - centroids_B : [(x_B0,y_B0),(x_B1,y_B1), ]
    
    Output:
     - [(index_BtoA0,distance_BtoA0), ]

    ---
    """
    
    length_A = len(centroids_A);
    length_B = len(centroids_B);
    
    x_A,y_A = zip(*centroids_A);
    x_B,y_B = zip(*centroids_B);
    
    Mat = np.zeros((length_A,length_B));
    
    for i in xrange(length_A):
        Mat[i] = (x_A[i] - x_B[:])**2 + (y_A[i] - y_B[:])**2;
    
    dist_min_BtoA = np.sqrt(np.amin(Mat,axis=1));
    indx_min_BtoA = np.argmin(Mat,axis=1);
    
    ind_dist_BtoA = zip(indx_min_BtoA,dist_min_BtoA);
    

    return (ind_dist_BtoA);

# ---

def match_positions(tbhdu_check, tbhdu_truth, radius):
    """
    Check if positions in tables match within a given radius
    
    "tbhdu"s 'x' and 'y' fields are expected for the processing.
    
    Input:
     - tbhdu_check : tbHDU
        Binary table HDU, from FITS catalog
     - tbhdu_truth : tbHDU
        Binary table HDU, from FITS catalog
     - radius : float
        Distance parameters for objects position matching
    
    Output:
     - new_tbhdu : tbHDU
        New table (check) with the matching information added
        
    ---
    """
    
    logging.debug("Input params (tbC,tbT,radius): %s,%s,%s" % (tbhdu_check.name,tbhdu_truth.name,radius));
    
    tbhdu_A = tbhdu_check;
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
    [ logging.debug("Checktable nearest neighbours (indx,dist): %s",each) for each in A_near_neibors ];
    
    A_nn_indx,A_nn_dist = zip(*A_near_neibors);
    
    logging.info("Matching radius: %s", radius);
    A_match_neibors = [ float(_d) < float(radius) for _d in A_nn_dist ];
    logging.debug("Checktable matched objects: %s", A_match_neibors);
    
    A_neibor_BX = [ x_B[i] for i in A_nn_indx ];
    A_neibor_BY = [ y_B[i] for i in A_nn_indx ];
    
    # Alias:
    naa = np.asarray;

    tbname = tbname_A+"_matched_objs";
    dict_Aout = {'X':naa(x_A),'Y':naa(y_A),'XY_MATCH':naa(A_match_neibors),'X_NEAR':naa(A_neibor_BX),'Y_NEAR':naa(A_neibor_BY)};
    new_tbhdu_A = fits_data.dict_to_tbHDU(dict_Aout,tbname);
    
    Ntrue = A_match_neibors.count(True);
    Nfalse = A_match_neibors.count(False);
    Ntotal_truth = tbhdu_truth.header['naxis2'];
    
    new_tbhdu_A.header.update('hierarch truthtbname',tbhdu_truth.name,"Name of Truth table used");
    new_tbhdu_A.header.update('hierarch ntottruth',Ntotal_truth,"Total number of arcs in Truth table");
    new_tbhdu_A.header.update('ntrue',Ntrue,"Number of found objs in Check table");
    new_tbhdu_A.header.update('nfalse',Nfalse,"Number of not-found objs in Check table");

    
    return (new_tbhdu_A);

# ---

def run(tbhdu_check, tbhdu_truth, radius=1):
    """
    """
    
    tbhdu_A_wneib = match_positions(tbhdu_check,tbhdu_truth,radius);

    Ntrue = int(tbhdu_A_wneib.header['ntrue']);
    Nfalse = int(tbhdu_A_wneib.header['nfalse']);
    Ntotal = Ntrue + Nfalse;
    
    # Completeza:
    Comp = float(Ntrue)/Ntotal;
    logging.info("Completeness (N_true/N_total): %s" % (Comp));
    
    # Contaminacao:
    Cont = float(Nfalse)/Ntotal;
    logging.info("Contamination (N_false/N_total): %s" % (Cont));

    return tbhdu_A_wneib;
    
def run_fit(filename_A, filename_B, radius=1):
    """
    Input:
     - filename_A : FITS catalog
     - filename_B : FITS catalog
     - radius : matching distance in pixels
     
    Output:
     - tbHDU
     
    ---
    """
    tbhdu_check = pyfits.open(filename_A)[1];
    tbhdu_truth = pyfits.open(filename_B)[1];
    
    return run(tbhdu_check,tbhdu_truth,radius);

def run_ds9(filename_A, filename_B, radius=1):
    """
    Input:
     - filename_A : DS9 region file
     - filename_B : DS9 region file
     
    Output:
     - tbhdu
    
    ---
    """
    
    # Open two tables, for now, ds9 region files:
    #
    DA_in = ascii_data.read_ds9cat(filename_A);
    DB_in = ascii_data.read_ds9cat(filename_B);
    
    # Translate them to table HDU (FITS tables Header D Unit)
    tbhdu_check = fits_data.dict_to_tbHDU(DA_in,tbname=filename_A);
    tbhdu_truth = fits_data.dict_to_tbHDU(DB_in,tbname=filename_B);
    
    return run(tbhdu_check,tbhdu_truth,radius);

# ==========================

if __name__ == '__main__' :

    import optparse;
    import re;

    usage="\n  %prog [options] <ds9_CheckTable.reg> <ds9_TruthTable.reg>"
    parser = optparse.OptionParser(usage=usage);

    parser.add_option('-r',
                        dest='radius', default=10,
                        help='Distance for the objects matching in both tables [10]');
    parser.add_option('-o',
                        dest='filename', default='matching_points',
                        help='Output tables (DS9 and FITS) filename [matching_points]');
    parser.add_option('--no_DS9reg', action='store_false',
                        dest='write_ds9reg', default=True,
                        help='Avoid DS9 region file output?');
                        
    (opts,args) = parser.parse_args();

    radius = opts.radius;
    outfilename = opts.filename;
    write_ds9reg = opts.write_ds9reg;

    if len(args) < 2 :
        parser.print_help();
        sys.exit(1);


    # Logging handler:
    logfile = outfilename+'.log'
    logging = log.init(logfile,debug=True,verbose=True);


    regionfile_A = args[0];
    regionfile_B = args[1];

    tbhdu_A = run_ds9(regionfile_A, regionfile_B, radius);

    try:
        tbhdu_A.writeto(outfilename+'.fit');
    except IOError:
        os.remove(outfilename+'.fit');
        tbhdu_A.writeto(outfilename+'.fit');

    # Set the dictionary to output to a ds9 region file
    # the points that matched with the truth table will
    # be written as yellow in ds9 region file.

    if write_ds9reg:
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
