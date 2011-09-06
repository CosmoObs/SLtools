#!/usr/bin/env python

"""Module to do tables matching, functions to deal with it"""

##@package catalog
##@file table_matching

import sys
import os
import logging

import numpy as np
import pyfits

from sltools.catalog import fits_data,ascii_data
from sltools.io import log


# ---
def nearest_neighbour(xy_check, xy_truth):
    """
    Function to compute and return the set of nearest point
        
    This function computes for each entry of 'xy_check', the
    nearest point in 'xy_truth'.
    
    Is returned, respective to each entry, the index and the
    distance measured correspondig to each identified nearest
    point.
    
    Input:
     - xy_check  [(int,int),...] : List of (x,y) positions for matching against "B"
     - xy_truth  [(int,int),...] : List of (x,y) positions to check "A" for
    
    Output:
     - near_to  [(int,float),...] : List of (index,distance) of nearby object from "A"


    * hint: [(a,b),...] is the output of zip(A,B) command

    ---
    """
    
    x_check,y_check = zip(*xy_check)
    x_truth,y_truth = zip(*xy_truth)
    logging.debug("X,Y of Check table: \n(X):%s \n(Y):%s",x_check,y_check)
    logging.debug("X,Y of Truth table: \n(X):%s \n(Y):%s",x_truth,y_truth)

    length_check = len(xy_check)
    length_truth = len(xy_truth)    
    logging.debug("Size of given point lists: %d(check), %d(truth)",length_check,length_truth)

    Mat = np.zeros((length_check,length_truth))
    
    for i in xrange(length_check):
        Mat[i] = (x_check[i] - x_truth[:])**2 + (y_check[i] - y_truth[:])**2
    
    dist_min_truth2check = np.sqrt(np.amin(Mat,axis=1))
    indx_min_truth2check = np.argmin(Mat,axis=1)
    logging.info("Nearst neibour distance: %s",dist_min_truth2check)
    logging.info("Nearst neibour index: %s",indx_min_truth2check)
    
    ind_dist_truth2check = zip(indx_min_truth2check,dist_min_truth2check)
    logging.debug("Index and Distante of nearest neibours: %s",ind_dist_truth2check)
    
    near_to = ind_dist_truth2check
    return (near_to)

# ---
def match_positions(xy_check, xy_truth, radius):
    """
    Check if positions in tables match within a given radius
    
    Required fieldnames inside tbhdu's:
     'X' and 'Y' or 'X_IMAGE' and 'Y_IMAGE'
    
    Input:
     - xy_check  [(int,int),...] : Being-checked centroids
     - xy_truth  [(int,int),...] : Truth positions
     - radius              float : Distance parameters for objects position matching
    
    Output:
     - dict_check {} : New table with matching information (fields):
            'X'        = [int] x position of input table
            'Y'        = [int] y position of input table
            'XY_MATCH' = [bool] was X,Y object matched within 'radius'?
            'X_NEAR'   = [int] x position of nearest object
            'Y_NEAR'   = [int] y position of nearest object


    * hint: [(a,b),...] is the output of zip(A,B) command

    ---
    """
    
    x_check,y_check = zip(*xy_check)
    x_truth,y_truth = zip(*xy_truth)
    logging.debug("X,Y of Check table: \n(X):%s \n(Y):%s",x_check,y_check)
    logging.debug("X,Y of Truth table: \n(X):%s \n(Y):%s",x_truth,y_truth)

    len_check = len(xy_check)
    len_truth = len(xy_truth)
    logging.debug("Size of given point lists: %d(check), %d(truth)",len_check,len_truth)
    
    Check_near_neibors = nearest_neighbour(xy_check,xy_truth)
    Check_nn_indx,Check_nn_dist = zip(*Check_near_neibors)
    logging.debug("Check-table (from Truth) nearest neighbours (index,distance): %s",Check_near_neibors)
    
    logging.info("Matching radius: %s", radius)
    Check_match_neibors = [ float(_d) < float(radius) for _d in Check_nn_dist ]
    logging.info("Check-table matched objects: %s", Check_match_neibors)
    
    Check_neibor_truthX = [ x_truth[i] for i in Check_nn_indx ]
    Check_neibor_truthY = [ y_truth[i] for i in Check_nn_indx ]
    logging.info("Position of nearest objects: %s,%s",Check_neibor_truthX,Check_neibor_truthY)

    # Alias:
    naa = np.asarray

    dict_check = {'X_NEAR':naa(Check_neibor_truthX),'Y_NEAR':naa(Check_neibor_truthY)}
    dict_check.update({'XY_MATCH':naa(Check_match_neibors)})
    dict_check.update({'X':naa(x_check),'Y':naa(y_check)})
    logging.info("Matching infos: \n(Did match?): %s \n(X): %s \n(Y): %s \n(X_NEAR): %s \n(Y_NEAR): %s",
        dict_check['XY_MATCH'],dict_check['X'],dict_check['Y'],dict_check['X_NEAR'],dict_check['Y_NEAR'])
    
    return (dict_check);

# ---
def run(tbhdu_check, tbhdu_truth, radius=1):
    """
    Check if positions in tables match within a given radius
    
    Required fieldnames inside tbhdu's:
     'X' and 'Y' or 'X_IMAGE' and 'Y_IMAGE'
    
    Input:
     - tbhdu_check  tbHDU : Being-checked Table (FITS)
     - tbhdu_truth  tbHDU : Truth Table (FITS)
     - radius       float : Distance parameters for objects position matching
    
    Output:
     - new_tbhdu tbHDU : New table with matching information (fields):
            'X'        = [int] x position of input table
            'Y'        = [int] y position of input table
            'XY_MATCH' = [bool] was X,Y object matched within 'radius'?
            'X_NEAR'   = [int] x position of nearest object
            'Y_NEAR'   = [int] y position of nearest object


    * obs: tbHDU = pyfits.BinTableHDU
    
    ---
    """
    
    Ntotal_truth = tbhdu_truth.header['naxis2']
    Ntotal_check = tbhdu_check.header['naxis2']
    logging.debug("Input params (tbC,tbT,radius): %s,%s,%s" % (tbhdu_check.name,tbhdu_truth.name,radius))
    logging.debug("Check table length: %d",Ntotal_check)
    logging.debug("Truth table length: %d",Ntotal_truth)

    tbname_check = tbhdu_check.name
    tbname_truth = tbhdu_truth.name
    tbname_out = tbname_check+"_matched_objs"

    try:
        x_A = tbhdu_check.data.field('X')
        y_A = tbhdu_check.data.field('Y')
    except:
        x_A = tbhdu_check.data.field('X_IMAGE')
        y_A = tbhdu_check.data.field('Y_IMAGE')
    try:
        x_B = tbhdu_truth.data.field('X')
        y_B = tbhdu_truth.data.field('Y')
    except:
        x_B = tbhdu_truth.data.field('X_IMAGE')
        y_B = tbhdu_truth.data.field('Y_IMAGE')
    
    xy_check = zip(x_A,y_A)
    xy_truth = zip(x_B,y_B)
    
    dict_Matchout = match_positions(xy_check,xy_truth,radius)
    
#    Ntrue = int(tbhdu_check_wneib.header['ntrue']);
#    Nfalse = int(tbhdu_check_wneib.header['nfalse']);

    new_tbhdu = fits_data.dict_to_tbHDU(dict_Matchout,tbname_out)
    new_tbhdu.header.update('hierarch truthtbname',tbname_truth,"Name of Truth table used")
    new_tbhdu.header.update('hierarch ntottruth',Ntotal_truth,"Total number of arcs in Truth table")
    new_tbhdu.header.update('ntrue',Ntrue,"Number of found objs in Check table")
    new_tbhdu.header.update('nfalse',Nfalse,"Number of not-found objs in Check table")
    Ntrue = np.where(new_tbhdu.data.field('XY_MATCH') == True).size
    Nfalse = np.where(new_tbhdu.data.field('XY_MATCH') == False).size
    Ntotal = Ntrue + Nfalse
    logging.debug("Matched table columns/data: %s",new_tbhdu.data)
    
    # Completeza:
    Comp = float(Ntrue)/Ntotal
    logging.info("Completeness (N_true/N_total): %s" % (Comp))
    
    # Contaminacao:
    Cont = float(Nfalse)/Ntotal
    logging.info("Contamination (N_false/N_total): %s" % (Cont))

    return new_tbhdu

# ---

# =====
# \cond
# =====

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

"""
# ==========================
# Shell interface
# 
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
        sys.exit(2);


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

"""
# ========
# \endcond
# ========
