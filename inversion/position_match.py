#!/usr/bin/env python
# ==================================
# Author: Gabriel Bartosch Caminha - gbcaminha@gmail.com
# ==================================

"""
FIXME - package to make positions match
"""
import numpy as np
import math

import functions as fc

def max_diff(ra_ref, dec_ref, ra_vec, dec_vec):
    """
    Function to compute the maximun distance between a point and a set of points
    Input
     - ra_ref  float : 
     - def_ref  float : 
     - ra_vec  float : 
     - dec_vec  float : 
    Output
     -
    #FIXME - finish documentation
    """
    if (type(ra_vec) != np.array):
        ra_vec = np.array(ra_vec)
    if (type(dec_vec) != np.array):
       dec_vec = np.array(dec_vec)

    diff_pos = []

    for i in range(len(ra_vec)):
        ra_diff = (ra_ref - ra_vec[i]) * math.cos(dec_ref - dec_vec[i])
        dist = fc.modulus( ra_diff, (dec_ref - dec_vec[i]) )
        diff_pos.append( math.fabs(dist) )

    max_index = diff_pos.index( max(diff_pos) )
    max_index1 = np.where( ra_vec == ra_vec[max_index] )

    return max(diff_pos)*60.0, max_index1[0][0]

def min_diff_mask(x_ref, y_ref, x_vec, y_vec, tolerance = 1.5):
    """
    Optimized functions to compute the mininum distance between a point and a array
    Input
     - x_ref  float : 
     - y_ref  float : 
     - x_vec  float : 
     - y_vec  float : 
     - toll   float : tolerance for the mask, in arc minutes 
    Output
     -
    #FIXME - finish documentation
    """
    x_min = x_ref - tolerance/60.0
    x_max = x_ref + tolerance/60.0
    y_min = y_ref - tolerance/60.0
    y_max = y_ref + tolerance/60.0
    if (type(x_vec) != np.array):
        x_vec = np.array(x_vec)
    if (type(y_vec) != np.array):
       y_vec = np.array(y_vec)

    msk = (x_vec >= x_min) & (x_vec <= x_max) & \
          (y_vec >= y_min) & (y_vec <= y_max)

    x_vec_msk = x_vec[msk]
    y_vec_msk = y_vec[msk]
    #print len(x_vec_msk), len(y_vec_msk)

    if len(x_vec_msk) == 0:
        return 1000, -1

    diff_pos = []
    for i in range(len(x_vec_msk)):
        ra_diff = (x_ref - x_vec_msk[i]) * math.cos(y_ref - y_vec_msk[i])
        dist = fc.modulus( ra_diff, (y_ref - y_vec_msk[i]) )
        diff_pos.append( math.fabs(dist) )

    min_index = diff_pos.index( min(diff_pos) )
    min_index1 = np.where( x_vec == x_vec_msk[min_index] )
    #min_index2 = np.where( y_vec == y_vec_msk[min_index] )

    #print min_index1[0][0], min_index2[0][0]

    #print min(diff_pos)*60.0, min_index1[0][0]
    return min(diff_pos)*60.0, min_index1[0][0]



def main_position_match():
    """
    Main function
    Input
     -
    Output
     -
    #FIXME - finish documentation
    """
    print "Hello world! from main_position_match()"

if __name__ == '__main__':

    print "Hello world!"
    main_position_match()
