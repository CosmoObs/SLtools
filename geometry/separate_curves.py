#!/usr/bin/env python
# ==================================
# Authors:
# Angelo Fausti Neto, Pedro Ferreira
# ==================================
  
"""Module to separate curves defined by connected line segments"""


##@package separate_curves
#
# Executable package: NO


from matplotlib.pyplot import *
import numpy as np
from numpy import * 
from random import *
import logging

#=======================================================================================================
def separate_curves(x1, y1, x2, y2, delta=None):

    ''' 
    Separates curves defined by connected line segments.
    
    The arrays (x1,y1) and (x2,y2) define line segments, that when connected, form curves. The 
    requirement is that the input points are initialy sorted. The parameter delta is the criterium to 
    separate one curve from the other and is an upper limit to the size of a line segment. If the 
    distance between two consecutive points (x1[i],y1[i]) and (x1[i+1], y1[i+1]) is larger than delta it
    means we "jumped to another curve". By default, we set delta as 1/10 of a caracteristic curve size 
    (the diagonal of the curve boundary rectangle).

    This function allows cases where only 2 curves are present in the input coordinates.

    Input:
     - x1  <list> : x coordinate of a point
     - y1  <list> : y coordinate of a point
     - x2  <list> : x coordinate of the point consecutive to (x1,y1)
     - y2  <list> : y coordinate of the point consecutive to (x1,y1)

    Output:
     - x1  <list> : x coordinates of points from one curve (in gravlens, is the radial curve)
     - y1  <list> : y coordinates of points from one curve (in gravlens, is the radial curve)
     - x2  <list> : x coordinates of points from the other curve (in gravlens, is the tangential curve)
     - y2  <list> : y coordinates of points from the other curve (in gravlens, is the tangential curve)

    '''

    start=[] # the start point and the end point of the separated curves
    end=[]

    # Default delta

    if delta==None:
        delta=0.1*sqrt(pow((max(x1)-min(x1)),2)+pow((max(y1)-min(y1)),2))

    i=0 

    # if the distance of two consecutive points is larger than delta, we start a new curve

    while i < (len(x1)):

        start.append(i)

        while (i < (len(x1) -1)) and ((pow(x1[i+1]-x1[i],2)+pow(y1[i+1]-y1[i],2)) < (pow(delta,2))): 
            i=i+1

        i=i+1
        end.append(i)

    ncurves = len(start)
    curves = [] # curves[0] is the radial curve and curves[1] is the tangential curve


    for index in range(ncurves):
        curves.append([x1[start[index]:end[index]], y1[start[index]:end[index]], x2[start[index]:end[index]], y2[start[index]:end[index]]])

    if len(curves) == 1:
        logging.warning('Only one caustic (and CC) was found (usually its the tangential). Maybe you are approaching gravlens precision. Try changing units (ex., from arcsec to miliarcsec).')


    radial_curve = curves[0] # [ [x1_i], [y1_i], [x2_i], [y2_i] ]
    tang_curve = curves[1] # [ [x1_i], [y1_i], [x2_i], [y2_i] ]
    # connecting each (x1_i,y1_i) to (x2_i,y2_i) will draw a closed curve, but the (x1_i,y1_i)
    # already contains all the curve points

    x_rad, y_rad, x_tg, y_tg = radial_curve[0], radial_curve[1], tang_curve[0], tang_curve[1]

    # repeating the 1st element in the end of each array will make the plot easier
    # first, convert the array to a list
    x_rad, y_rad, x_tg, y_tg = list(x_rad), list(y_rad), list(x_tg), list(y_tg)

    x_rad.append(x_rad[0])
    y_rad.append(y_rad[0])
    x_tg.append(x_tg[0])
    y_tg.append(y_tg[0])

    x_rad, y_rad, x_tg, y_tg = np.array(x_rad), np.array(y_rad), np.array(x_tg), np.array(y_tg)

    return x_rad, y_rad, x_tg, y_tg



