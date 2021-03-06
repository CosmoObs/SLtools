#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =================================================
# Authors:
# Angelo Fausti Neto, Pedro Ferreira, Carlos Brandt
# =================================================
  
"""Module to identify multiple curves among ordered points (line segments)"""


##@package separate_curves
#
#
# Module to separate curves defined by connected line segments.
# Its single function, separate_curves, receives as inputs sets of consecutive points (see on line 
# documentation), forming line segments. The code tries to find where the line segments joins, closing 
# the curve and possibly starting another one.
# Executable package: NO


from matplotlib.pyplot import *
from numpy import * 
from random import *
import numpy as np

#=======================================================================================================
def separate_curves(x1, y1, x2, y2, delta=None):

    ''' 
    Separates curves defined by connected line segments.
    
    The arrays (x1,y1) and (x2,y2) define line segments that when connected form curves. The 
    requirement is that the input points are initialy sorted. The parameter delta is the criterium to 
    separate one curve from the other and is an upper limit to the size of a line segment. If the 
    distance between two consecutive points (x1[i],y1[i]) and (x1[i+1], y1[i+1]) is larger than delta it
    means we "jumped to another curve". By default, we set delta as 1/10 of a caracteristic curve size 
    (the diagonal of the curve boundary rectangle).

    Input:
     - x1  <list> : x coordinate of a point
     - y1  <list> : y coordinate of a point
     - x2  <list> : x coordinate of the point consecutive to (x1,y1)
     - y2  <list> : y coordinate of the point consecutive to (x1,y1)

    Output:
     - curves  <list> : list of curves. Each element curves[i] is a closed curve with 4 arrays. The 
                        first 2 arrays being the (x,y) coordinates of the curve points and the 3rd and 
			4th collumns are the consecutive points to the 1st and 2nd collumns
     - start   <list> : indexes (integers) of x1, y1, x2, y2 indicating when a curve starts
     - end     <list> : indexes (integers) of x1, y1, x2, y2 indicating when a curve ends

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
    curves = [] # in the SL case, curves[0] is the radial curve and curves[1] is the tangential curve


    for index in range(ncurves):
        curves.append([x1[start[index]:end[index]], y1[start[index]:end[index]], x2[start[index]:end[index]], y2[start[index]:end[index]]])


    return curves, start, end


