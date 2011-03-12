#!/usr/bin/env python
# ==================================
# Authors:
# Pedro Ferreira - pferreira@dfte.ufrn.br
# ==================================

"""Function to determine the two extrema of a set of points through a reference point"""


##@package get_extrema_2loops
#
# Find the 2 furthest points of an image. It uses a reference point to find the first extreme and then 
# finds the furthest point of this extreme through the furthest-of-furthest method.


from numpy import argmax, sum
import numpy as np;


# ======================================================================================================
def get_extrema_2loops( ximg, yimg, ref_position ):
	"""
	Function to determine the two extrema of a set of points through a reference point. 

	The first extreme point (P1) is found by searching for the point furthest from the reference 
	point (usually is some definition of center). The second extreme point is the furthest one from 
	point P1.


	Input:
	 - ximg        <list> : x coordinates of the image
	 - yimg        <list> : y coordinates of the image
	 - ref_position <int> : the (Python) position of the reference point in the ximg and yimg lists

	Output: 
	 - <int> : position in the lists ximg and yimg of one extreme point
	 - <int> : position in the lists ximg and yimg of the other extreme point

	"""


	# define the reference point coordinates
	x_ref = ximg[ref_position] 
	y_ref = yimg[ref_position] 

	# find the furthest point from the reference point (1st extreme)
	furthest_pt1 = 0
	distmax = -1
	for i in range(0,len(ximg)): # loops over all points
		if (ximg[i] - x_ref)**2 + (yimg[i] - y_ref)**2 > distmax:
			distmax = (ximg[i] - x_ref)**2 + (yimg[i] - y_ref)**2
			furthest_pt1 = i

	# find the furthest point from the first extreme (2nd extreme)
	furthest_pt2 = 0
	distmax = -1	
	for j in range(0,len(ximg)): # loops over all points
		if (ximg[j] - ximg[furthest_pt1])**2 + (yimg[j] - yimg[furthest_pt1])**2 > distmax:
			distmax = (ximg[j] - ximg[furthest_pt1])**2 + (yimg[j] - yimg[furthest_pt1])**2
			furthest_pt2 = j


	return furthest_pt1,furthest_pt2;

