#!/usr/bin/env python
# ==================================
# Authors:
# Pedro Ferreira - pferreira@dfte.ufrn.br
# ==================================

"""Function to determine (exactly) the two extrema of a set of points"""


##@package get_extrema_bf
#
# Finds the 2 furthest points of an image comparing the distances between each pair of points ("brute 
# force" method).


from numpy import argmax, sum


# ======================================================================================================
def get_extrema_bf(ximg, yimg):

	"""
	Function to determine the two extrema of a set of points through calculating the distance 
	between every pair of points.

	Input:
	 - ximg <list> : x coordinates of the image
	 - yimg <list> : y coordinates of the image

	Output: 
	 - <int> : position in the lists ximg and yimg of one extreme point
	 - <int> : position in the lists ximg and yimg of the other extreme point

	"""


	# calculate the distance of each pair of points recording the greatest 
	# distance and point positions
	furthest_pt1 = 0
	furthest_pt2 = 0
	distmax = -1
	for i in range(0,len(ximg)):
		for j in range(i+1,len(ximg)):
			if (ximg[i] - ximg[j])**2 + (yimg[i] - yimg[j])**2 > distmax:
				distmax = (ximg[i] - ximg[j])**2 + (yimg[i] - yimg[j])**2
				furthest_pt1 = i
				furthest_pt2 = j	


	return furthest_pt1, furthest_pt2
