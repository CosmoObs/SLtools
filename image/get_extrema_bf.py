from __future__ import division
from numpy import argmax, sum


def get_extrema_bf(ximg, yimg):

	furthest_pt1 = 0
	furthest_pt2 = 0
	distmax = -1 # lets find the 2 furthest points in an image
	for i in range(0,len(ximg)): # goes through all image pixels
		for j in range(i+1,len(ximg)):
			if (ximg[i] - ximg[j])**2 + (yimg[i] - yimg[j])**2 > distmax:
				distmax = (ximg[i] - ximg[j])**2 + (yimg[i] - yimg[j])**2
				furthest_pt1 = i
				furthest_pt2 = j	

	return furthest_pt1, furthest_pt2
