#!/usr/bin/env python

"""Function to determine the points near the bisectrix of the line joining two points"""


##@package get_bisectrix_pts
#
# Determines the points of an object closer than (2^0.5)/2 from the bisectrix of the line joining two chosen points of the object
#


#-----------------
import numpy as np


def get_bisectrix_pts( ximg, yimg, ref_position_1, ref_position_2 ):
	"""
	Function to determine the points near the bisectrix of the line joining
	two points.

	Input:
	 - ximg : x coordinates of the image (list)
	 - yimg : y coordinates of the image (list)
	 - ref_position_1 : the (Python) position of one reference point in the 
			    ximg and yimg lists
	 - ref_position_2 : the (Python) position of the other reference point 
			    in the ximg and yimg lists

	Output: 
	 - (list) : x positions of the near bisectrix points
	 - (list) : y positions of the near bisectrix points

	"""


	distmax = ( (ximg[ref_position_1] - ximg[ref_position_2])**2 + (yimg[ref_position_1] - yimg[ref_position_2])**2 )**0.5

	if ximg[ref_position_2] >= ximg[ref_position_1]:
		costheta = (ximg[ref_position_2] - ximg[ref_position_1])/distmax
		sentheta = (yimg[ref_position_2] - yimg[ref_position_1])/distmax

	if ximg[ref_position_1] > ximg[ref_position_2]:
		costheta = (ximg[ref_position_1] - ximg[ref_position_2])/distmax
		sentheta = (yimg[ref_position_1] - yimg[ref_position_2])/distmax
	xl = [] # list that contains the x coordinates on the local system
	yl = [] # list that contains the y coordinates on the local system
	midx = (ximg[ref_position_1] + ximg[ref_position_2])/2.
	midy = (yimg[ref_position_1] + yimg[ref_position_2])/2.
	for i in range (0,len(ximg)):
		xl.append(costheta*(ximg[i] - midx) + sentheta*(yimg[i] - midy))
		yl.append(-sentheta*(ximg[i] - midx) + costheta*(yimg[i] - midy))
	yl = np.array(yl)
	yl = yl - yl[ref_position_1] # so the line that conects the 2 furthest points will be at y=0
	# search for the points closer than (2^0.5)/2 from the y axis
	xlq = []
	ylq = []
	for i in range(0,len(xl)):
		if abs(xl[i]) <=0.70711: # 0.70711 = (2**0.5)/2.		
			xlq.append(ximg[i])
			ylq.append(yimg[i])	
			# xlq.append(xl[i])
			# ylq.append(yl[i])	


	return xlq, ylq

