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


# given 3 points on a plane, returns circunference radius and center
# see http://en.wikipedia.org/wiki/Circumcircle
def circ_3pts(p1, p2, p3):

	x1, y1 = p1
	x2, y2 = p2
	x3, y3 = p3

	# definitions
	P1_P2 = ( (x1 - x2)**2 + (y1 - y2)**2 )**0.5
	P2_P3 = ( (x2 - x3)**2 + (y2 - y3)**2 )**0.5
	P3_P1 = ( (x3 - x1)**2 + (y3 - y1)**2 )**0.5
	P1_P2_vec_P2_P3 = abs( (x1 - x2)*(y2 - y3) - (x2 - x3)*(y1 - y2) )
	# ==================================================================
	# radius

	R = (P1_P2 * P2_P3 * P3_P1)/( 2*P1_P2_vec_P2_P3 )

	# definitions
	alpha = ( P2_P3**2 * ( ( (x1 - x2)* (x1 - x3) ) + ( (y1 - y2)* (y1 - y3) ) ) )/ (2*P1_P2_vec_P2_P3**2)
	beta  = ( P3_P1**2 * ( ( (x2 - x1)* (x2 - x3) ) + ( (y2 - y1)* (y2 - y3) ) ) )/ (2*P1_P2_vec_P2_P3**2)
	gamma = ( P1_P2**2 * ( ( (x3 - x1)* (x3 - x2) ) + ( (y3 - y1)* (y3 - y2) ) ) )/ (2*P1_P2_vec_P2_P3**2)
	# ==================================================================
	# center (x0, y0)

	x0 = alpha*x1 + beta*x2 + gamma*x3
	y0 = alpha*y1 + beta*y2 + gamma*y3


	return R, x0, y0


## length_arc_circ
# given the properties of the object, this function decides if the object is an arc or not
#
#@param  L1 (greater distance between 2 points)
#@param R (radius of the circunference)
#@return length of the arc of cincunference that crosses the 2 furthest points of the arc plus a 3rd one (defined implicitly on R)
def get_length_arc_circ(L1,R):
	L3 = np.arccos(1 - (L1**2)/(2*(R**2)) )*R # this formula is directly shown with the "cosine law" (lei dos cossenos)
	return L3

