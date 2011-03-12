#!/usr/bin/env python
# ==================================
# Authors:
# Pedro Ferreira - pferreira@dfte.ufrn.br
# ==================================

""" Package to determine the points near the bisectrix of the line joining two points"""


##@package get_bisectrix_pts
#
#
# Determines the points of an object closer than \f$ \sqrt{2}/2 \f$ from the bisectrix of the line 
# joining two chosen points of the object. 
#



import numpy as np


# ======================================================================================================
def get_bisectrix_pts(ximg, yimg, ref_position_1, ref_position_2):
	"""
	Function to determine the points near the bisectrix of the line joining two points.

	Given a set of an object coordinates (ximg,yimg) and 2 reference points, this function 
	calculates all points that are closer than \f$ \sqrt{2}/2 \f$ from the bisectrix line (line 
	perpendicular to the line joinning the 2 input reference points).

	Input:
	 - ximg          <list> : list with x coordinates of the image
	 - yimg          <list> : list with y coordinates of the image
	 - ref_position_1 <int> : the (Python) position of one reference point in the ximg and yimg lists
	 - ref_position_2 <int> : the (Python) position of the other reference point in the ximg and yimg 
			    lists

	Output: 
	 - <list> : x positions of the near bisectrix points
	 - <list> : y positions of the near bisectrix points

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



# ======================================================================================================
def circ_3pts(p1, p2, p3):
	"""
	Determines the circunference that passes through 3 points.

	Given 3 points (x,y) on a plane, returns circunference radius and center. The method implemented
	in this funcion is the one described in <A HREF="http://en.wikipedia.org/wiki/Circumcircle"> 
	CircumcircleWikipedia</A> .

	Input:
	 - p1 <list> : pair (x,y) of coordinates of the 1st point
	 - p2 <list> : pair (x,y) of coordinates of the 2nd point
	 - p3 <list> : pair (x,y) of coordinates of the 3rd point

	Output: 
	 - <float> : circunference radius
	 - <float> : x position of the circunference center
	 - <float> : y position of the circunference center

	"""

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



# ======================================================================================================
def get_length_arc_circ(L1,R):
	"""
	Calculates the length of an arc of a circle.

	Given a circle of radius R, calculates the length of the arc of circle between 2
	points L1 apart from each other.
	The length L of the arc is calculated through the cosine law:
	\f$ L = R*arccos(1 - (L1^2)/(2R^2)) \f$

	Input:
	 - L1 <float> : shortest distance between the 2 extrema points of the desired arc of circle
	 - R  <float> : radius of the circle

	Output: 
	 - <float> : length of the arc of circle

	"""

	L3 = np.arccos(1 - (L1**2)/(2*(R**2)) )*R # this formula is directly shown with the "cosine law" (lei dos cossenos)
	return L3



