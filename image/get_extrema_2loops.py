from __future__ import division
from numpy import argmax, sum

##@package get_extrema_2loops
# finds the 2 furthest points of an image through the furthest-of-furthest method (using input reference point)
#
#@param ximg (x coordinates of the image)
#@param yimg (y coordinates of the image)
#@param ref_coordinates (a list: [x,y])
#@return L, max_dist, furthest_pt_1, furthest_pt_2
def get_extrema_2loops(ximg, yimg, ref_coordinates):

	x_ref = ref_coordinates[0] 
	y_ref = ref_coordinates[1] 

	furthest_pt1 = 0
	furthest_pt2 = 0
	distmax = -1 # vamos encontrar os 2 pts mais distantes entre si em uma imagem
	for i in range(0,len(ximg)): # varre todos os pixels da imagem
		if (ximg[i] - x_ref)**2 + (yimg[i] - y_ref)**2 > distmax:
			distmax = (ximg[i] - x_ref)**2 + (yimg[i] - y_ref)**2
			furthest_pt1 = i
	distmax = -1 # vamos encontrar os 2 pts mais distantes entre si em uma imagem		
	for j in range(0,len(ximg)): # varre todos os pixels da imagem
		if (ximg[j] - ximg[furthest_pt1])**2 + (yimg[j] - yimg[furthest_pt1])**2 > distmax:
			distmax = (ximg[j] - ximg[furthest_pt1])**2 + (yimg[j] - yimg[furthest_pt1])**2
			furthest_pt2 = j
#	L1 =  ( (ximg[furthest_pt1] - x_ref)**2 + (yimg[furthest_pt1] - y_ref)**2 )**0.5
#	L2 =  ( (ximg[furthest_pt2] - x_ref)**2 + (yimg[furthest_pt2] - y_ref)**2 )**0.5
#	L = L1+ L2
#	distmax = (ximg[furthest_pt1] - ximg[furthest_pt2])**2 + (yimg[furthest_pt1] - yimg[furthest_pt2])**2
	return furthest_pt1, furthest_pt2


