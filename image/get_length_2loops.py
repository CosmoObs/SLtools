from __future__ import division
from numpy import argmax, sum

##@package get_length_2loops
# finds the 2 furthest points of an image through the furthest-of-furthest method and computes the length, defined as the sum of the distances from the "center" to each end
#
#@param ximg (x coordinates of the image)
#@param yimg (y coordinates of the image)
#@param cont (counts in each pixel)
#@param source_type
#@return L, max_dist, furthest_pt_1, furthest_pt_2
def get_length_2loops(ximg,yimg,cont, source_type):
	if source_type == 'sersic':
		xcentro = ximg[argmax(cont)]
		ycentro = yimg[argmax(cont)]
	if source_type == 'uniform':
		xcentro = sum(ximg)/len(ximg) 
		ycentro = sum(yimg)/len(yimg)
	ptmaisdist1 = 0
	ptmaisdist2 = 0
	distmax = -1 # vamos encontrar os 2 pts mais distantes entre si em uma imagem
	for i in range(0,len(ximg)): # varre todos os pixels da imagem
		if (ximg[i] - xcentro)**2 + (yimg[i] - ycentro)**2 > distmax:
			distmax = (ximg[i] - xcentro)**2 + (yimg[i] - ycentro)**2
			ptmaisdist1 = i
	distmax = -1 # vamos encontrar os 2 pts mais distantes entre si em uma imagem		
	for j in range(0,len(ximg)): # varre todos os pixels da imagem
		if (ximg[j] - ximg[ptmaisdist1])**2 + (yimg[j] - yimg[ptmaisdist1])**2 > distmax:
			distmax = (ximg[j] - ximg[ptmaisdist1])**2 + (yimg[j] - yimg[ptmaisdist1])**2
			ptmaisdist2 = j
	L1 =  ( (ximg[ptmaisdist1] - xcentro)**2 + (yimg[ptmaisdist1] - ycentro)**2 )**0.5
	L2 =  ( (ximg[ptmaisdist2] - xcentro)**2 + (yimg[ptmaisdist2] - ycentro)**2 )**0.5
	L = L1+ L2
	distmax = (ximg[ptmaisdist1] - ximg[ptmaisdist2])**2 + (yimg[ptmaisdist1] - yimg[ptmaisdist2])**2
	return L, distmax, ptmaisdist1, ptmaisdist2
