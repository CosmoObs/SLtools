from __future__ import division
from numpy import argmax, sum


def get_length_brute(ximg,yimg,cont, source_type):
	if source_type == 'sersic':
		xcentro = ximg[argmax(cont)]
		ycentro = yimg[argmax(cont)]
	if source_type == 'uniform':
		xcentro = sum(ximg)/len(ximg) 
		ycentro = sum(yimg)/len(yimg)
	ptmaisdist1 = 0
	ptmaisdist2 = 0
	distmax = -1 # lets find the 2 furthest points in an image
	for i in range(0,len(ximg)): # goes through all image pixels
		for j in range(i+1,len(ximg)):
			if (ximg[i] - ximg[j])**2 + (yimg[i] - yimg[j])**2 > distmax:
				distmax = (ximg[i] - ximg[j])**2 + (yimg[i] - yimg[j])**2
				ptmaisdist1 = i
				ptmaisdist2 = j	
	#print "Os pts mais dist entre si sao (%d,%d) e (%d,%d) "  % (ximg[ptmaisdist1],yimg[ptmaisdist1], ximg[ptmaisdist2],yimg[ptmaisdist2]) 
	L1 =  ( (ximg[ptmaisdist1] - xcentro)**2 + (yimg[ptmaisdist1] - ycentro)**2 )**0.5
	L2 =  ( (ximg[ptmaisdist2] - xcentro)**2 + (yimg[ptmaisdist2] - ycentro)**2 )**0.5
	L = L1 + L2 # distance from the central point till a extreme + dist from the central point till the other extreme
	return L, distmax, ptmaisdist1, ptmaisdist2
