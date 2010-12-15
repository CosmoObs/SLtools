#!/usr/bin/python

""" Module to fit arcs parameters using ArcEllipse analytical expression."""

##@module fit_arcellipse

"""
 This module fits arcs parameters using ArcEllipse analytical expression combined with a Sersic profile to describe the radial 
 light distribution of the arc. The ArcEllipse expression is the result of a deformation of an ellipse such 
 that one of its main axes becomes a circle segment. 

 
 It uses chi2 minimization to find the best fit parameters values. The chi2 minimization is done using Pyminuit 
 (see http://code.google.com/p/pyminuit/).

"""

from __future__ import division
import numpy
import math
import pyfits
import sltools
from functions import *
import minuit
from sltools.image import get_header_parameter


#========================================================================
def f(x,y,x0,y0,rc,ellip,I0,re,n,theta0):
    """ Function to compute the value of the surface brightness of the
	(x,y) pixel following the ArcEllipse prescription.

	Input:
	- x: x position (pixel)
	- y: y position (pixel)
	- x0: x position of the cluster center or the point where the arc is centered (pixel)
	- y0: y position of the cluster center or the point where the arc is centered (pixel)
	- r_c: distance of the arc center related to the cluster center (pixel)
	- ellip: ellipticity of the arc
	- I_0: central surface brightness of the Sersic profile
	- r_e: radius that encloses half of the total luminosity of the Sersic profile (pixel)
	- n: index that describes the shape of the Sersic profile 
	- theta0: orientation of the arc center related to the cluster center or the point where the arc is centered (radians)

	Output:
	- the value of the surface brightness of the (x,y) point following the ArcEllipse prescription
	"""
	
    x = float(x)
	y = float(y)	
	#x_0 = 10.
	#y_0 = 10.
	#ellip = 0.857142857143
	#n = 1.0
	#r_e = 2.22222222222
	#r_c = 74.0740740741
	#theta_0 = 0.785398163397
	#I_0 = 32.974346 

	theta = math.atan2((y - y0),(x - x0))
		
	r = math.sqrt((x - x0)**2 + (y - y0)**2)
	theta_diff = abs(theta - theta0)
	if (theta_diff >= (math.pi)):
		theta_diff = 2. * math.pi - theta_diff

	r_prof = math.sqrt((rc * theta_diff * (1. - ellip))**2 + (r - rc)**2 )
        bn = 1.9992*float(n) - 0.3271
	exp_arg =  bn * ((r_prof / re) ** (1. / float(n)) )
	#print rc,theta_diff,ellip,r,r_prof,re
	return I0 * math.exp(-1. * exp_arg)	


#========================================================================
def chi2(x0,y0,rc,ellip,I0,re,n,theta0):
    """ Function to compute the chi2 of the f(x,y,x0,y0,rc,ellip,I0,re,n,theta0) function related to the data.
    
	Since we do not have an error associated with the data, we use Poissonic error (sigma[i][j] = sqrt(data[i][j]) )

	Input:
	- x0: x position of the cluster center or the point where the arc is centered (pixel)
	- y0: y position of the cluster center or the point where the arc is centered (pixel)
	- r_c: distance of the arc center related to the cluster center (pixel)
	- ellip: ellipticity of the arc
	- I_0: central surface brightness of the Sersic profile
	- r_e: radius that encloses half of the total luminosity of the Sersic profile (pixel)
	- n: index that describes the shape of the Sersic profile 
	- theta0: orientation of the arc center related to the cluster center or the point where the arc is centered (radians)

	Output:
	- the value of the chi2
	"""

	c2 = 0.

	for i in range(dim_x):
		for j in range(dim_y):
			c2 += ((f(i,j,x0,y0,rc,ellip,I0,re,n,theta0)-data[j][i])**2)/data[j][i]
	return c2


###########################################################################
def run():
    #Get image data
    data = pyfits.getdata("arc.fits")
    #Get image parameters
    dim_x,dim_y = get_header_parameter("arc.fits",'NAXIS1','NAXIS2') 

    #Using Pyminuit

    m = minuit.Minuit(chi2)

    #Initial guesses
    m.values["x0"] = 10.
    m.values["y0"] = 10.
    m.values["rc"] = 74.5
    m.values["ellip"] = 0.85
    #m.values["I0"] = 28.
    m.values["re"] = 2.22
    m.values["n"] = 1.0
    m.values["theta0"] = 0.78

    #m.limits["I0"] = (10., 50.)

    m.printMode = 1

    #print m.limits

    #Minimization
    m.migrad() 
    #m.simplex()
    #m.hesse()
    m.scan(("I0",10,0,50))


    print m.values

if __name__ == '__main__':
    run();
