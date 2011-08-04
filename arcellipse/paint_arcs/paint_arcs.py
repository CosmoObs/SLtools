#!/usr/bin/python
# =====================================================
# Authors:
# Cristina Furlanetto - furlanetto.cristina@gmail.com
# =====================================================


""" Module to draw geometric figures that mimic arcs and add them to images using the prescription of ArcEllipse. """

##@module paint_arcs
#
# This module draws geometric figures that mimic arcs using the prescription of
# ArcEllipse, which consists in a deformation of an ellipse such that one of its
# main axes becomes a circle segment.
#
# It uses Sersic profiles as the radial light distribution of the arcs. It also 
# performs Gaussian convolution and adds Poissonic noise to arc images.
# 
# The main goal of this code is that it allows the user to control all the input
# parameters of the arcs that are going to be added. 
#
#
import re
import sys
import os
import numpy
import math
import pyfits
import sltools
import string
from sltools.image import gauss_convolution
from sltools.image import add_noise
#import add_noise
from sltools.image import get_header_parameter
from sltools.image import imcp
from sltools.io import config_parser as cp
from sltools.catalog import open_fits_catalog
from sltools.catalog.halos import *
#from sltools.image import get_image_limits
from scipy.integrate import dblquad, quad

sys.path.append( os.path.abspath( '.' ))
cwd = os.getcwd()+'/'

ps_factor = 5.
pix_scale = 0.27

def convert_arcsec_pix(x,pix_scale):
	'''Convert the units from arcsec to pixels.
	@param x value in arcsec
	@param pix_scale pixel scale (arcsec/pix)
	@return value in pixels
	'''
	return x / pix_scale

def convert_pix_arcsec(x,pix_scale):
	'''Convert the units from pixels to arcsec.
	@param x value in pixel
	@param pix_scale pixel scale (arcsec/pix)
	@return value in arcsec
	'''
	return x * pix_scale

def get_b_n(m):

	''' Finds the b_n coefficient of the Sersic profile.
	This expression is more accurate thant the usual b_n = 1.9992*n - 0.3271
	Reference: Ciotti, L., Bertin, G. 1999, A&A, 352, 447
	MacArthur, L.A., Courteau, S., & Holtzman, J.A. 2003, ApJ, 582, 689
	b_n is chosen so that half the total luminosity predicted by the Sersic profile comes from r <= r_e .

	@param m Sersic index 

	@return b_n coefficient of the Sersic profile 	
	'''

	if m == 0.0: return -0.3271
	b_n = 2.0*m - 1.0/3.0 + 4.0/m/405.0 + 46.0/m/m/25515.0 \
		+ 131.0/m/m/m/1148175.0  \
                - 2194697.0/m/m/m/m/30690717750.0

	return b_n


def arcellipse_sbmap_r_theta(r,theta_diff,n,b_n,r_c,r_e,I_0,ellip):

	""" 
	Function to compute the value of the surface brightness of the
	(r,theta) point following the ArcEllipse prescription.

	Input:
	- r: radius (pixel)
	- theta_diff: difference between the position angle (theta) and the arc center position angle (theta_0) (radians)
		the angles are measured conterclockwise starting in the x-axis
	- n: index that describes the slope of the Sersic profile 
	- b_n: coefficient of the Sersic profile that depends on the n index
		b_n is chosen so that half the total luminosity predicted by the Sersic profile comes from r <= r_e 
	- r_c: distance of the arc center to the point where the arc is centered (pixel)
	- r_e: radius that encloses half of the total luminosity of the Sersic profile (pixel)
	- I_0: central surface brightness of the Sersic profile
	- ellip: ellipticity of the arc

	Output:
	- the value of the surface brightness of the (r,theta) point following the ArcEllipse prescription
	"""

	r_prof = math.sqrt((r_c * theta_diff * (1. - ellip))**2 + (r - r_c)**2 )

	exp_arg = - b_n * (r_prof / r_e) ** (1. / float(n)) 

	return I_0 * math.exp(exp_arg)


def arcellipse_sbmap_x_y(x,y,x_0,y_0,theta_0,n,b_n,r_c,r_e,I_0,ellip):

	""" 
	Function to compute the value of the surface brightness of the
	(x,y) pixel following the ArcEllipse prescription.

	Input:
	- x: x position (pixel)
	- y: y position (pixel)
	- x_0: x position of the point where the arc is centered (pixel)
	- y_0: y position of the point where the arc is centered (pixel)
	- theta_0: orientation of the arc center related to the x-axis (radians)
	- n: index that describes the slope of the Sersic profile 
	- b_n: coefficient of the Sersic profile that depends on the n index
		b_n is chosen so that half the total luminosity predicted by the Sersic profile comes from r <= r_e 
	- r_c: distance of the arc center related to the point where the arc is centered (pixel)
	- r_e: radius that encloses half of the total luminosity of the Sersic profile (pixel)
	- I_0: central surface brightness of the Sersic profile
	- ellip: ellipticity of the arc

	Output:
	- the value of the surface brightness of the (x,y) point following the ArcEllipse prescription
	"""

	theta = math.atan2((y - y_0),(x - x_0))
		
	r = math.sqrt((x - x_0)**2 + (y - y_0)**2)

	intens = 0.


       	theta_diff = abs(theta - theta_0)
       	if (theta_diff >= (math.pi)):
	       	theta_diff = 2. * math.pi - theta_diff




	r_prof = math.sqrt((r_c * theta_diff * (1. - ellip))**2 + (r - r_c)**2 )

	exp_arg = - b_n * (r_prof / r_e) ** (1. / float(n))

	intens =  I_0 * math.exp(exp_arg)

		

	return intens


def derivative(f):
	""" 
	Function to compute the numerical derivative of the arcellipse_sbmap_r_theta function at (r,theta). 

	Input:
	- r: radius (pixel)
	- theta_diff: difference between the position angle (theta) and the arc center position angle (theta_0) - (radians)
	- n: index that describes the slope of the Sersic profile 
	- b_n: coefficient of the Sersic profile that depends on the n index
		b_n is chosen so that half the total luminosity predicted by the Sersic profile comes from r <= r_e 
	- r_c: distance of the arc center related to the point where the arc is centered (pixel)
	- r_e: radius that encloses half of the total luminosity of the Sersic profile (pixel)
	- I_0: central surface brightness of the Sersic profile
	- ellip: ellipticity of the arc
	- h: factor used to define dr as a fraction of the r (dr = r/h)

	Output:
	- the value of the derivative of the arcellipse_sbmap_r_theta function at (r,theta)
	"""
	
	def df(r,theta_diff,n,b_n,r_c,r_e,I_0,ellip,h=1.e4):
		dr = r / h
		return ( f(r+dr/2.,theta_diff,n,b_n,r_c,r_e,I_0,ellip) - f(r-dr/2.,theta_diff,n,b_n,r_c,r_e,I_0,ellip) )/dr

	return df


def I_0_integral(n,r_e,r_c,theta_0,ellip):

	'''
	Computes the definite double integral of the Sersic profile (where I_0 = 1) used to define 
	the profile constant I_0 in terms of the total signal of the arc.
	The Sersic profile is defined as I = I_0 * exp(-b_n * (r_prof/r_e)**(1/n)), where
	I_0 is the central surface brightness and r_prof comes from the ArcEllipse prescription.

	Input:
	- n: Sersic index
	- r_e: radius that encloses half of the total luminosity of the Sersic profile (pix)
	- r_c: distance of the arc center related to the point where the arc is centered (pix)
	- theta_0: orientation of the arc center related to x-axis 
	- ellip: ellipticity of the arc

	Output:
	- value of the definite double integral of the Sersic profile (where I_0 = 1)
	'''

	b_n = get_b_n(n)
	
	return dblquad(lambda r, theta: r * math.exp(- b_n * ( (math.sqrt((r_c * (theta - theta_0) * \
		(1. - ellip))**2 + (r - r_c)**2 ) ) / r_e)**(1. / n)),theta_0 - math.pi, theta_0 + math.pi, lambda r: 0, lambda r: 7*r_c)[0]



def create_arcellipse_sbmap(x_0,y_0,theta_0,n,b_n,r_c,r_e,I_0,ellip,mag_zpt,pix_size,x_arc,y_arc):
	"""	
	Function to create the surface brightness map following the ArcEllipse prescription 
	and digitalize ("pixelize") it to an image with the same size as the input image. 

	Input:
	- x_0: x position of the point where the arc is centered (pixel)
	- y_0: y position of the point where the arc is centered (pixel)
	- theta_0: orientation of the arc center related to the x-axis (radians)
	- n: index that describes the shape of the Sersic profile 
	- b_n: coefficient of the Sersic profile that depends on the n index
		b_n is chosen so that half the total luminosity predicted by the Sersic profile comes from r <= r_e 
	- r_c: distance of the arc center related to the point where the arc is centered (pixel)
	- r_e: radius that encloses half of the total luminosity of the Sersic profile (pixel)
	- I_0: central surface brightness of the Sersic profile
	- ellip: ellipticity of the arc
	- mag_zpt: zero-point magnitude of the image
	- pix_size: size of the "new" pixel for the digitalization
	- dim_x: x dimension of the image
	- dim_y: y dimension of the image

	Output:
	- (ndarray,s_total): resultant array of the digitalization of the surface brightness map and the
		value of the total signal distributed over the new array
	"""
	# Determine the size of the postage stamp of the arc, which is equal to a factor (ps_factor) of the lenght of the arc
	# l = w / (1-ellip); w = 2b; b=r_e
	ps_size = int(ps_factor *  2.* r_e / (1. - ellip))
	print ps_size,ps_factor,r_e / (1. - ellip)
	s = numpy.zeros((ps_size+1,ps_size+1))
	s_total = []

	m = int(1. / pix_size)
	if m == 0:
		m = 1

	h = 1. / m
	
	ps_xmin = int(x_arc - ps_size/2.)
	ps_xmax = int(x_arc + ps_size/2.)
	ps_ymin = int(y_arc - ps_size/2.)
	ps_ymax = int(y_arc + ps_size/2.)	

	x_new = -1. 

	
	for x in range(ps_xmin,ps_xmax):
		x_new = x_new + 1
		y_new = -1.

		for y in range(ps_ymin,ps_ymax):	
			y_new = y_new + 1

	
			if x == 0:
				x_min = 0
			else:
				x_min = x #- 0.5  # shift -0.5 to take the lower value of the pixel, considering that the pixel value is in the center of the pixel
			if y == 0:
				y_min = 0
			else:
				y_min = y #- 0.5  # shift -0.5 to take the lower value of the pixel
			if x == ps_xmax:
				x_max = ps_xmax
			else:
				x_max = x #+ 0.5 # shift +0.5 to take the upper value of the pixel
			if y == ps_ymax:
				y_max = ps_ymax
			else:
				y_max = y #+ 0.5 # shift +0.5 to take the upper value of the pixel


		       	theta = math.atan2((y - y_0),(x - x_0))
		
			r = math.sqrt((x - x_0)**2 + (y - y_0)**2)


			theta_diff = abs(theta - theta_0)
			if (theta_diff >= (math.pi)):
				theta_diff = 2. * math.pi - theta_diff


		       	if (r_c - 10.*r_e)<= r <= (r_c+10.*r_e): #and theta_diff <= (math.pi/4.):
				
			       	s[y_new][x_new] = sum([sum([arcellipse_sbmap_x_y((x_min+i*h),(y_min+j*h),x_0,y_0,theta_0,n,b_n,r_c,r_e,I_0,ellip) for j in range(0,m)]) for i in range(0,m)])/float(m)/float(m)
				s_total.append(s[y_new][x_new])	
				

		       	else:
				
				s[y_new][x_new]= 0.
	
				s_total.append(s[y_new][x_new])
	
	s_sum = sum(s_total)

	return s,s_sum,ps_size


def paint_arcs(params):
	"""	
	Function to digitalize the arc surface brightness map and write it on a HDU.
    
	It computes the arcellipse_sbmap_r_theta derivative to find the size of the "new" pixel for 
	digitalization. Then it calls the function create_arcellipse_sbmap to create the surface 
	brightness map and digitalize it. 
    
 	This function returns the HDU resultant of the digitalization of the surface brightness map.
 	The parameters of the arc are added to the HDU header. 

	Input:
	- params: list with all necessary parameters (ellip,n,r_e,r_c,theta_0,mag,seeing,x_0,y_0,mag_zpt,dim_x,dim_y)
	
	Output:
	- hdu: resultant hdu of the digitalization of the surface brightness map 
	"""

	ellip = 1. - (1./float(params[0]))
	n = float(params[1])
	r_e = convert_arcsec_pix(float(params[2]),pix_scale)
	r_c = convert_arcsec_pix(float(params[3]),pix_scale)
	theta_0 = math.radians(float(params[4]))
	mag = float(params[5])
	x_0 = float(params[6])
	y_0 = float(params[7])
	mag_zpt = float(params[8])

	x_arc = x_0 + r_c * math.cos(theta_0)
	y_arc = y_0 + r_c * math.sin(theta_0)

	print ellip,n,r_e,r_c,theta_0,mag,x_0,y_0,mag_zpt,x_arc,y_arc

	signal_total = 10. ** (-0.4 * (mag - mag_zpt))

	I_0 = signal_total / I_0_integral(n,r_e,r_c,theta_0,ellip)

	b_n = get_b_n(n)

	d_arcellipse = derivative(arcellipse_sbmap_r_theta)(r_c+r_e,0,n,b_n,r_c,r_e,I_0,ellip)
	
	epsilon = 0.5

	pix_size = abs(arcellipse_sbmap_r_theta(r_c+r_e,0,n,b_n,r_c,r_e,I_0,ellip)/d_arcellipse * epsilon)

	sbmap,s_sum,ps_size = create_arcellipse_sbmap(x_0,y_0,theta_0,n,b_n,r_c,r_e,I_0,ellip,mag_zpt,pix_size,x_arc,y_arc)
	
	error = abs(s_sum - signal_total)/signal_total


	hdu=pyfits.PrimaryHDU(sbmap)

	hdr = hdu.header

	hdr.update('ALW',params[0],comment='Arc lenght-width ratio')
	hdr.update('ARC',params[3],comment='Arc curvature radius (arcsec)')
	hdr.update('ATHETA',params[4],comment='Arc theta (degree)')
	hdr.update('AMAG',mag,comment='Arc magnitude')
	hdr.update('MAG_ZPT',mag_zpt,comment='Zero point magnitude')	
	hdr.update('NINDEX',params[1],comment='Arc Sersic profile index')
	hdr.update('RE',params[2],comment='Arc Sersic effective radius (arcsec)')
	hdr.update('XARC',x_arc,comment='Arc x center')
	hdr.update('YARC',y_arc,comment='Arc y center')
	
	hdr.update('ASIG',signal_total,comment='Total signal of the arc (from the input)')
	hdr.update('ASIGDIST',s_sum,comment='Signal of the arc distributed by PaintArcs')
	hdr.update('FLUXERR',error,comment='Error on the distributed arc flux')
	
	return hdu



def run_paint_arcs(params,seeing,background):

	hdu = paint_arcs(params)
	
	arc_array = hdu.data
	arc_header = hdu.header
	
	pyfits.writeto("arc.fits",arc_array,arc_header)
	
	dim_x = arc_header["NAXIS1"]
	dim_y = arc_header["NAXIS2"]
	
	# Convolving the arc
	seeing = convert_arcsec_pix(seeing,pix_scale)	
	sigma = seeing / (2.*math.sqrt(2.*math.log(2.)))
	arc_header.update('APSF',seeing,comment='Arc psf (arcsec)')	
	
	arc_conv_array = gauss_convolution.fftconvolve_image(arc_array,3.,sigma)
	pyfits.writeto("arc_conv.fits",arc_conv_array,arc_header)
		
	# Adding Poisson noise
	arc_conv_noise_array = add_noise.add_poisson_noise(arc_conv_array)
	pyfits.writeto("arc_conv_noise.fits",arc_conv_noise_array,arc_header)	
	
	# Creating background image with noise
	bkg_array = background * numpy.ones((dim_y,dim_x)) 
	
	bkg_noise_array = add_noise.add_poisson_noise(bkg_array)
	
	pyfits.writeto("background.fits",bkg_noise_array)
	
	# Adding arc to the background image
	arc_conv_noise_bkg = arc_conv_noise_array + bkg_noise_array
	arc_header.update('BKG',background,comment='Background mean value')	
	
	pyfits.writeto("arc_conv_noise_bkg.fits",arc_conv_noise_bkg,arc_header)	

	








