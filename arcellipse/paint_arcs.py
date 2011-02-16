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

import sys
import os
import numpy
import math
import pyfits
import sltools
from sltools.image import gauss_convolution
from sltools.image import add_noise
from sltools.image import get_header_parameter
from sltools.image import imcp
from sltools.io import config_parser as cp
from sltools.catalog import open_fits_catalog
from sltools.catalog.halos import *
import string
from functions import get_b_n,convert_arcsec_pix,convert_pix_arcsec # substituir por get_b_n, convert_arcsec_pix, convert_pix_arcsec
from sltools.image import get_image_limits_pywcs
from sltools.coordinate import *
import pywcs


sys.path.append( os.path.abspath( '.' ))
cwd = os.getcwd()+'/'

def read_input_file(input_file):

	""" Function to read the input file and get the images list and the parameters list.
	
	Input:
	- input_file: input file name
	
	Output:
	- images_list: list of images where arc will be added
	- params_list: list of parameters of the arcs that will be added
	"""

	file_in = open(input_file,'r')
	images_list = []
	params_list = []

	for line in file_in.readlines():

		if ( re.search('^$',line) ) : continue
		if ( re.search('\s*#',line) ) : continue
	
		in_list = string.split( line,sep=None )
		images_list.append(in_list[0])
		params_list.append(in_list[1:len(in_list)])


	return images_list,params_list



def arcellipse_sbmap_r_theta(r,theta_diff,n,b_n,r_c,r_e,I_0,ellip):
    """ Function to compute the value of the surface brightness of the
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
    """ Function to compute the value of the surface brightness of the
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
	theta_diff = abs(theta - theta_0)
	if (theta_diff >= (math.pi)):
		theta_diff = 2. * math.pi - theta_diff

	r_prof = math.sqrt((r_c * theta_diff * (1. - ellip))**2 + (r - r_c)**2 )

	exp_arg = - b_n * (r_prof / r_e) ** (1. / float(n)) 

	return I_0 * math.exp(exp_arg)	



def derivative(f):
    """ Function to compute the numerical derivative of the arcellipse_sbmap_r_theta function at (r,theta). 

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

	'''Computes the definite double integral of the Sersic profile (where I_0 = 1) used to define 
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
		(1. - ellip))**2 + (r - r_c)**2 ) ) / r_e)**(1. / n)),theta_0 - math.pi, theta_0 + math.pi, lambda r: 0, lambda r: 5*r_c)[0]



def create_arcellipse_sbmap(x_0,y_0,theta_0,n,b_n,r_c,r_e,I_0,ellip,mag_zpt,pix_size,dim_x,dim_y):
    """	Function to create the surface brightness map following the ArcEllipse prescription 
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
	
	s = numpy.zeros((dim_y,dim_x))
	s_total = []

	m = int(1. / pix_size)
	if m == 0:
		m = 1

	h = 1. / m

	for x in range(dim_x):
		for y in range(dim_y):

			if x == 0:
				x_min = 0
			else:
				x_min = x - 0.5  # shift -0.5 to take the lower value of the pixel, considering that the pixel value is in the center of the pixel
			if y == 0:
				y_min = 0
			else:
				y_min = y - 0.5  # shift -0.5 to take the lower value of the pixel
			if x == dim_x:
				x_max = dim_x
			else:
				x_max = x + 0.5 # shift +0.5 to take the upper value of the pixel
			if y == dim_y:
				y_max = dim_y
			else:
				y_max = y + 0.5 # shift +0.5 to take the upper value of the pixel

			s[y][x] = sum([sum([arcellipse_sbmap_x_y((x_min+i*h),(y_min+j*h),x_0,y_0,theta_0,n,b_n,r_c,r_e,I_0,ellip) for j in range(0,m)]) for i in range(0,m)])/float(m)/float(m)
				
			s_total.append(s[y][x])
	
	s_sum = sum(s_total)

	return s,s_sum


def digitalize_arc_sbmap(params,image_name):
    """	Function to create the arc image. 
    
    It computes the arcellipse_sbmap_r_theta derivative to find the size of the "new" pixel for 
    digitalization. Then it calls the function create_arcellipse_sbmap to create the surface 
    brightness map and digitalize it. 
    
    The arc image is written in a fits file with suffix "_arc.fits" and the arc properties are 
    written in a ascii file with suffix "_arc.txt". This function returns the array resultant of 
    the digitalization of the surface brightness map.

	The arc image has the same size as the input image.

	Input:
	- params: list with all necessary parameters (ellip,n,r_e,r_c,theta_0,mag,seeing,x_0,y_0,mag_zpt,dim_x,dim_y)
	- image_name: name of the image where the arc will be added
	
	Output:
	- ndarray: resultant array of the digitalization of the surface brightness map 
	"""

	ellip = 1. - (1./float(params[0]))
	n = float(params[1])
	r_e = convert_arcsec_pix(float(params[2]),pix_scale)
	r_c = convert_arcsec_pix(float(params[3]),pix_scale)
	theta_0 = math.radians(float(params[4]))
	mag = float(params[5])
	seeing = convert_arcsec_pix(float(params[6]),pix_scale)	
	x_0 = float(params[7])
	y_0 = float(params[8])
	mag_zpt = float(params[9])
	dim_x = params[10]
	dim_y = params[11]	
	
	x_arc = x_0 + r_c * math.cos(theta_0)
	y_arc = y_0 + r_c * math.sin(theta_0)

	print ellip,n,r_e,r_c,theta_0,mag,x_0,y_0,mag_zpt,x_arc,y_arc

	signal_total = 10. ** (-0.4 * (mag - mag_zpt))

	I_0 = signal_total / I_0_integral(n,r_e,r_c,theta_0,ellip)

	b_n = get_b_n(n)

	d_arcellipse = derivative(arcellipse_sbmap_r_theta)(r_c+r_e,0,n,b_n,r_c,r_e,I_0,ellip)
	
	epsilon = 0.5

	pix_size = abs(arcellipse_sbmap_r_theta(r_c+r_e,0,n,b_n,r_c,r_e,I_0,ellip)/d_arcellipse * epsilon)

	sbmap,s_sum = create_arcellipse_sbmap(x_0,y_0,theta_0,n,b_n,r_c,r_e,I_0,ellip,mag_zpt,pix_size,dim_x,dim_y)
	
	error = abs(s_sum - signal_total)/signal_total

	img_name,img_ext = os.path.splitext(image_name)

	#Writing arc properties in a file
	arc_prop_file = open(paths['output'] + img_name + "_arc.txt", "w")
	arc_prop_file.write("ellip=%f, n=%f, r_e=%f, r_c=%f, theta=%f, mag=%f, seeing=%f \n" %(ellip,n,r_e,r_c,theta_0,mag,seeing))
	arc_prop_file.write("x_arc=%f, y_arc=%f \n" %(x_arc,y_arc))
	arc_prop_file.write("Signal total= %f, I_0= %f \n" %(signal_total,I_0))
	arc_prop_file.write("epsilon= %f \n" %epsilon)
	arc_prop_file.write("Derivative= %f \nn=1/pix_size= %f \nSignal_distributed= %f \nerror %f \n \n" %(d_arcellipse,int(1./pix_size),s_sum,error))

	#Writing the arc to a fits file
	arc_file = paths['output'] + img_name + "_arc" + img_ext

	pyfits.writeto(arc_file,sbmap)
	
	return sbmap


def run_paint_arcs(input_image,params_list):
	"""	Function to run PaintArcs and add an arc to an image.
    
	It calls the function digitalize_arc_sbmap to digitalize the arc surface brightness map and 
	create the arc image (with suffix "_arc.fits"). Then it convolves the arc image with a gaussian
	kernel, to add the seeing effect. After that, it add Poissonic noise to the arc image. Finally, 
	it adds the arc to the input image. The final image with the added arc has the suffix "_final.fits". 
    
 
	Input:
	- image_name: name of the image where the arc will be added
	- params_list: list with arc parameters (ellip,n,r_e,r_c,theta_0,mag,seeing,x_0,y_0,mag_zpt,dim_x,dim_y)

	
	Output:
	- ndarray: final arc array
	"""
	
	image_name = string.split(input_image,sep='/')[-1]
	image_root,image_extension = os.path.splitext(image_name) 

	image_array,image_header = pyfits.getdata(input_image,header=True)	
	
	#Digitalize and create arc image
	sbmap = digitalize_arc_sbmap(params,image_name)

	#Convolving arc image with a gaussian kernel
	seeing = convert_arcsec_pix(float(params[6]),pix_scale)	
	sigma = seeing / (2.*math.sqrt(2.*math.log(2.)))
	
	fftconv_image = gauss_convolution.fftconvolve_image(sbmap,3.,sigma)
	
	#Adding convolved arc to image array
	#img_array_conv = image_array + fftconv_image

	#Adding noise
	pyfits.writeto("arc_conv2.fits",fftconv_image)	
	noise_file_list = ["arc_conv2.fits"]
	noisy_files_list = add_noise_2_image.add_noise_2_image(noise_file_list)
	
	for filename in noisy_files_list:
		arc_array = pyfits.getdata(filename)
				
		img_array = img_array + arc_array
	

	#Writing final image to a fits file
	final_image = paths['output'] + image_root + "_final" + image_extension

	pyfits.writeto(final_image,img_array_conv,image_header)
	
	os.system("rm arc_conv2.fits")
	
	return img_array


#-----------------------------------------------------------------------------		
#Reading configuration file and input file to get parameters

configfile = "paint_arc_config.cfg"
paint_arcs_config = cp.read_config( configfile )

input = paint_arcs_config['input']
paths = paint_arcs_config['paths']
output = paths['output']


images_path = os.path.expanduser( paths['input_images'])
catalogs_path = os.path.expanduser( paths['input_catalogs'])

# Scalars:
input_file = str( input['input'])	
pix_scale = float(input['dimpix'])


images_list,params_list = read_input_file(images_path + input_file)


for i in range(len(images_list)):

	input_image = paths['input'] + images_list[i]

	params = params_list[i]	

	mag_zpt,dim_x,dim_y = get_header_parameter(input_image,'SEXMGZPT','NAXIS1','NAXIS2') 

	if len(params) == 7:  
	#Add the position where the arc is centered (x_0,y_0) as the image center to params list in the 
	#case these parameters are not in the input file.
	
		params.insert(7,dim_x / 2.)
		params.insert(8,dim_y / 2.)
	
	#Add mag_zpt, dim_x and dim_y to parameters list
	params.insert(9,mag_zpt)
	params.insert(10,dim_x)
	params.insert(11,dim_y)	

	
	#Run PaintArcs
	img_array = run_paint_arcs(input_image,params_list)