#!/usr/bin/python
#
#
##
# Module to draw geometric figures that mimic arcs and add them to images using the prescription of ArcEllipse.

##@module paint_arcs
#
#
# This module draws geometric figures that mimic arcs using the prescription of
# ArcEllipse, which consists in a deformation of an ellipse such that one of its# main axes becomes a circle segment.
#
# It uses Sersic profiles as the radial light distribution of the arcs. It also 
# performs Gaussian convolution and adds Poissonic noise to arc images.
# 
# The main goal of this code is that it allows the user to control all the input
# parameters of the arcs that are going to be added. 
#
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
from functions import *
from sltools.image import get_image_limits_pywcs
from sltools.coordinate import *
import pywcs


sys.path.append( os.path.abspath( '.' ))
cwd = os.getcwd()+'/'

#========================================================================
##
#	Function to compute the value of the surface brightness of the
#	(r,theta) point following the ArcEllipse prescription.
#
#	Input:
#	- r: radius (pixel)
#	- theta_diff: difference between the position angle and the arc center position angle - theta-theta_0 (radians)
#	- n: index that describes the shape of the Sersic profile 
#	- b_n: coefficient of the Sersic profile that depends on the n index; b_n is chosen so that half the total luminosity predicted by the Sersic 
#	profile comes from r <= r_e 
#	- r_c: distance of the arc center related to the cluster center (pixel)
#	- r_e: radius that encloses half of the total luminosity of the Sersic profile (pixel)
#	- I_0: central surface brightness of the Sersic profile
#	- ellip: ellipticity of the arc
#
#	Output:
#	- the value of the surface brightness of the (r,theta) point following the ArcEllipse prescription
#
#
def arcellipse_sbmap_r_theta(r,theta_diff,n,b_n,r_c,r_e,I_0,ellip):

	r_prof = math.sqrt((r_c * theta_diff * (1. - ellip))**2 + (r - r_c)**2 )

	exp_arg = - b_n * (r_prof / r_e) ** (1. / float(n)) 

	return I_0 * math.exp(exp_arg)


#========================================================================
##
#	Function to compute the value of the surface brightness of the
#	(x,y) pixel following the ArcEllipse prescription.
#
#	Input:
#	- x: x position (pixel)
#	- y: y position (pixel)
#	- x_0: x position of the cluster center or the point where the arc is centered (pixel)
#	- y_0: y position of the cluster center or the point where the arc is centered (pixel)
#	- theta_0: orientation of the arc center related to the cluster center or the point where the arc is centered (radians)
#	- n: index that describes the shape of the Sersic profile 
#	- b_n: coefficient of the Sersic profile that depends on the n index; b_n is chosen so that half the total luminosity predicted by the Sersic 
#	profile comes from r <= r_e 
#	- r_c: distance of the arc center related to the cluster center (pixel)
#	- r_e: radius that encloses half of the total luminosity of the Sersic profile (pixel)
#	- I_0: central surface brightness of the Sersic profile
#	- ellip: ellipticity of the arc
#
#	Output:
#	- the value of the surface brightness of the (x,y) point following the ArcEllipse prescription
#	
#
def arcellipse_sbmap_x_y(x,y,x_0,y_0,theta_0,n,b_n,r_c,r_e,I_0,ellip):

	theta = math.atan2((y - y_0),(x - x_0))
		
	r = math.sqrt((x - x_0)**2 + (y - y_0)**2)
	theta_diff = abs(theta - theta_0)
	if (theta_diff >= (math.pi)):
		theta_diff = 2. * math.pi - theta_diff

	r_prof = math.sqrt((r_c * theta_diff * (1. - ellip))**2 + (r - r_c)**2 )

	exp_arg = - b_n * (r_prof / r_e) ** (1. / float(n)) 

	return I_0 * math.exp(exp_arg)	


#========================================================================
##
#	Function to compute the numerical derivative of the arcellipse_sbmap_r_theta function at (r,theta). 
#
#	Input:
#	- r: radius (pixel)
#	- theta_diff: difference between the position angle (theta) and the arc center position angle (theta_0) - (radians)
#	- n: index that describes the shape of the Sersic profile 
#	- b_n: coefficient of the Sersic profile that depends on the n index; b_n is chosen so that half the total luminosity predicted by the Sersic 
#	profile comes from r <= r_e 
#	- r_c: distance of the arc center related to the cluster center (pixel)
#	- r_e: radius that encloses half of the total luminosity of the Sersic profile (pixel)
#	- I_0: central surface brightness of the Sersic profile
#	- ellip: ellipticity of the arc
#	- h: factor used to define dr as a fraction of the r (dr = r/h)
#
#	Output:
#	- the value of the derivative of the arcellipse_sbmap_r_theta function at (r,theta)
#	
#
def derivative(f):

	
	def df(r,theta_diff,n,b_n,r_c,r_e,I_0,ellip,h=1.e4):
		dr = r / h
		return ( f(r+dr/2.,theta_diff,n,b_n,r_c,r_e,I_0,ellip) - f(r-dr/2.,theta_diff,n,b_n,r_c,r_e,I_0,ellip) )/dr

	return df


#========================================================================
##
#	Function to create the surface brightness map following the ArcEllipse prescription and digitalize ("pixelize") it.
#
#	Input:
#	- x_0: x position of the cluster center or the point where the arc is centered (pixel)
#	- y_0: y position of the cluster center or the point where the arc is centered (pixel)
#	- theta_0: orientation of the arc center related to the cluster center or the point where the arc is centered (radians)
#	- n: index that describes the shape of the Sersic profile 
#	- b_n: coefficient of the Sersic profile that depends on the n index; b_n is chosen so that half the total luminosity predicted by the Sersic 
#	profile comes from r <= r_e 
#	- r_c: distance of the arc center related to the cluster center (pixel)
#	- r_e: radius that encloses half of the total luminosity of the Sersic profile (pixel)
#	- I_0: central surface brightness of the Sersic profile
#	- ellip: ellipticity of the arc
#	- mag_zpt: zero-point magnitude of the image
#	- pix_size: size of the "new" pixel for the digitalization
#	- dim_x: x dimension of the image
#	- dim_y: y dimension of the image
#
#	Output:
#	- (ndarray,s_total): resultant array of the digitalization of the surface brightness map and the value of the total signal distributed over the #	new array
#	
#	
#
def create_arcellipse_sbmap(x_0,y_0,theta_0,n,b_n,r_c,r_e,I_0,ellip,mag_zpt,pix_size,dim_x,dim_y):

	
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


#========================================================================
##
#	Function to create the arc image. It computes the arcellipse_sbmap_r_theta derivative to find the size of the "new" pixel 
#	for the digitalization. Then it calls create_arcellipse_sbmap to create the surface brightness map and digitalize it. The arc image 
#	is written in a fits file and the arc properties are written in a ascii file. 
#
#	Input:
#	- params: list with all necessary parameters
#	- image_name: name of the image where the arc will be added.
#	
#	Output:
#	- ndarray: resultant array of the digitalization of the surface brightness map 
#	
#	




def run_arcellipse(params,image_name):


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

	
	return sbmap

#-----------------------------------------------------------------------------		
#Reading configuration file and input file to get parameters

configfile = "paint_arc_config.cfg"
paint_arcs_config = cp.read_config( configfile )

input = paint_arcs_config['input']
paths = paint_arcs_config['paths']
output = paths['output']
#cosmol_params = paint_arcs_config['cosmological_background_params']

images_path = os.path.expanduser( paths['input_images'])
catalogs_path = os.path.expanduser( paths['input_catalogs'])

# Scalars:
input_file = str( input['input'])	
#halocat_file = str( input['halo_catalogs'])
#RAshift = float( input['rashift'])
#DECfactor = float( input['decfactor'])
#minimum_halo_mass = float( input['minimum_halo_mass'])
pix_scale = float(input['dimpix'])
#omega_m = float(cosmol_params['omega_m'])
#omega_l = float(cosmol_params['omega_l'])

#arc_loop = paint_arcs_config['loop_arc_properties_params']


images_list,params_list = read_input_file(images_path + input_file)


#-----------------------------------------------------------------------------	
#Run PaintArcs to each image in the image list:
for i in range(len(images_list)):

	image = paths['input'] + images_list[i]

	params = params_list[i]	

	mag_zpt,dim_x,dim_y = get_header_parameter(image,'SEXMGZPT','NAXIS1','NAXIS2') 

	dim_x = 100
	dim_y = 100

	if len(params) == 7:  
	#Add the position where the arc is centered (x_0,y_0) as the image center to params list in the case these
	#parameters are not in the input file
	
		#params.insert(7,dim_x / 2.)
		#params.insert(8,dim_y / 2.)
		params.insert(7,10)
		params.insert(8,10)
	
	
	params.insert(9,mag_zpt)
	params.insert(10,dim_x)
	params.insert(11,dim_y)	

	image_name, image_extension = os.path.splitext(images_list[i])	

	image_array,image_header = pyfits.getdata(image,header=True)	
	
	#Run arcellipse module
	sbmap = run_arcellipse(params,images_list[i])

	#Writing the arc to a fits file
	arc_file = paths['output'] + image_name + "_arc" + image_extension

	pyfits.writeto(arc_file,sbmap)

	#Convolving image
	seeing = convert_arcsec_pix(float(params[6]),pix_scale)	
	sigma = seeing / (2.*math.sqrt(2.*math.log(2.)))
	
	fftconv_image = gauss_convolution.fftconvolve_image(sbmap,3.,sigma)
	
	#Adding convolved arc to image array
#	img_array_conv = image_array + fftconv_image

#	Adding noise
#	pyfits.writeto("arc_conv2.fits",fftconv_image)	
#	noise_file_list = ["arc_conv2.fits"]
#	noisy_files_list = add_noise_2_image.add_noise_2_image(noise_file_list)
#
#	for filename in noisy_files_list:
#		arc_array = pyfits.getdata(filename)
#			
#		img_array = img_array + arc_array
#
#	out_conv = "arc_conv_dc4_" + str(i+1) + ".fits" 
#	out_noise = "arc_noise_dc4_" + str(i+1) + ".fits" 
#
#	os.system("mv arc_conv2.fits %s" %out_conv)
#	os.system("mv arc_convns_2.fits %s" %out_noise)	

	#Writing final image to a fits file
#	final_image = paths['output'] + image_name + "_final" + image_extension

#	pyfits.writeto(final_image,img_array_conv,image_header)
	



			


