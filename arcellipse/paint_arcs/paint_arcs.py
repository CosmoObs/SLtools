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


def paint_arcs(params,image_name):
	"""	
	Function to create the arc image. 
    
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
#	seeing = convert_arcsec_pix(float(params[6]),pix_scale)	
	x_0 = float(params[6])
	y_0 = float(params[7])
	mag_zpt = float(params[8])

	
	x_arc = x_0 + r_c * math.cos(theta_0)
	y_arc = y_0 + r_c * math.sin(theta_0)

	print ellip,n,r_e,r_c,theta_0,mag,x_0,y_0,mag_zpt,x_arc,y_arc

	signal_total = 10. ** (-0.4 * (mag - mag_zpt))

	I_0 = signal_total / I_0_integral(n,r_e,r_c,theta_0,ellip)
	print I_0,mag,mag_zpt,signal_total

	b_n = get_b_n(n)

	d_arcellipse = derivative(arcellipse_sbmap_r_theta)(r_c+r_e,0,n,b_n,r_c,r_e,I_0,ellip)
	
	epsilon = 0.5

	pix_size = abs(arcellipse_sbmap_r_theta(r_c+r_e,0,n,b_n,r_c,r_e,I_0,ellip)/d_arcellipse * epsilon)

	sbmap,s_sum,ps_size = create_arcellipse_sbmap(x_0,y_0,theta_0,n,b_n,r_c,r_e,I_0,ellip,mag_zpt,pix_size,x_arc,y_arc)
	
	error = abs(s_sum - signal_total)/signal_total

	img_name,img_ext = os.path.splitext(image_name)

	#Writing arc properties in a file
	arc_prop_file = open( img_name + "_arc.txt", "w")
	#arc_prop_file.write("ellip=%f, n=%f, r_e=%f, r_c=%f, theta=%f, mag=%f, seeing=%f \n" %(ellip,n,r_e,r_c,theta_0,mag,seeing))
	arc_prop_file.write("x_arc=%f, y_arc=%f \n" %(x_arc,y_arc))
	arc_prop_file.write("Signal total= %f, I_0= %f \n" %(signal_total,I_0))
	arc_prop_file.write("epsilon= %f \n" %epsilon)
	arc_prop_file.write("Derivative= %f \nn=1/pix_size= %f \nSignal_distributed= %f \nerror %f \n \n" %(d_arcellipse,int(1./pix_size),s_sum,error))

	#Writing the arc to a fits file
	arc_file = img_name + "_arc" + img_ext
	
	#image_header.update('ALW',l_w_header,comment='Arc lenght-width ratio')
	#image_header.update('ARC',r_c_header,comment='Arc curvature radius (arcsec)')
	#image_header.update('ATHETA',theta_0_header,comment='Arc theta (degree)')
	#image_header.update('AMAG',mag_header,comment='Arc magnitude')
	#image_header.update('ASEEING',seeing_header,comment='Arc seeing (arcsec)')
	#image_header.update('NINDEX',n_header,comment='Arc Sersic profile index')
	#image_header.update('RE',r_e_header,comment='Arc Sersic effective radius (arcsec)')
	#image_header.update('XARC',x_arc_header,comment='Arc x center')
	#image_header.update('YARC',y_arc_header,comment='Arc y center')

	pyfits.writeto(arc_file,sbmap)
	
	return sbmap















#-------------------------------------Funcoes preservadas para consulta futura ------------------


def run_paint_arcs(input_image,params_list):
	"""	
	Function to run PaintArcs and add an arc to an image.
    
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
	sbmap,ps_size,x_arc,y_arc = paint_arcs(params,image_name)

	ps_xmin = int(x_arc - ps_size/2.)
	ps_xmax = int(x_arc + ps_size/2.)+1
	ps_ymin = int(y_arc - ps_size/2.)
	ps_ymax = int(y_arc + ps_size/2.)+1



	#Convolving arc image with a gaussian kernel
	seeing = convert_arcsec_pix(float(params[6]),pix_scale)	
	sigma = seeing / (2.*math.sqrt(2.*math.log(2.)))
	
	fftconv_image = gauss_convolution.fftconvolve_image(sbmap,3.,sigma)
	
	#Adding convolved arc to image array
	#img_array = image_array + fftconv_image
	background_array = numpy.zeros((ps_size+1,ps_size+1))

	#Adding noise
	pyfits.writeto("arc_conv2.fits",fftconv_image)	
	pyfits.writeto("background.fits",background_array)
	noise_file_list = ["arc_conv2.fits","background.fits"]
	noisy_files_list = add_noise.add_noise(noise_file_list)
	
	for filename in noisy_files_list:
		arc_array = pyfits.getdata(filename)
				
		image_array[ps_ymin:ps_ymax,ps_xmin:ps_xmax] = image_array[ps_ymin:ps_ymax,ps_xmin:ps_xmax] + arc_array
	

	#Writing final image to a fits file
	final_image = paths['output'] + image_root + "_final" + image_extension

	pyfits.writeto(final_image,image_array,image_header)
	
	os.system("rm arc_conv2.fits")
	os.system("rm arc_convns_2.fits")	

	return image_array


#-----------------------------------------------------------------------------		
#Reading configuration file and input file to get parameters

#configfile = "paint_arc_config.cfg"
#paint_arcs_config = cp.read_config( configfile )

#input = paint_arcs_config['input']
#paths = paint_arcs_config['paths']
#output = paths['output']


#images_path = os.path.expanduser( paths['input_images'])
#catalogs_path = os.path.expanduser( paths['input_catalogs'])

# Scalars:
#input_file = str( input['input'])	
#pix_scale = float(input['dimpix'])


#images_list,params_list = read_input_file(images_path + input_file)


#for i in range(len(images_list)):

#	input_image = paths['input'] + images_list[i]

#	params = params_list[i]	

#	mag_zpt,dim_x,dim_y = get_header_parameter.get_header_parameter(input_image,'SEXMGZPT','NAXIS1','NAXIS2') 


#	if len(params) == 7:  
	#Add the position where the arc is centered (x_0,y_0) as the image center to params list in the 
	#case these parameters are not in the input file.
	
#		params.insert(7,dim_x / 2.)
#		params.insert(8,dim_y / 2.)
	
	#Add mag_zpt, dim_x and dim_y to parameters list
#	params.insert(9,mag_zpt)
#	params.insert(10,dim_x)
#	params.insert(11,dim_y)	

	
	#Run PaintArcs
#	img_array = run_paint_arcs(input_image,params_list)
