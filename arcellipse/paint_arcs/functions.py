from scipy.integrate import dblquad, quad
from scipy import constants
import math
import re
import string
import pywcs
import numpy

def read_input_file(input_file):

	"""
	input_fp : input filename, with the list of DC-images/sky-coordinates.
	bands_ls : list with image bands (filters) to be deal with.
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


def I_0_integral(n,r_e,r_c,theta_0,ellip):

	'''Computes the definite double integral of the Sersic profile (where I_0 = 1) used to define 
	the profile constant I_0 in terms of the total signal of the arc.
	The Sersic profile is defined as I = I_0 * exp(-b_n * (r_prof/r_e)**(1/n)), where
	I_0 is the central surface brightness and r_prof comes from the ArcEllipse prescription.

	@param n Sersic index
	@param r_e radius that encloses half of the total luminosity of the Sersic profile (pix)
	@param r_c distance of the arc center related to the cluster center (pix)
	@param theta_0 orientation of the arc center related to the cluster center 
	@param ellip ellipticity of the arc


	@return value of the definite double integral of the Sersic profile (where I_0 = 1)
	'''

	b_n = get_b_n(n)
	
	return dblquad(lambda r, theta: r * math.exp(- b_n * ( (math.sqrt((r_c * (theta - theta_0) * (1. - ellip))**2 + (r - r_c)**2 ) ) / r_e)**(1. / n)),theta_0 - math.pi, theta_0 + math.pi, lambda r: 0, lambda r: 5*r_c)[0]



def radec2xy(ra,dec,wcs):

	skycrd = numpy.array([[ra,dec]])

	pixcrd = wcs.wcs_sky2pix(skycrd,1)

	x = pixcrd[0][0]

	y = pixcrd[0][1]

	return x,y



def proper_dist_int(z_1,omega_m,omega_lambda):
	'''Computes the integrand for the proper distance integral for a LCDM model with omega_k = 0.

	'''
	
	return 1./math.sqrt(omega_m*((1.+z_1)**3)+omega_lambda)


def proper_dist(z,h_0,omega_m,omega_lambda):

	'''Computes the angular distance for a LCDM model with omega_k = 0

	'''

	c_h_0 = constants.c / (1000.* h_0) #(c in m/s; h_0 in km/s/Mpc; c_h_0 in Mpc)


	return c_h_0 * quad(lambda z_1: proper_dist_int(z_1,omega_m,omega_lambda),0,z)[0]

def angular_diameter_distance(z,h_0,omega_m,omega_lambda):
	

	return proper_dist(z,h_0,omega_m,omega_lambda)/ (1.+z)


