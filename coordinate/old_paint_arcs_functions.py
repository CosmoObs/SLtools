from scipy.integrate import dblquad, quad
from scipy import constants
import math
import re
import string
import numpy

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
