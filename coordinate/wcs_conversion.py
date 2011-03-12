#!/usr/bin/python
# ====================================================
# Authors:
# Cristina Furlanetto - furlanetto.cristina@gmail.com
# ====================================================

"""Module to convert sky coordinates to pixel coordinates and vice-versa."""

##@package wcs_conversion
#
# This package contains functions to transform sky coordinates (RA and Dec) to pixel coordinates 
# (x and y) and vice-versa. It is a simple implementation of pywcs functions.
#
# 'radec2xy' transforms sky coordinates (RA and Dec) to pixel coordinates (x and y). 
#
# 'xy2radec' transforms pixel coordinates (x and y) to sky coordinates (RA and Dec).
#
# Both radec2xy and xy2radec have three mandatory arguments: the FITS image name and the coordinate 
# values (RA and Dec our x and y). RA and Dec values must be in degrees.
#
#

import os
import sys
import pyfits
import pywcs
import numpy


def radec2xy(image,ra,dec):

	"""Transforms sky coordinates (RA and Dec) to pixel coordinates (x and y).
	
	Input:
	- image: FITS image name
	- ra : Right ascension value in degrees
	- dec: Declination value in degrees
	
	Output:
	- (x,y) : pixel coordinates

	"""

	hdulist = pyfits.open(image)

	wcs = pywcs.WCS(hdulist[0].header)
	
	skycrd = numpy.array([[ra,dec]])

	pixcrd = wcs.wcs_sky2pix(skycrd,1)

	x = pixcrd[0][0]

	y = pixcrd[0][1]

	print x,y

	return (x,y)
	
	


def xy2radec(image,x,y):

	""" Transforms pixel coordinates (x and y) to sky coordinates (RA and Dec).
	
	Input:
	- image: FITS image name
	- x: x value
	- y: y value
	
	Output:
	- (ra,dec) : Right ascension and declination values in degrees 

	"""

	hdulist = pyfits.open(image)
	
	wcs = pywcs.WCS(hdulist[0].header)

	pixcrd = numpy.array([[x,y]], numpy.float_)

	skycrd = wcs.wcs_pix2sky(pixcrd,1)

	ra = skycrd[0][0]

	dec = skycrd[0][1]

	print ra,dec

	return (ra,dec)	




