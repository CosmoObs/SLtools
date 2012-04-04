#!/usr/bin/python
# ==================================
# Authors:
# Cristina Furlanetto - furlanetto.cristina@gmail.com
# ==================================

"""Module to convert RA and DEC formats."""

##@package convert_ra_dec
#
# This package contains functions to convert RA and DEC formats, 
# from sexagesimal to decimal format and vice-versa.
#
# 'convert_sexag2deg' converts from sexagesimal format (hh:mm:ss.ss or dd:mm:ss.ss) to decimal degrees format.
#
# 'convert_deg2sexag' converts from decimal degrees format to sexagesimal format (hh:mm:ss.ss or dd:mm:ss.ss).
#
# Both, convert_sexag2deg and convert_deg2sexag, have two arguments, RA and DEC, which can be given in 
# hr:min:sec or deg:min:sec format or decimal degree format. Depending on the format, the package calls 
# the appropriate function for format conversion.
#
# Executable package: YES
#
# To run the code, just type:
#
# > python convert_ra_dec.py --ra=10 --dec=20
#
# or
#
# > python convert_ra_dec.py --ra=00:40:00 --dec=20:00:00
#
# To see the package help message, type:
#
# > python convert_ra_dec.py --help
#

import string
import os
import sys


def convert_sexag2deg(ra,dec):

	"""Converts from sexagesimal format (hh:mm:ss.ss or dd:mm:ss.ss) to decimal degrees format.
	
	Input:
	- ra : (ra_h,ra_min,ra_sec) list with hr, min and sec components of RA in sexagesimal format 
	- dec: (dec_d,dec_min,dec_sec) list with deg, min and sec components of DEC in sexagesimal format 
	
	Output:
	- ra,dec : RA and DEC in decimal degree format.

	"""
	
	ra_h,ra_min,ra_sec = ra
	dec_d,dec_min,dec_sec = dec

	dec_signal = "+"
	if "-" in dec_d:
		dec_signal = "-"
		dec_d = -1. *float(dec_d)

	

	ra_deg = (((((float(ra_sec) / 60.0) + float(ra_min)) / 60.0) + float(ra_h)) / 24.0) * 360.0
	dec_deg = ((((float(dec_sec) / 60.0) + float(dec_min)) / 60.0) + float(dec_d))
	
	if dec_signal == "-":
		dec_deg = -1. * dec_deg
	
#	print "%.6f  %.6f\n" %(ra_deg,dec_deg)
	return ra_deg,dec_deg


def convert_deg2sexag(ra,dec):

	"""Converts from decimal degree format to sexagesimal format (hh:mm:ss.ss or dd:mm:ss.ss).
	
	Input:
	- ra : RA in decimal degree format.
	- dec: DEC in decimal degree format.
	
	Output:
	- ra : (ra_h,ra_min,ra_sec) list with hr, min and sec components of RA in sexagesimal format 
	- dec: (dec_d,dec_min,dec_sec) list with deg, min and sec components of DEC in sexagesimal format 
	
	"""

	ra = float(ra)
	dec = float(dec)
	
	signal = "+"
	if dec < 0:
		signal = "-"
		dec = dec*(-1)
	
	if ra < 0:
		ra = ra + 360.
	
	aux1 = ra/15.
	ra_h = int(aux1)
	aux2 = (aux1- float(ra_h))*60
	ra_min = int(aux2)
	ra_sec = (aux2-float(ra_min))*60
	
	ra_conv = (ra_h,ra_min,ra_sec)

	dec_d = int(dec)	
	aux3 = (dec-float(dec_d))*60.
	dec_min = int(aux3)
	dec_sec = (aux3-float(dec_min))*60.

	if signal == "-":
		dec_d = dec_d*(-1)


	dec_conv = (dec_d,dec_min,dec_sec)
	
#	print "%d:%d:%f %d:%d:%f" %(ra_h,ra_min,ra_sec,dec_d,dec_min,dec_sec) 
	return ra_conv,dec_conv

#\cond	
if __name__ == "__main__" :

	from optparse import OptionParser

	parser = OptionParser()

	parser.add_option('--ra',
			  dest='ra', default=None,
			  help='Right ascension value in hh:mm:ss.ss or decimal degree')
	parser.add_option('--dec',
			  dest='dec', default=None,
			  help='Declination value in dd:mm:ss.ss or decimal degree')


	(opts,args) = parser.parse_args()

	ra = opts.ra
	dec = opts.dec
	


	if ( opts=={} ):
		print""
		parser.print_help()
		print ""
		sys.exit(1)
	if (ra == None or dec == None ):
		print ""
		parser.print_help()
		print ""
		sys.exit(1)

	# Look if the input is in sexagesimal fomat (xx:xx:xx)
	if ":" in ra and dec:
		ra = string.split(ra,sep=":")
		dec = string.split(dec,sep=":")
		convert_sexag2deg(ra,dec)
	else:
		
		convert_deg2sexag(ra,dec)
	

	sys.exit(0)

#\endcond 


