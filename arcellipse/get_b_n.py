#!/usr/bin/python
# ==================================
# Authors:
# Cristina Furlanetto - furlanetto.cristina@gmail.com
# ==================================

"""Module to finde the b_n coefficient of the Sersic profile."""


##@package get_b_n
#
# This module finds the b_n coefficient of the Sersic profile. b_n is 
# chosen so that half the total luminosity predicted by the Sersic profile 
# comes from r <= r_e .
#
# This module uses the asymptotic expression given in the following references:
# Ciotti, L., Bertin, G. 1999, A&A, 352, 447
# MacArthur, L.A., Courteau, S., & Holtzman, J.A. 2003, ApJ, 582, 689
#
# This expression is more accurate thant the usual b_n = 1.9992*n - 0.3271. 
# For values of n > 0.36, the asymptotic expression is accurate to better than 10^-4.
#
# Executable package: YES
#
# To run the code, just type:
#
# > python get_b_n.py --n=1
#
# To see the package help message, type:
#
# > python get_b_n.py --help
#

import os
import sys


def get_b_n(n):

	''' Finds the b_n coefficient of the Sersic profile.
	This expression is more accurate thant the usual b_n = 1.9992*n - 0.3271
	Reference: Ciotti, L., Bertin, G. 1999, A&A, 352, 447
	MacArthur, L.A., Courteau, S., & Holtzman, J.A. 2003, ApJ, 582, 689
	b_n is chosen so that half the total luminosity predicted by the Sersic 
	profile comes from r <= r_e . 

	Input:
	-n: Sersic index 
	
	Output:
	- b_n: coefficient of the Sersic profile 	
	'''
	
	n = float(n)
	
	if n < 0.36:
		print "Sorry. This module just computes the b_n coefficient for n > 0.36."
		return 

	#if n == 0.0: 
	#	return -0.3271
		
	b_n = 2.0*n - 1.0/3.0 + 4.0/(405.0*n) + 46.0/(25515.0*(n**2) ) + 131.0/(1148175.0*(n**3)) \
		- 2194697.0/(30690717750.0*(n**4))

	print b_n
	
	return b_n

#\cond	
if __name__ == "__main__" :

	from optparse import OptionParser

	parser = OptionParser()

	parser.add_option('--n',
			  dest='n', default=None,
			  help='Sersic index')


	(opts,args) = parser.parse_args()

	n = opts.n


	if ( opts=={} ):
		print""
		parser.print_help()
		print ""
		sys.exit(1)
	if (n == None):
		print ""
		parser.print_help()
		print ""
		sys.exit(1)
	
	get_b_n(n)
	
	sys.exit(0)

#\endcond 
