import numpy as np
import matplotlib.pyplot as pyplot
import string
import random
import astropy.io.fits as pyfits

from sltools.image.imcp import cutout

# get parameters from file


def plot_arcellipse_contour(params_file):
	"""
	Function to plot contour of an arcellipse.
	
	Input:
	 - params_file  <file>    : Output file of the function sltools.arcellipse.arcfitting. Contains
				    the arcellipse parameters
	 - input_files.txt <file> : each line contains the name of a output file from arcellipse 
				    (params_file)

	Output:
	 - <file> : the plots of each
	"""

	image_name, x_0, y_0, r_c, ellip, I_0, r_e, theta_0, n, chi_2_min = np.loadtxt(params_file, dtype='str')

	x_0, y_0, r_c, ellip, theta_0, r_e = float(x_0), float(y_0), float(r_c), float(ellip), float(theta_0), float(r_e)

	#    * x0: x position of the point where the arc is centered (pixel)
	#    * y0: y position of the point where the arc is centered (pixel)
	#    * r_c: distance of the arc center related to the (x0,y0) (pixel)
	#    * theta0: orientation of the arc center related x axis (radians)
	#    * ellip: ellipticity of the arc # e = 1 - 1/(L/W) => 
	#    * I_0: central surface brightness of the Sersic profile
	#    * r_e: radius that encloses half of the total luminosity of the Sersic profile (pixel)
	#    * n: index that describes the shape of the Sersic profile
	#    * I_0: intensity of the Sersic profile at the center of the arc (r = 0) 

	# CORTE EM ellip
	if ellip > 0.99:
		ellip = 0.99


	# determine 'a' and 'b' from 'ellip'

	b = r_e
	a = b/(1. - ellip)

	# make sure we are using the correct limits
#	if a > np.pi * r_c:
#		a = np.pi * r_c

#	if b > r_c:
#		b = r_c

	if theta_0 > 2*np.pi:
		n = int( theta_0/(2*np.pi) )
		theta_0 = theta_0 - (n*2*np.pi)

	# =========================================

	theta_i = theta_0 - a/r_c
	theta_f = theta_0 + a/r_c
	theta = np.linspace(theta_i, theta_f, 1000)

	# arcellipse equation
	r_ext = r_c + b * np.sqrt( 1 - (r_c *(theta_0 - theta)/a)**2 )
	r_int = r_c - b * np.sqrt( 1 - (r_c *(theta_0 - theta)/a)**2 )

	pyplot.clf()
	f1 = pyplot.figure() # figsize = (15,15) )

#	# plot de todos os arcos
#	x_rand, y_rand = 1000*random.random(), 1000*random.random()

#	pyplot.plot(x_rand + r_ext*np.cos(theta), y_rand + r_ext*np.sin(theta) )
#	pyplot.plot( x_rand + r_int*np.cos(theta), y_rand + r_int*np.sin(theta) )
#	#======================

	pyplot.plot( x_0 + r_ext*np.cos(theta), y_0 + r_ext*np.sin(theta), 'r' )
	pyplot.plot( x_0 + r_int*np.cos(theta), y_0 + r_int*np.sin(theta), 'r' )

#	pyplot.plot( [r_ext[-1]*np.cos(theta[-1]), r_int[0]*np.cos(theta[0]) ], [ r_ext[-1]*np.sin(theta[-1]), r_int[0]*np.sin(theta[0])  ], 'k' )
#	pyplot.plot( [r_ext[0]*np.cos(theta)[0], r_int[-1]*np.cos(theta[-1]) ], [ r_ext[0]*np.sin(theta[0]), r_int[-1]*np.sin(theta[-1])  ], 'k' )
#	pyplot.plot( [r_ext[-1]*np.cos(theta[0]), r_int[-1]*np.cos(theta[-1]) ], [ r_ext[-1]*np.sin(theta[0]), r_int[-1]*np.sin(theta[-1])  ], 'k' )
#	pyplot.plot( [r_ext[0]*np.cos(theta)[-1], r_int[0]*np.cos(theta[0]) ], [ r_ext[0]*np.sin(theta[-1]), r_int[0]*np.sin(theta[0])  ], 'k' )

	pyplot.axis('equal')
	#pyplot.polar( theta , r_ext, 'r' )
	#pyplot.polar( theta , r_int, 'g' )
	#figname = params_file[:-4] + '.png'
	#pyplot.savefig(figname)

	fits_file = params_file[:-18] + '.fits'
	img = pyfits.getdata(fits_file)
	img_nonzero = np.where(img)
	pyplot.plot(img_nonzero[0], img_nonzero[1], 'g.'	)
	figname = params_file[:-18] + '_data_arcellipse.png'
	pyplot.savefig(figname)


#plot_arcellipse_contour('105_8647475120351217669_S82m35m_0_segm_118_fitted_params.txt')
#plot_arcellipse_contour('2_8647474693013963015_S82m9p__0_segm_61_fitted_params.txt')

input_files = open('input_files.txt', 'r').readlines()
for line in input_files:
        params_file = string.split(line)[0]
	print "Starting file %s" % params_file
	if len( open(params_file, 'r').readlines()  ) > 1: # so chamas os arquivos q foram gerados corretamente pelo arcfitting
		plot_arcellipse_contour(params_file) 









