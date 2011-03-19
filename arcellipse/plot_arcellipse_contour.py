import numpy as np
import matplotlib.pyplot as pyplot
import string
import random

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

	# determine 'a' and 'b' from 'ellip'

	b = r_e
	a = b/(1. - ellip)


	theta = np.linspace(0,2*np.pi, 1000)
	# arcellipse equation
	r_ext = r_c + b * np.sqrt( 1 - (r_c *(theta_0 - theta)/a)**2 )
	r_int = r_c - b * np.sqrt( 1 - (r_c *(theta_0 - theta)/a)**2 )

	f1 = pyplot.figure() # figsize = (15,15) )

#	# plot de todos os arcos
#	x_rand, y_rand = 1000*random.random(), 1000*random.random()

#	pyplot.plot(x_rand + r_ext*np.cos(theta), y_rand + r_ext*np.sin(theta) )
#	pyplot.plot( x_rand + r_int*np.cos(theta), y_rand + r_int*np.sin(theta) )
#	#======================

	pyplot.plot( r_ext*np.cos(theta), r_ext*np.sin(theta), 'r' )
	pyplot.plot( r_int*np.cos(theta), r_int*np.sin(theta), 'r' )
	pyplot.axis('equal')
	#pyplot.polar( theta , r_ext, 'r' )
	#pyplot.polar( theta , r_int, 'g' )
	figname = params_file[:-4] + '.png'
	pyplot.savefig(figname)

#plot_arcellipse_contour('105_8647475120351217669_S82m35m_0_segm_118_fitted_params.txt')


#input_files = open('input_files.txt', 'r').readlines()
#for line in input_files:
#         params_file = string.split(line)[0]
#	 plot_arcellipse_contour(params_file)












