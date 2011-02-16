#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ==================================
# Authors:
# Pedro Ferreira - pferreira@dfte.ufrn.br
# Comments cleanup by MSSG
# Habib Dumet Montoya - habibdumet@gmail.com
# ==================================

""" """


##@package find_CC_new
# Determines the caustics and CC of a given lens
#
# Determines the lens curves using a iterating method using gravlens


import os
import logging
#from lens_parameters_new import lens_parameters_new
from sltools.gravlens.lens_parameters_new import lens_parameters_new
from sltools.geometry import separate_curves

import numpy as np
import matplotlib.pyplot as pyplot

# given the gravlens 'config file' runs gravlens to generate the 
# file with caustics and CC
#
def _critcurves_new(gravinput, caustic_CC_file, gravlens_input_file='gravlens_CC_input.txt'):
	fmag = open(gravlens_input_file, 'w')
	fmag.write(gravinput)
	fmag.write('plotcrit %s\n' % caustic_CC_file )
	fmag.close()
	os.system('gravlens %s > /dev/null' % gravlens_input_file )


## find_CC_new
# Determines the caustics and CC of a given lens
#
# Determines the lens curves using a iterating method using gravlens


#=======================================================================================================
def find_CC_new(lens_model, mass_scale, model_param_8, model_param_9, model_param_10, galaxy_position=[0,0], e_L=0, theta_L=0, shear=0, theta_shear=0, gravlens_params={}, caustic_CC_file='crit.txt', gravlens_input_file='gravlens_CC_input.txt'):
	""" 
	Determines the caustics and critical curves for a single component lens model. Uses an iterative
	procedure to determine the grid size and the minimum number of points on the curves.

	The CC size may vary considerably depending on the lens-source configuration. Because of this, 
	we developed an iterative method to find the CC, that involves the size of the grid used by 
	gravlens (gridhi1). The initial value of gridhi1 must be an upper limit for the size of the CC, 
	and find_CC_new will try to find the CC. Each time gravlens can't find the CC, we lower	gridhi1 
	by a factor of grid_factor (=5), in order to increase precision. This continues on until 
	gravlens finds the CC or the number of max_iterations_number (=20) iterations is reached. After 
	finding the CC, if it is composed of less than min_n_lines (=200) points, the code redefines 
	gridhi1 to a third of it (grid_factor2=3), usually increasing the number of points. With this 
	first determination of the CC, the code uses its scale (gridhi1_CC_factor(=2) times the distance
	of the furthest CC point to the origin) to reobtain the CC with the apropriate value for the 
	grid. If the CC were not found after all these attempts, this function returns "False" and an 
	error message. If all caustic points have coordinates < acceptable_res_limit (=2E-4), a warning 
	message will be generated, meaning that the precision of gravlens (limited by the number of 
	digits in the output file) is being reached. 

	Example:
	         x_caustic, y_caustic, x_CC, y_CC, gravlens_config = find_CC_new('nfw',1,1,0,0)

	Input:
	 - lens_model       <str> : Lens name (see gravlens manual table 3.1)
	 - mass_scale     <float> : Mass scale of the lens - "parameter 1"
	 - model_param_8  <float> : misc. lens parameter - often scale radio
	 - model_param_9  <float> : misc. lens parameter - often scale radio
	 - model_param_10 <float> : misc. lens parameter - often a power law index
	 - galaxy_position <list> : [x,y] position of the lens
	 - e_L            <float> : lens ellipticity (default=0)
	 - theta_L        <float> : lens position angle (in degrees) with respect to the vertical 
					  (counterclockwise)
	 - shear          <float> : external shear amplitude
	 - theta_shear    <float> : external shear direction (in degrees)
	 - gravlens_params  <dic> : Contains the keys and values of the gravlens configuration 
					  (see default parameters at function set_gravlens_default,
                                          inside lens_parameters_new)
	 - caustic_CC_file        <str> : name of the output file with the caustic and CC positions
	 - gravlens_input_file    <str> : name of the input file used to run gravlens

	Output:
	 - <list>     : [x_caustic, ycaustic], with x_caustic and ycaustic being the lists with the 
			caustic coordinates
	 - <list>     : [x_CC, y_CC], with x_CC and y_CC being the lists with the CC coordinates
	 - <dict>     : all configuration variables used for running gravlens (including gridhi1)
	 - <file>     : file named 'caustic_CC_file' with the caustic and CC positions

	"""

	# definitions of some control params
	grid_factor = 5.
	max_iterations_number = 20
	grid_factor2 = 3.
	min_n_lines = 200
	gridhi1_CC_factor = 2.
	acceptable_res_limit = 2E-4
	# ==================================

	inputlens, setlens, gravlens_params_updated = lens_parameters_new(lens_model, mass_scale, model_param_8, model_param_9, model_param_10, galaxy_position, e_L, theta_L, shear, theta_shear, gravlens_params) # inputlens is the gravlens input (sets general parameters and the lens parameters)	# setlens is the gravlens input line that concerns the lens (ex: nfw 1 0 ...)
	logging.debug('Determined the strings that defines the mass model in gravlens through the \'lensparameter\' function')
	#--------------------------------------

	#-----------------------------
	_critcurves_new(inputlens, caustic_CC_file, gravlens_input_file) # gets the critical curves (crit.txt file)
	logging.debug('Got the critical curves (%s file) with function \"critcurves\"' % caustic_CC_file)
	#-----------------------------

	# ITERATION ON gridhi1 (TO GET BETTER RESOLUTION ON CC)
	counter = 0
	while os.path.isfile('./' + caustic_CC_file) == False and counter < max_iterations_number: # looks for the C.C. file (crit.txt) - if this file does not exist, add a note to arccatalog and quit
		gravlens_params_updated['gridhi1'] = float( gravlens_params_updated['gridhi1'] ) / grid_factor; # gridhi1 /= 5. 

		lens_par_out = lens_parameters_new(lens_model, mass_scale, model_param_8, model_param_9, model_param_10, galaxy_position, e_L, theta_L, shear, theta_shear, gravlens_params_updated) # inputlens is the gravlens input (sets general parameters and the lens parameters)

		inputlens, setlens, gravlens_params_updated = lens_par_out
		_critcurves_new(inputlens, caustic_CC_file, gravlens_input_file) # gets the critical curves (crit.txt file)
		counter += 1
	logging.debug( '(Number of iterations on gridhi1 = %d)' % (counter) )

	if os.path.isfile('./' + caustic_CC_file) == False: # looks for the C.C. file (crit.txt) - if this file does not exist, returns 'False'
		logging.error( 'Gravlens could not determine the caustics and critical curves. Try changing the input units (for example, from arcseconds to miliarcsecond), because gravlens outputs the critical curves with only 6 digits.' % (float( gravlens_params_updated['gridhi1'] ),  counter) )
		return False

 	if len(open('./' + caustic_CC_file).readlines() ) < min_n_lines: # these correspond to critical cases when the CC have too few points.
		gravlens_params_updated['gridhi1'] = float( gravlens_params_updated['gridhi1'] ) / grid_factor2 # gridhi1 /= 3. 

		lens_par_out = lens_parameters_new(lens_model, mass_scale, model_param_8, model_param_9, model_param_10, galaxy_position, e_L, theta_L, shear, theta_shear, gravlens_params_updated) # inputlens is the gravlens input (sets general parameters and the lens parameters)
		inputlens, setlens, gravlens_params_updated = lens_par_out

		_critcurves_new(inputlens, caustic_CC_file, gravlens_input_file) # gets the critical curves (crit.txt file)
	logging.debug( 'gridhi1 = %f and number of iterations on gridhi1 = %d' % (float( gravlens_params_updated['gridhi1'] ),  counter) )


	x1, y1, u1, v1 = np.loadtxt(caustic_CC_file, usecols = (0,1,2,3), unpack=True) # CC_x, CC_y, caustic_x, caustic_y

	#-----------------------------------------------------------------------------------------------------------
	# redefine gridhi1 according to the CC size

	index_CC = np.argmax(x1**2 + y1**2)

	gravlens_params_updated['gridhi1'] =  grid_factor * ( (x1[index_CC]**2 + y1[index_CC]**2)**0.5 )
	#-----------------------------------------------------------------------------------------------------------
	# get CC with best value for the grid
	lens_par_out = lens_parameters_new(lens_model, mass_scale, model_param_8, model_param_9, model_param_10, galaxy_position, e_L, theta_L, shear, theta_shear, gravlens_params_updated) # inputlens is the gravlens input (sets general parameters and the lens parameters)
	inputlens, setlens, gravlens_params_updated = lens_par_out

	_critcurves_new(inputlens, caustic_CC_file, gravlens_input_file) # gets the critical curves (crit.txt file)

	x1, y1, u1, v1 = np.loadtxt(caustic_CC_file, usecols = (0,1,2,3), unpack=True) # CC_x, CC_y, caustic_x, caustic_y

	# check if the precision is ok (gravlens outputs coordinates with only 6 decimal places)
	if max( max(np.abs(u1)), max(np.abs(v1)) ) < acceptable_res_limit:
		logging.warning('The caustics seem to be determined with low resolution')

	return x1, y1, u1, v1, gravlens_params_updated




#===========================================================================



## plot_CC
# Plots the caustics and critical curves
#
def plot_CC(tan_caustic_x, tan_caustic_y, rad_caustic_x, rad_caustic_y, tan_CC_x, tan_CC_y, rad_CC_x, rad_CC_y, plot_filename, show_plot=0):
	pyplot.clf()
	f1 = pyplot.figure(1,figsize = (15,15) )
	pyplot.subplot(121) # (numRows,numCols, plotNum)
	pyplot.plot(tan_caustic_x, tan_caustic_y, color='red',   marker='.', markersize=2 )
	pyplot.plot(rad_caustic_x, rad_caustic_y, color='black', marker='.', markersize=2 )
	pyplot.axis('equal')
	pyplot.xlabel('x', fontsize = 15)
	pyplot.ylabel('y', fontsize = 15)
	pyplot.title('Caustics'  , fontsize = 20)

	pyplot.subplot(122) # (numRows,numCols, plotNum)
	pyplot.plot(tan_CC_x, tan_CC_y, color='red',   marker='', markersize=2 )
	pyplot.plot(rad_CC_x, rad_CC_y, color='black', marker='', markersize=2 )
	pyplot.axis('equal')
	pyplot.xlabel('x', fontsize = 15)
	pyplot.ylabel('y', fontsize = 15)
	pyplot.title('Critical Curves'  , fontsize = 20)
	if show_plot == 0:
		pyplot.savefig(plot_filename)
	if show_plot == 1:
		pyplot.show() 
 
## run_find_CC_new
# Finds the caustics and CC for a given lens model, separates the radial from the tangential and plots the curves
#
def run_find_CC(lens_model, mass_scale, model_param_8, model_param_9, model_param_10, galaxy_position=[0,0], e_L=0, theta_L=0, shear=0, theta_shear=0, gravlens_params={}, caustic_CC_file='crit.txt', gravlens_input_file='gravlens_CC_input.txt', rad_curves_file='crit_rad.txt', tan_curves_file='crit_tan.txt', curves_plot='crit_curves.png', show_plot=0, write_to_file=0):
	""" 
	This is a pipeline that runs 'find_CC_new', 'separate_CC' and
	'plot_CC'. For details of these functions, see their documentation.

	Input:
	 - lens_model             <str> : Lens name (see gravlens manual table 3.1)
	 - mass_scale           <float> : Mass scale of the lens - "parameter 1"
	 - model_param_8        <float> : misc. lens parameter - often scale radio
	 - model_param_9        <float> : misc. lens parameter - often scale radio
	 - model_param_10       <float> : misc. lens parameter - often a power law index
	 - galaxy_position <list float> : [x,y] position of the lens
	 - e_L                  <float> : lens ellipticity (default=0)
	 - theta_L              <float> : lens position angle (in degrees) with respect to the vertical 
					  (counterclockwise)
	 - shear                <float> : external shear amplitude
	 - theta_shear          <float> : external shear direction (in degrees)
	 - gravlens_params        <dic> : Contains the keys and values of the gravlens configuration 
					  (see default parameters at function set_gravlens_default,
                                          inside lens_parameters_new)
	 - caustic_CC_file        <str> : name of the output file with the caustic and CC positions
	 - gravlens_input_file    <str> : name of the input file used to run gravlens

	Output:
	 - <list>     : [x_caustic, ycaustic], with x_caustic and ycaustic being the lists with the caustic coordinates
	 - <list>     : [x_CC, y_CC], with x_CC and y_CC being the lists with the CC coordinates
	 - <dict>     : all configuration variables used for running gravlens (including gridhi1)
	 - <file>     : file named 'caustic_CC_file' with the caustic and CC positions
	 - <file>     : separated radial curves (CC + caustics) - rad_curves_file
	 - <file>     : separated tangential curves (CC + caustics) - tan_curves_file
	 - <file>     : curves_plot
	"""

	x_caustic, y_caustic, x_CC, y_CC, gravlens_params_updated = find_CC_new(lens_model, mass_scale, model_param_8, model_param_9, model_param_10, galaxy_position, e_L, theta_L, shear, theta_shear, gravlens_params, caustic_CC_file, gravlens_input_file)


	x1, y1, u1, v1, x2, y2, u2, v2 = np.loadtxt(caustic_CC_file, comments='#', unpack=True)

	# Separating the critical cuves

	rad_CC_x,rad_CC_y, tan_CC_x, tan_CC_y = separate_curves(x1, y1, x2, y2)

	# separating the caustics curves

	rad_caustic_x, rad_caustic_y, tan_caustic_x, tan_caustic_y  = separate_curves(u1, v1, u2, v2)

	if curves_plot != 0: # curves_plot = 0 means you don't want any plots.
		plot_CC(tan_caustic_x, tan_caustic_y, rad_caustic_x, rad_caustic_y, tan_CC_x, tan_CC_y, rad_CC_x, rad_CC_y, curves_plot, show_plot)

#	if write_to_file == 1
#	"imprimir os arquivos tan_cc_x, tan_cc_y, tan_caustic_x, tan_caustic_y in to a file tang_curves.txt"
#	"imprimir os arquivos rad_cc_x, rad_cc_y, rad_caustic_x, rad_caustic_y in to a file rad_curves.txt"







