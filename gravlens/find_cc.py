#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ==================================
# Authors:
# Pedro Ferreira - pferreira@dfte.ufrn.br
# Comments cleanup by MSSG
# Habib Dumet Montoya - habibdumet@gmail.com
# ==================================

""" Determines, separates and plots the caustics and the critical curves of a given lens model """


##@package find_cc
# 
#
# This package has functions to i) get the caustic and the critical curves of a given model (see find_cc) and 
# ii) plot the caustics and the critical curves (plot_CC).
# The curves are determined using an iterating method that uses gravlens (see find_cc).
# The curves are identifyed as radial and tangencial using the function sltools.geometry.separate_curves.

import os
import logging
import numpy as np
import matplotlib.pyplot as pyplot
import string
import time

from sltools.gravlens.lens_parameters import lens_parameters
from sltools.geometry.separate_curves import separate_curves
from sltools.coordinate.translation_and_rotation import translate_and_rotate_coord_system



def critcurves(gravinput, caustic_CC_file, gravlens_input_file='gravlens_CC_input.txt'):
	"""
	Given the gravlens 'config file', runs gravlens to generate the file with caustics and critical 
	curves.

	Input:
	 - gravinput           <str> : a string with all the lines gravlens needs to its configuration
	 - caustic_CC_file     <str> : name of the output file with the caustic and CC positions
	 - gravlens_input_file <str> : name of the input file used to run gravlens

	Output:
	 - <str> : name of the output file with the caustic and CC positions (caustic_CC_file)

	"""
	
	fmag = open(gravlens_input_file, 'w')
	fmag.write(gravinput)
	fmag.write('plotcrit %s\n' % caustic_CC_file )
	fmag.close()
	os.system('gravlens %s > /dev/null' % gravlens_input_file )

	return caustic_CC_file


#=======================================================================================================
def find_cc(lens_model, mass_scale, model_param_8, model_param_9, model_param_10, galaxy_position=[0,0], 
		e_L=0, theta_L=0, shear=0, theta_shear=0, gravlens_params={}, caustic_CC_file='crit.txt', 
		gravlens_input_file='gravlens_CC_input.txt', grid_factor=5., grid_factor2=3., max_iter_number=20,
		min_n_lines=200, gridhi1_CC_factor=2., accept_res_limit=2E-4):
	""" 
	Determines the caustics and critical curves (CC) for a single component lens model. Uses an iterative
	procedure to determine the grid size and the minimum number of points on the curves.

	The CC size may vary considerably depending on the lens-source configuration. Because of this, 
	we developed an iterative method to find the CC, that involves the size of the grid used by 
	gravlens (gridhi1). The initial value of gridhi1 must be an upper limit for the size of the CC, 
	and find_cc will try to find the CC. Each time gravlens can't find the CC, we lower	gridhi1 
	by a factor of grid_factor (=5), in order to increase precision. This continues on until 
	gravlens finds the CC or the number of iterations, max_iter_number (=20), is reached. After 
	finding the CC, if the file is composed of less than min_n_lines (=200) lines, the code redefines 
	gridhi1 to gridhi1/grid_factor2 (grid_factor2=3), usually increasing the number of points. With this 
	first determination of the CC, the code uses its scale (gridhi1_CC_factor(=2) times the distance
	of the furthest CC point to the lens center) to reobtain the CC with the apropriate value for the 
	grid, ensuring best precision. If the CC were not found after these attempts, this function returns "False" and an 
	error message. If all caustic points have coordinates < accept_res_limit (=2E-4), a warning 
	message will be generated, meaning that the precision of gravlens (limited by the number of 
	digits in the output file) is being reached. 
 
	Example:
	         x_caustic, y_caustic, x_CC, y_CC, gravlens_config = find_cc('nfw',1,1,0,0)

	Input:
	 - lens_model          <str> : lens name (see gravlens manual table 3.1)
	 - mass_scale        <float> : mass scale of the lens - "parameter 1"
	 - model_param_8     <float> : misc. lens parameter - often scale radio (depends on the lens model)
	 - model_param_9     <float> : misc. lens parameter - often scale radio (depends on the lens model)
	 - model_param_10    <float> : misc. lens parameter - often a power law index  (depends on the lens model)
	 - galaxy_position    <list> : [x,y] position of the lens
	 - e_L               <float> : lens ellipticity (default=0)
	 - theta_L           <float> : lens position angle (in degrees) with respect to the vertical 
				       (counterclockwise)
	 - shear             <float> : external shear amplitude
	 - theta_shear       <float> : external shear direction (in degrees)
	 - gravlens_params     <dic> : contains the keys and values of the gravlens configuration (see 
				       default parameters at function set_gravlens_default,inside lens_parameters)
	 - caustic_CC_file     <str> : name of the output file with the caustic and CC positions
	 - gravlens_input_file <str> : name of the input file used to run gravlens
	 - grid_factor       <float> : used in the iteration to find the critical curves. Each time no 
				       curves were found, the grid size is divided by grid_factor
	 - min_n_lines         <int> : minimum number of lines demanded in the CC file to decrease the grid by grid_factor_2 
	 - grid_factor2      <float> : if the CC file is composed of less than min_n_lines, the code divides
				       the grid size by grid_factor2 
	 - max_iter_number     <int> : maximum number of time the grid size will be divided by grid_factor
	 - gridhi1_CC_factor <float> : the final size of the grid will be gridhi1_CC_factor * (the furthest CC point)
	 - accept_res_limit  <float> : if the furthest CC point is closer than accept_res_limit from the lens 
				       center, display a warning about lack of precision in determining the curves

	Output:
	 - <list>     : [x_caustic, ycaustic], with x_caustic and ycaustic being the lists with the 
			caustic coordinates
	 - <list>     : [x_CC, y_CC], with x_CC and y_CC being the lists with the CC coordinates
	 - <dict>     : all configuration variables used for running gravlens (including gridhi1)
	 - <str>      : name of the file with the caustic and CC positions ('caustic_CC_file')
	 - <str>      : name of the input file used to run gravlens (gravlens_input_file)

	"""


	inputlens, setlens, gravlens_params_updated = lens_parameters(lens_model, mass_scale, model_param_8, model_param_9, model_param_10, galaxy_position, e_L, theta_L, shear, theta_shear, gravlens_params) # inputlens is the gravlens input (sets general parameters and the lens parameters)	# setlens is the gravlens input line that concerns the lens (ex: nfw 1 0 ...)
	logging.debug('Determined the strings that defines the mass model in gravlens through the \'lens_parameters\' function')
	#--------------------------------------

	#-----------------------------
	os.system('rm -f %s' % caustic_CC_file) # if the file previously exists, delete it
	critcurves(inputlens, caustic_CC_file, gravlens_input_file) # gets the critical curves (crit.txt file)
	logging.debug('Got the critical curves (%s file) with function \"critcurves\"' % caustic_CC_file)
	#-----------------------------

	# ITERATION ON gridhi1 (TO GET BETTER RESOLUTION ON CC)
	counter = 0
	while os.path.isfile('./' + caustic_CC_file) == False and counter < max_iter_number: # looks for the C.C. file (crit.txt) - if this file does not exist, add a note to arccatalog and quit
		gravlens_params_updated['gridhi1'] = float( gravlens_params_updated['gridhi1'] ) / grid_factor; # gridhi1 /= 5. 

		lens_par_out = lens_parameters(lens_model, mass_scale, model_param_8, model_param_9, model_param_10, galaxy_position, e_L, theta_L, shear, theta_shear, gravlens_params_updated) # inputlens is the gravlens input (sets general parameters and the lens parameters)

		inputlens, setlens, gravlens_params_updated = lens_par_out
		critcurves(inputlens, caustic_CC_file, gravlens_input_file) # gets the critical curves (crit.txt file)
		counter += 1
	logging.debug( '(Number of iterations on gridhi1 = %d)' % (counter) )

	if os.path.isfile('./' + caustic_CC_file) == False: # looks for the C.C. file (crit.txt) - if this file does not exist, returns 'False'
		logging.error( 'Gravlens could not determine the caustics and critical curves. Try changing the input units (for example, from arcseconds to miliarcsecond), because gravlens outputs the critical curves with only 6 digits. The grid size used was %f and it was obtained through %d iteration.' % (float( gravlens_params_updated['gridhi1'] ),  counter) )
		return False

 	if len(open('./' + caustic_CC_file).readlines() ) < min_n_lines: # these correspond to critical cases when the CC have too few points.
		gravlens_params_updated['gridhi1'] = float( gravlens_params_updated['gridhi1'] ) / grid_factor2 # gridhi1 /= 3. 

		lens_par_out = lens_parameters(lens_model, mass_scale, model_param_8, model_param_9, model_param_10, galaxy_position, e_L, theta_L, shear, theta_shear, gravlens_params_updated) # inputlens is the gravlens input (sets general parameters and the lens parameters)
		inputlens, setlens, gravlens_params_updated = lens_par_out

		critcurves(inputlens, caustic_CC_file, gravlens_input_file) # gets the critical curves (crit.txt file)


	x1, y1, u1, v1 = np.loadtxt(caustic_CC_file, usecols = (0,1,2,3), unpack=True) # CC_x, CC_y, caustic_x, caustic_y

	#-----------------------------------------------------------------------------------------------------------
	# redefine gridhi1 according to the CC size

	index_CC = np.argmax(x1**2 + y1**2) # it is not necessary to use the lens center since we want the furthest point anyway
	logging.debug( 'The CC furthest point is at (%s, %s), a distance %s of the origin .' % (x1[index_CC], y1[index_CC], (x1[index_CC]**2 + y1[index_CC]**2)**0.5 ) )

	gravlens_params_updated['gridhi1'] =  gridhi1_CC_factor * ( (x1[index_CC]**2 + y1[index_CC]**2)**0.5 )

	#-----------------------------------------------------------------------------------------------------------
	# get CC with best value for the grid
	lens_par_out = lens_parameters(lens_model, mass_scale, model_param_8, model_param_9, model_param_10, galaxy_position, e_L, theta_L, shear, theta_shear, gravlens_params_updated) # inputlens is the gravlens input (sets general parameters and the lens parameters)
	inputlens, setlens, gravlens_params_updated = lens_par_out

	logging.debug( 'gridhi1 = %f and number of iterations on gridhi1 = %d' % (float( gravlens_params_updated['gridhi1'] ),  counter) )

	critcurves(inputlens, caustic_CC_file, gravlens_input_file) # gets the critical curves (crit.txt file)

	x1, y1, u1, v1 = np.loadtxt(caustic_CC_file, usecols = (0,1,2,3), unpack=True) # CC_x, CC_y, caustic_x, caustic_y

#	index_CC = np.argmax(x1**2 + y1**2)
#	logging.debug( 'The CC furthest point is at (%s, %s), a distance %s of the origin .' % (x1[index_CC], y1[index_CC], (x1[index_CC]**2 + y1[index_CC]**2)**0.5 ) )

	# check if the precision is ok (gravlens outputs coordinates with only 6 decimal places)
	if max( max(np.abs(u1)), max(np.abs(v1)) ) < accept_res_limit:
		logging.warning('The caustics seem to be determined with low resolution')

	return x1, y1, u1, v1, gravlens_params_updated




#=======================================================================================================
def plot_cc(tan_caustic_x, tan_caustic_y, rad_caustic_x, rad_caustic_y, tan_CC_x, tan_CC_y, rad_CC_x, rad_CC_y, plot_filename, show_plot=0, src_x=[], src_y=[], img_x=[], img_y=[] ):
	"""
	Plots given caustics and critical curves (CC).

	Plots caustics and CC side by side. It can be chosen if the plot will be displayed in the screen
	instead of saved to a file.

	Input:
	 - tan_caustic_x  <list> : list of x coordinates of points from the tangencial caustic
         - tan_caustic_y  <list> : list of y coordinates of points from the tangencial caustic
         - rad_caustic_x  <list> : list of x coordinates of points from the radial caustic
         - rad_caustic_y  <list> : list of y coordinates of points from the radial caustic
         - tan_CC_x       <list> : list of x coordinates of points from the tangencial CC
         - tan_CC_y       <list> : list of y coordinates of points from the tangencial CC
         - rad_CC_x       <list> : list of x coordinates of points from the radial CC
         - rad_CC_y       <list> : list of y coordinates of points from the radial CC
	 - plot_filename   <str> : the name of the generated plot file (including the format, ex. 
			           'curves_plot.png'). To see the formats available see matplotlib.pyplot help
	 - show_plot       <int> : use 0 for 'no screen display' (default) and 1 for display on the 
				   screen

	Output:
	 - <str> : name of the plot file (plot_filename). It will be generated only if show_plot=0 (default)

	"""

	pyplot.clf()
	f1 = pyplot.figure(figsize = (20,12) )
	pyplot.subplot(121) # (numRows,numCols, plotNum)
	pyplot.plot(tan_caustic_x, tan_caustic_y, color='black', marker='', linewidth=2, label='Tangential')
	pyplot.plot(rad_caustic_x, rad_caustic_y, color='red', marker='', linewidth=2, label='Radial')
	if len(src_x) != 0:
		pyplot.plot(src_x, src_y, 'b.', label='Source')
	pyplot.legend(loc='upper right', shadow=True)
	pyplot.axis('equal')
	pyplot.xlabel('x', fontsize = 15)
	pyplot.ylabel('y', fontsize = 15)
	pyplot.title('Caustics'  , fontsize = 20)

	pyplot.subplot(122) # (numRows,numCols, plotNum)
	pyplot.plot(tan_CC_x, tan_CC_y, color='black', marker='', linewidth=2, label='Tangential' )
	pyplot.plot(rad_CC_x, rad_CC_y, color='red', marker='', linewidth=2, label='Radial' )
	if len(img_x) != 0:
		pyplot.plot(img_x, img_y, 'g.', label='Images')
	pyplot.legend(loc='upper right', shadow=True)	
	pyplot.axis('equal')
	pyplot.xlabel('x', fontsize = 15)
	pyplot.ylabel('y', fontsize = 15)
	pyplot.title('Critical Curves'  , fontsize = 20)
	if show_plot == 0:
		pyplot.savefig(plot_filename)
	if show_plot == 1:
		pyplot.show() 

        pyplot.close()

	return True


