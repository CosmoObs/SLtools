# Output of coordinates of the CC curves (lens curves)
# By Pedro F.
# Nov 2010
# Comment cleanup by MSSG

##@package find_CC_new
# Determines the caustics and CC of a given lens
#
# Determines the lens curves using a iterating method using gravlens
#

import os
import logging
#from lens_parameters_new import lens_parameters_new
from sltools.gravlens.lens_parameters_new import lens_parameters_new
import numpy as np

# given the gravlens 'config file' runs gravlens to generate the 
# file with caustics and CC
#
def _critcurves_new(gravinput, caustic_CC_file, gravlens_input_file='gravlens_CC_input.txt'):
	fmag = open(gravlens_input_file, 'w')
	fmag.write(gravinput)
	fmag.write('plotcrit %s\n' % caustic_CC_file )
	fmag.close()
	os.system('gravlens %s > /dev/null' % gravlens_input_file )

## get_CC_from_file
# Reads critical curves file outputed by gravlens (with plotmode=1)
#
def get_CC_from_file(caustic_CC_file):
	"""
	Reads critical curves file outputed by gravlens (with 
	plotmode=1). 
	
	Input:
	 - caustic_CC_file   <str> : name of the file to be read

	Output:
	 - <list>     : [x_caustic, ycaustic], with x_caustic and ycaustic being the lists with the caustic coordinates
	 - <list>     : [x_CC, y_CC], with x_CC and y_CC being the lists with the CC coordinates

	"""
	fcrit = open(caustic_CC_file, 'r').readlines()
	xcritsrc = []
	ycritsrc = []
	xcritimg = []
	ycritimg = []
	for i in range (10,len(fcrit)): # because the first 10 lines of this file are only comments
		xcritimg.append(float(fcrit[i].split()[0]))
		ycritimg.append(float(fcrit[i].split()[1]))
		xcritsrc.append(float(fcrit[i].split()[2]))
		ycritsrc.append(float(fcrit[i].split()[3]))

	return [xcritsrc, ycritsrc], [xcritimg, ycritimg]

## find_CC_new
# Determines the caustics and CC of a given lens
#
# Determines the lens curves using a iterating method using gravlens
#
def find_CC_new(lens_model, mass_scale, model_param_8, model_param_9, model_param_10, galaxy_position=[0,0], e_L=0, theta_L=0, shear=0, theta_shear=0, gravlens_params={}, caustic_CC_file='crit.txt', gravlens_input_file='gravlens_CC_input.txt'):
	""" 
	Determines the caustics and critical curves for a single 
	component lens model. Uses an iterative procedure to determine 
	the grid size and the minimum number of points on the curves.

	The CC size may vary considerably depending on the lens-source
	configuration. Because of this, we developed an iterative
	method to find the CC, that involves the size of the grid used
	by gravlens (gridhi1). The initial value of gridhi1 must be an
	upper limit for the size of the CC, and find_CC_new will try
	to find the CC. Each time gravlens can't find the CC, we lower
	gridhi1 by a factor of grid_factor (=5), in order to increase 
	precision. This continues on until gravlens finds the CC or 
	the number of max_iterations_number (=20) iterations is reached. 
	After finding the CC, if it is composed of less than min_n_lines
	(=200) points, the code redefines gridhi1 to a third of it 
	(grid_factor2=3), usually increasing the number of points. With 
	this first determination of the CC, the code uses its scale 
	(gridhi1_CC_factor(=2) times the distance of the furthest CC 
	point to the origin) to reobtain the CC with the apropriate 
	value for the grid.

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

	"""

	# definitions of some control params
	grid_factor = 5.
	max_iterations_number = 20
	grid_factor2 = 3.
	min_n_lines = 200
	gridhi1_CC_factor = 2.
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
		return False

 	if len(open('./' + caustic_CC_file).readlines() ) < min_n_lines: # these correspond to critical cases when the CC have too few points.
		gravlens_params_updated['gridhi1'] = float( gravlens_params_updated['gridhi1'] ) / grid_factor2 # gridhi1 /= 3. 

		lens_par_out = lens_parameters_new(lens_model, mass_scale, model_param_8, model_param_9, model_param_10, galaxy_position, e_L, theta_L, shear, theta_shear, gravlens_params_updated) # inputlens is the gravlens input (sets general parameters and the lens parameters)
		inputlens, setlens, gravlens_params_updated = lens_par_out

		_critcurves_new(inputlens, caustic_CC_file, gravlens_input_file) # gets the critical curves (crit.txt file)
	logging.debug( 'gridhi1 = %f and number of iterations on gridhi1 = %d' % (float( gravlens_params_updated['gridhi1'] ),  counter) )


	caustics_xy, crit_curves_xy = get_CC_from_file(caustic_CC_file)

	#-----------------------------------------------------------------------------------------------------------
	# redefine gridhi1 according to the CC size
	tan_CC_x = np.array(crit_curves_xy[0])
	tan_CC_y = np.array(crit_curves_xy[1])
	index_CC = np.argmax(tan_CC_x**2 + tan_CC_y**2)

	gravlens_params_updated['gridhi1'] =  grid_factor * ( (tan_CC_x[index_CC]**2 + tan_CC_y[index_CC]**2)**0.5 )
	#-----------------------------------------------------------------------------------------------------------
	# get CC with best value for the grid
	lens_par_out = lens_parameters_new(lens_model, mass_scale, model_param_8, model_param_9, model_param_10, galaxy_position, e_L, theta_L, shear, theta_shear, gravlens_params_updated) # inputlens is the gravlens input (sets general parameters and the lens parameters)
	inputlens, setlens, gravlens_params_updated = lens_par_out

	_critcurves_new(inputlens, caustic_CC_file, gravlens_input_file) # gets the critical curves (crit.txt file)

	caustics_xy, crit_curves_xy = get_CC_from_file(caustic_CC_file)

	return caustics_xy, crit_curves_xy, gravlens_params_updated
