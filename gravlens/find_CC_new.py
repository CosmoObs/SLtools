# Output of coordinates of the CC curves (lens curves)
# By Pedro F.
# Nov 2010
# Comment cleanup by MSSG

import os
import logging
#from lens_parameters_new import lens_parameters_new
from sltools.gravlens.lens_parameters_new import lens_parameters_new


def critcurves_new(gravinput, filename):
	fmag = open('gravlens_CC_input.txt', 'w')
	fmag.write(gravinput)
	fmag.write('plotcrit %s\n' % filename )
	fmag.close()
	os.system('gravlens gravlens_CC_input.txt > /dev/null')


def get_CC_from_file(filename):
	fcrit = open(filename, 'r').readlines()
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



def find_CC_new(lens, mass_scale, model_param_8, model_param_9, model_param_10, galaxy_position=[0,0], e_L=0, theta_L=0, shear=0, theta_shear=0, gravlens_params={}, filename='crit.txt'):
	""" Determines the caustics and critical curves of a given lens model

	The CC size may vary considerably depending on the lens-source
	configuration. Because of this, we developed an iterative
	method to find the CC, that involves the size of the grid used
	by gravlens (gridhi1). The initial value of gridhi1 must be an
	upper limit for the size of the CC, and find_CC_new will try
	to find the CC. Each time gravlens can't find the CC, we lower
	gridhi1 by a factor of 5, in order to increase precision. This
	continues on until gravlens finds the CC or the number of 20
	iterations is reached. After finding the CC, if it is composed
	of less than 200 points, the code redefines gridhi1 to a third
	of it, usually increasing the number of points.

	Input:
	 - lens                   <str> : Lens name (see gravlens manual table 3.1)
	 - mass_scale           <float> : Mass scale of the lens - "parameter 1"
	 - model_param_8        <float> : Vertical size (in pixels) of output image
	 - model_param_9        <float> : Vertical size (in pixels) of output image
	 - model_param_10       <float> : Vertical size (in pixels) of output image
	 - galaxy_position <list float> : [x,y] position of the lens
	 - e_L                  <float> : Horizontal central position for output (cut) image
	 - theta_L              <float> : Vertical central position for output (cut) image
	 - shear                <float> : 'pixel' or 'degrees' for size (x_size,y_size) values
	 - theta_shear          <float> : Horizontal size (in pixels) of output image
	 - gravlens_params        <dic> : Contains the keys and values of the gravlens configuration (see default 
				          parameters at function set_gravlens_default (inside lens_parameters_new)
	 - filename               <str> : name of the output file with the castic and CC positions

	Output:
	 - <list>     : [x_caustic, ycaustic], with x_caustic and ycaustic being the lists with the caustic coordinates
	 - <list>     : [x_CC, y_CC], with x_CC and y_CC being the lists with the CC coordinates
	 - <file>     : file named 'filename' with the caustic and CC positions

	"""


	inputlens, setlens, gravlens_params_updated = lens_parameters_new(lens, mass_scale, model_param_8, model_param_9, model_param_10, galaxy_position, e_L, theta_L, shear, theta_shear, gravlens_params) # inputlens is the gravlens input (sets general parameters and the lens parameters)	# setlens is the gravlens input line that concerns the lens (ex: nfw 1 0 ...)
	logging.debug('Determined the strings that defines the mass model in gravlens through the \'lensparameter\' function')
	#--------------------------------------

	#-----------------------------
	critcurves_new(inputlens, filename) # gets the critical curves (crit.txt file)
	logging.debug('Got the critical curves (%s file) with function \"critcurves\"' % filename)
	#-----------------------------

	# ITERATION ON gridhi1 (TO GET BETTER RESOLUTION ON CC)
	counter = 0
	while os.path.isfile('./' + filename) == False and counter < 20: # looks for the C.C. file (crit.txt) - if this file does not exist, add a note to arccatalog and quit
		gravlens_params_updated['gridhi1'] = float( gravlens_params_updated['gridhi1'] ) / 5.; # gridhi1 /= 5. 

		lens_par_out = lens_parameters_new(lens, mass_scale, model_param_8, model_param_9, model_param_10, galaxy_position, e_L, theta_L, shear, theta_shear, gravlens_params_updated) # inputlens is the gravlens input (sets general parameters and the lens parameters)

		inputlens, setlens, gravlens_params_updated = lens_par_out
		critcurves_new(inputlens, filename) # gets the critical curves (crit.txt file)
		counter += 1
	logging.debug( '(Number of iterations on gridhi1 = %d)' % (counter) )

	if os.path.isfile('./' + filename) == False: # looks for the C.C. file (crit.txt) - if this file does not exist, returns 'False'
		return False

 	if len(open('./' + filename).readlines() ) < 200: # these correspond to critical cases when the CC have too few points.
		gravlens_params_updated['gridhi1'] = float( gravlens_params_updated['gridhi1'] ) / 3. # gridhi1 /= 3. 

		lens_par_out = lens_parameters_new(lens, mass_scale, model_param_8, model_param_9, model_param_10, galaxy_position, e_L, theta_L, shear, theta_shear, gravlens_params_updated) # inputlens is the gravlens input (sets general parameters and the lens parameters)
		inputlens, setlens, gravlens_params_updated = lens_par_out

		critcurves_new(inputlens, filename) # gets the critical curves (crit.txt file)
	logging.debug( 'gridhi1 = %f and number of iterations on gridhi1 = %d' % (float( gravlens_params_updated['gridhi1'] ),  counter) )


	caustics_xy, crit_curves_xy = get_CC_from_file(filename)

	return caustics_xy, crit_curves_xy
