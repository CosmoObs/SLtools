#saidas das coordenadas das CC


import os
import logging
from lens_parameters_new import lens_parameters_new
#from sltools.gravlens.lens_parameters_new import lens_parameters_new

## get_source_model
# this function executes gravlens to get the file with the critical curves ("crit.txt")
#
#@param inputlens (the base text a input file of gravlens must have, "set gridhi1 = #", "setlens 1 1", etc)
#@return File "crit.txt" with the critical curves positions
def critcurves(gravinput, filename):
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


## find_CC_new
# this function iteratively runs "critcurves" (iteranting on gridhi1) to find the critical curves
#
#@param lens_model
#@param gravlens_params
#@return "inputlens" and "setlens": two strins that are used to compose gravlens' input files
def find_CC_new(lens, mass_scale, model_param_8, model_param_9, model_param_10, galaxy_position=[0,0], e_L=0, theta_L=0, shear=0, theta_shear=0, gravlens_params={}, filename='crit.txt'):

	inputlens, setlens, gravlens_params_updated = lens_parameters_new(lens, mass_scale, model_param_8, model_param_9, model_param_10, galaxy_position, e_L, theta_L, shear, theta_shear, gravlens_params) # inputlens is the gravlens input (sets general parameters and the lens parameters)	# setlens is the gravlens input line that concerns the lens (ex: nfw 1 0 ...)
	logging.debug('Determined the strings that defines the mass model in gravlens through the \'lensparameter\' function')
	#--------------------------------------

	#-----------------------------
	critcurves(inputlens, filename) # gets the critical curves (crit.txt file)
	logging.debug('Got the critical curves (%s file) with function \"critcurves\"' % filename)
	#-----------------------------

	# ITERATION ON gridhi1 (TO GET BETTER RESOLUTION ON CC)
	counter = 0
	while os.path.isfile('./' + filename) == False and counter < 20: # looks for the C.C. file (crit.txt) - if this file does not exist, add a note to arccatalog and quit
		gravlens_params_updated['gridhi1'] = float( gravlens_params_updated['gridhi1'] ) / 5.; # gridhi1 /= 5. 

		lens_par_out = lens_parameters_new(lens, mass_scale, model_param_8, model_param_9, model_param_10, galaxy_position, e_L, theta_L, shear, theta_shear, gravlens_params_updated) # inputlens is the gravlens input (sets general parameters and the lens parameters)

		inputlens, setlens, gravlens_params_updated = lens_par_out
		critcurves(inputlens, filename) # gets the critical curves (crit.txt file)
		counter += 1
	logging.debug( '(Number of iterations on gridhi1 = %d)' % (counter) )

	if os.path.isfile('./' + filename) == False: # looks for the C.C. file (crit.txt) - if this file does not exist, returns 'False'
		return False

 	if len(open('./' + filename).readlines() ) < 200: # these correspond to critical cases when the CC have too few points.
		gravlens_params_updated['gridhi1'] = float( gravlens_params_updated['gridhi1'] ) / 3. # gridhi1 /= 3. 

		lens_par_out = lens_parameters_new(lens, mass_scale, model_param_8, model_param_9, model_param_10, galaxy_position, e_L, theta_L, shear, theta_shear, gravlens_params_updated) # inputlens is the gravlens input (sets general parameters and the lens parameters)
		inputlens, setlens, gravlens_params_updated = lens_par_out

		critcurves(inputlens, filename) # gets the critical curves (crit.txt file)
	logging.debug( 'gridhi1 = %f and number of iterations on gridhi1 = %d' % (float( gravlens_params_updated['gridhi1'] ),  counter) )


	caustics_xy, crit_curves_xy = get_CC_from_file(filename)

	return caustics_xy, crit_curves_xy
