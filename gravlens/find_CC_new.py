#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ==================================
# Authors:
# Pedro Ferreira - pferreira@dfte.ufrn.br
# Comments cleanup by MSSG
# Habib Dumet Montoya - habibdumet@gmail.com
# ==================================

""" Determines, separates and plots the caustics and the critical curves of a given lens model """


##@package find_CC_new
# 
#
# This package has functions to i) get the caustic and the critical curves of a given model (see find_CC_new) and 
# ii) plot the caustics and the critical curves (plot_CC).
# The curves are determined using an iterating method that uses gravlens (see find_CC_new).
# The curves are identifyed as radial and tangencial using the function sltools.geometry.separate_curves.
# There is a 3rd function that concatenates find_CC_new, separate_curvs and plot_CC, optimizing the curves 
# separation (see run_find_CC for details)
# For high-level documentation, see http://twiki.linea.gov.br/bin/view/StrongLensing/SLtoolsGravlens#Module_find_CC_new_py

import os
import logging
import numpy as np
import matplotlib.pyplot as pyplot
import string

from sltools.gravlens.lens_parameters_new import lens_parameters_new
from sltools.geometry.separate_curves import separate_curves_a
from sltools.coordinate.translation_followed_by_rotation import translation_and_rotation



def _critcurves_new(gravinput, caustic_CC_file, gravlens_input_file='gravlens_CC_input.txt'):
	"""
	Given the gravlens 'config file', runs gravlens to generate the file with caustics and critical 
	curves.

	"""
	
	fmag = open(gravlens_input_file, 'w')
	fmag.write(gravinput)
	fmag.write('plotcrit %s\n' % caustic_CC_file )
	fmag.close()
	os.system('gravlens %s > /dev/null' % gravlens_input_file )



#=======================================================================================================
def find_CC_new(lens_model, mass_scale, model_param_8, model_param_9, model_param_10, galaxy_position=[0,0], 
		e_L=0, theta_L=0, shear=0, theta_shear=0, gravlens_params={}, caustic_CC_file='crit.txt', 
		gravlens_input_file='gravlens_CC_input.txt', grid_factor=5., grid_factor2=3., max_iter_number=20,
		min_n_lines=200, gridhi1_CC_factor=2., accept_res_limit=2E-4):
	""" 
	Determines the caustics and critical curves (CC) for a single component lens model. Uses an iterative
	procedure to determine the grid size and the minimum number of points on the curves.

	The CC size may vary considerably depending on the lens-source configuration. Because of this, 
	we developed an iterative method to find the CC, that involves the size of the grid used by 
	gravlens (gridhi1). The initial value of gridhi1 must be an upper limit for the size of the CC, 
	and find_CC_new will try to find the CC. Each time gravlens can't find the CC, we lower	gridhi1 
	by a factor of grid_factor (=5), in order to increase precision. This continues on until 
	gravlens finds the CC or the number of iterations, max_iter_number (=20), is reached. After 
	finding the CC, if it is composed of less than min_n_lines (=200) points, the code redefines 
	gridhi1 to a third of it (grid_factor2=3), usually increasing the number of points. With this 
	first determination of the CC, the code uses its scale (gridhi1_CC_factor(=2) times the distance
	of the furthest CC point to the origin) to reobtain the CC with the apropriate value for the 
	grid. If the CC were not found after all these attempts, this function returns "False" and an 
	error message. If all caustic points have coordinates < accept_res_limit (=2E-4), a warning 
	message will be generated, meaning that the precision of gravlens (limited by the number of 
	digits in the output file) is being reached. 
 
	Example:
	         x_caustic, y_caustic, x_CC, y_CC, gravlens_config = find_CC_new('nfw',1,1,0,0)

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
				       default parameters at function set_gravlens_default,inside lens_parameters_new)
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


	inputlens, setlens, gravlens_params_updated = lens_parameters_new(lens_model, mass_scale, model_param_8, model_param_9, model_param_10, galaxy_position, e_L, theta_L, shear, theta_shear, gravlens_params) # inputlens is the gravlens input (sets general parameters and the lens parameters)	# setlens is the gravlens input line that concerns the lens (ex: nfw 1 0 ...)
	logging.debug('Determined the strings that defines the mass model in gravlens through the \'lens_parameters\' function')
	#--------------------------------------

	#-----------------------------
	os.system('rm -f %s' % caustic_CC_file) # if the file previously exists, delete it
	_critcurves_new(inputlens, caustic_CC_file, gravlens_input_file) # gets the critical curves (crit.txt file)
	logging.debug('Got the critical curves (%s file) with function \"critcurves\"' % caustic_CC_file)
	#-----------------------------

	# ITERATION ON gridhi1 (TO GET BETTER RESOLUTION ON CC)
	counter = 0
	while os.path.isfile('./' + caustic_CC_file) == False and counter < max_iter_number: # looks for the C.C. file (crit.txt) - if this file does not exist, add a note to arccatalog and quit
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


	x1, y1, u1, v1 = np.loadtxt(caustic_CC_file, usecols = (0,1,2,3), unpack=True) # CC_x, CC_y, caustic_x, caustic_y

	#-----------------------------------------------------------------------------------------------------------
	# redefine gridhi1 according to the CC size

	index_CC = np.argmax(x1**2 + y1**2) # it is not necessary to use the lens center since we want the furthest point anyway
	logging.debug( 'The CC furthest point is at (%s, %s), a distance %s of the origin .' % (x1[index_CC], y1[index_CC], (x1[index_CC]**2 + y1[index_CC]**2)**0.5 ) )

	gravlens_params_updated['gridhi1'] =  gridhi1_CC_factor * ( (x1[index_CC]**2 + y1[index_CC]**2)**0.5 )

	#-----------------------------------------------------------------------------------------------------------
	# get CC with best value for the grid
	lens_par_out = lens_parameters_new(lens_model, mass_scale, model_param_8, model_param_9, model_param_10, galaxy_position, e_L, theta_L, shear, theta_shear, gravlens_params_updated) # inputlens is the gravlens input (sets general parameters and the lens parameters)
	inputlens, setlens, gravlens_params_updated = lens_par_out

	logging.debug( 'gridhi1 = %f and number of iterations on gridhi1 = %d' % (float( gravlens_params_updated['gridhi1'] ),  counter) )

	_critcurves_new(inputlens, caustic_CC_file, gravlens_input_file) # gets the critical curves (crit.txt file)

	x1, y1, u1, v1 = np.loadtxt(caustic_CC_file, usecols = (0,1,2,3), unpack=True) # CC_x, CC_y, caustic_x, caustic_y

#	index_CC = np.argmax(x1**2 + y1**2)
#	logging.debug( 'The CC furthest point is at (%s, %s), a distance %s of the origin .' % (x1[index_CC], y1[index_CC], (x1[index_CC]**2 + y1[index_CC]**2)**0.5 ) )

	# check if the precision is ok (gravlens outputs coordinates with only 6 decimal places)
	if max( max(np.abs(u1)), max(np.abs(v1)) ) < accept_res_limit:
		logging.warning('The caustics seem to be determined with low resolution')

	return x1, y1, u1, v1, gravlens_params_updated




#=======================================================================================================
def plot_CC(tan_caustic_x, tan_caustic_y, rad_caustic_x, rad_caustic_y, tan_CC_x, tan_CC_y, rad_CC_x, rad_CC_y, plot_filename, show_plot=0):
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
	 - show_plot       <int> : use 0 for 'no screen display' (default) and 1 for diaplay on the 
				   screen

	Output:
	 - <str> : name of the plot file (plot_filename). It will be generated only if show_plot=0 (default)

	"""

	pyplot.clf()
	f1 = pyplot.figure(figsize = (20,12) )
	pyplot.subplot(121) # (numRows,numCols, plotNum)
	pyplot.plot(tan_caustic_x, tan_caustic_y, color='black',marker='', linewidth=2, label='Tangential')
	pyplot.plot(rad_caustic_x, rad_caustic_y, color='red', marker='', linewidth=2, label='Radial')
	pyplot.legend(loc='upper right', shadow=True)
	pyplot.axis('equal')
	pyplot.xlabel('x', fontsize = 15)
	pyplot.ylabel('y', fontsize = 15)
	pyplot.title('Caustics'  , fontsize = 20)

	pyplot.subplot(122) # (numRows,numCols, plotNum)
	pyplot.plot(tan_CC_x, tan_CC_y, color='black',marker='', linewidth=2, label='Tangential' )
	pyplot.plot(rad_CC_x, rad_CC_y, color='red', marker='', linewidth=2, label='Radial' )
	pyplot.legend(loc='upper right', shadow=True)	
	pyplot.axis('equal')
	pyplot.xlabel('x', fontsize = 15)
	pyplot.ylabel('y', fontsize = 15)
	pyplot.title('Critical Curves'  , fontsize = 20)
	if show_plot == 0:
		pyplot.savefig(plot_filename)
	if show_plot == 1:
		pyplot.show() 

	return True



## run_find_CC_new
# Finds the caustics and CC for a given lens model, separates the radial from the tangential and plots the curves
#
def run_find_CC(lens_model, mass_scale, model_param_8, model_param_9, model_param_10, galaxy_position=[0.,0.], e_L=0, theta_L=0, shear=0, theta_shear=0, gravlens_params={}, caustic_CC_file='crit.txt', gravlens_input_file='gravlens_CC_input.txt', rad_curves_file='lens_curves_rad.dat', tan_curves_file='lens_curves_tan.dat', curves_plot='crit-caust_curves.png', show_plot=0, write_to_file=0, max_delta_count=20, delta_increment=1.1, grid_factor=5., grid_factor2=3., max_iter_number=20, min_n_lines=200, gridhi1_CC_factor=2., accept_res_limit=2E-4):
	""" 
	This is a pipeline that runs 'find_CC_new', 'separate_curves' and 'plot_CC'. For details of these 
	functions, see their documentation.

	As we expect 2 critical curves, if the function separate_curves finds a different number, we 
	question it. If we find only one CC, then the separation code did not worked properly, and we 
	noticed that for some values of theta_L this occurs quite often. So we find the CC again for
	theta_L=0, separate them and rotate them back to its original angle, which proved to be a good
	solution.
	In the case that more than 2 CC are found, we try an iterating process with the control parameter
	of separate_curves 'delta' (see separate_curves documentation). So we increase the value of delta 
	until there are 2 CC.
	The combination of these 2 implementations above with the default gravlens parameters used in 
	'lens_parameters_new' allows for a 100% correct separation in a NFW profile with all possible 
	combination of parameters kappas = [0.03, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.87], 
	rs = [34., 54., 74., 94.] and theta_L = [0, 45, 90, 135, 180]. 
	For more details, see http://twiki.linea.gov.br/bin/view/StrongLensing/SLtoolsGravlens#Module_find_CC_new_py.

	We separate only the image plane curves (CC) and use the fact that the gravlens output has a 
	one-to-one correspondence, allowing to separate the caustics at the same time.

	Warning (known bug): 
	 - show_plot=1 only shows only one plot in the screen

	Input:
	 - lens_model             <str> : Lens name (see gravlens manual table 3.1)
	 - mass_scale           <float> : Mass scale of the lens - "parameter 1"
	 - model_param_8        <float> : misc. lens parameter - often scale radio (depends on the lens model)
	 - model_param_9        <float> : misc. lens parameter - often scale radio (depends on the lens model)
	 - model_param_10       <float> : misc. lens parameter - often a power law index (depends on the lens model)
	 - galaxy_position <list float> : [x,y] position of the lens
	 - e_L                  <float> : lens ellipticity (default=0)
	 - theta_L              <float> : lens position angle (in degrees) with respect to the vertical 
					  (counterclockwise)
	 - shear                <float> : external shear amplitude
	 - theta_shear          <float> : external shear direction (in degrees)
	 - gravlens_params        <dic> : Contains the keys and values of the gravlens configuration 
					  (see default parameters at function set_gravlens_default, inside lens_parameters_new)
	 - caustic_CC_file        <str> : name of the output file with the caustic and CC positions
	 - gravlens_input_file    <str> : name of the input file used to run gravlens
	 - rad_curves_file        <str> : name of the file contaning the radial curves (optional - see write_to_file) 
	 - tan_curves_file        <str> : name of the file contaning the radial curves (optional - see write_to_file) 
	 - curves_plot            <str> : the name of the generated plot file (including the format, ex. 
					'curves_plot.png'). To see the formats available see matplotlib.pyplot help 
					(default='crit_curves.png'). If curves_plot=0 means no plots.
	 - show_plot              <int> : use 0 for 'no screen display' (default) and 1 for display on 
					  the screen.
	 - write_to_file          <int> : Option to write the curves to a file (see rad_curves_file and 
					  tan_curves_file ). The default is 0, meaning not to write to a file.
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
	 - rad_CC_x       <list> : x coordinates of points from the radial CC
	 - rad_CC_y       <list> : y coordinates of points from the radial CC
	 - tan_CC_x       <list> : x coordinates of points from the tangential CC
	 - tan_CC_y       <list> : y coordinates of points from the tangential CC
	 - rad_caustic_x  <list> : x coordinates of points from the radial caustic
	 - rad_caustic_y  <list> : y coordinates of points from the radial caustic
	 - tan_caustic_x  <list> : x coordinates of points from the tangential caustic
	 - tan_caustic_y  <list> : y coordinates of points from the tangential caustic
	 - <dict>                : all configuration variables used for running gravlens (including gridhi1)
	 - <str>                : name of the file ('caustic_CC_file') with the caustic and CC positions 
				  (its a gravlens output file)
	 - <str>                : name of the file with separated radial curves (CC + caustics) - rad_curves_file
	 - <str>                : name of the file with separated tangential curves (CC + caustics) - tan_curves_file
	 - <str>                : name of the file with the plots (curves_plot)

	"""

	logging.info('Starting run_find_CC for the lens model %s with mass_scale (model_param_1) = %s, model_param_8 = %s, model_param_9 = %s, model_param_10 = %s, galaxy_position = %s, e_L = %s, theta_L = %s, shear = %s, theta_shear = %s' % (lens_model, mass_scale, model_param_8, model_param_9, model_param_10, galaxy_position, e_L, theta_L, shear, theta_shear) )

	x_caustic, y_caustic, x_CC, y_CC, gravlens_params_updated = find_CC_new(lens_model, mass_scale, model_param_8, model_param_9, model_param_10, galaxy_position, e_L, theta_L, shear, theta_shear, gravlens_params, caustic_CC_file, gravlens_input_file)

	logging.debug('Ran find_CC_new')

	x1, y1, u1, v1, x2, y2, u2, v2 = np.loadtxt(caustic_CC_file, comments='#', unpack=True)

	logging.debug('Starting separation of the curves...')

	# =====================================================================
	# Separating the critical curves
	#nodes=[]
	#curves = separate_curves(x1, y1, x2, y2,nodes)
	CC_curves, start_idx, end_idx = separate_curves_a(x1, y1, x2, y2)


	if len(CC_curves) != 2:
		logging.warning('Function separate_curves found %d critical curve(s) (expected 2). It is probable that the curves separation function (separate_curves) did not separat them properly. Maybe you are approaching gravlens precision. Try changing units (ex., from arcsec to miliarcsec). Note also that some angles are not very well dealt by separate_curves (ex. 90 degrees).' % len(CC_curves) )
	else:
		radial_curve = CC_curves[0]
		tang_curve = CC_curves[1]


	if len(CC_curves) == 1 and theta_L != 0:
		logging.debug('Since only one CC was found, I will try to check if setting the lens angle to 0 solves the separation problem.')
		# try to get the CC with theta_L=0 and rotate back to theta_L
		# 1- run find_CC_new again, but with theta_L=0
		x_caustic, y_caustic, x_CC, y_CC, gravlens_params_if_clause = find_CC_new(lens_model, mass_scale, model_param_8, model_param_9, model_param_10, galaxy_position, e_L, 0, shear, theta_shear, gravlens_params, caustic_CC_file, gravlens_input_file) # ran with theta_L = 0
		# 2- read data (loadtxt)
		x1, y1, u1, v1, x2, y2, u2, v2 = np.loadtxt(caustic_CC_file, comments='#', unpack=True)
		# 3- separate the curves
		CC_curves, start_idx, end_idx = separate_curves_a(x1, y1, x2, y2)
		# 4- check again if there are 2 curves 
		if len(CC_curves) == 2:
			logging.debug('Two CC found. It seems the separation was fixed by rotating the CC back and forth.')
			radial_curve = CC_curves[0] # [ [x1_i], [y1_i], [x2_i], [y2_i] ]
			tang_curve = CC_curves[1] # [ [x1_i], [y1_i], [x2_i], [y2_i] ]			
			# 5- rotate back the curves (radial_curve[0], radial_curve[1], tang_curve[0], tang_curve[1])
			out_trans_rot, out_trans_rot2  = translation_and_rotation(radial_curve[0], radial_curve[1], galaxy_position[0], galaxy_position[1], theta_L, angle='degree'), translation_and_rotation(radial_curve[2], radial_curve[3], galaxy_position[0], galaxy_position[1], theta_L, angle='degree')
			radial_curve = out_trans_rot[0], out_trans_rot[1], out_trans_rot2[0], out_trans_rot2[1]
			out_trans_rot, out_trans_rot2  = translation_and_rotation(  tang_curve[0],   tang_curve[1], galaxy_position[0], galaxy_position[1], theta_L, angle='degree'), translation_and_rotation(  tang_curve[2],   tang_curve[3], galaxy_position[0], galaxy_position[1], theta_L, angle='degree')
			tang_curve   = out_trans_rot[0], out_trans_rot[1], out_trans_rot2[0], out_trans_rot2[1]
			# 5b- rotate the caustics
			rot_caustic1, rot_caustic2 = translation_and_rotation( u1, v1, galaxy_position[0], galaxy_position[1], theta_L, angle='degree'), translation_and_rotation( u2, v2, galaxy_position[0], galaxy_position[1], theta_L, angle='degree')
			u1, v1, u2, v2 = rot_caustic1[0], rot_caustic1[1], rot_caustic2[0], rot_caustic2[1]

		# 6- if, again, we haven't the 2 curves, go to 'else'
		if len(CC_curves) < 2:
			radial_curve = [[],[],[],[]]
			tang_curve = CC_curves[0]

		if len(CC_curves) > 2: # in NFW example of kappas = 0.4, rs = 74. , eL=0.5 and theta_L=90 this 'if' is needed (to get the x1,y1, etc positions of the curves)
			back_rotation_x1, back_rotation_u1, back_rotation_x2, back_rotation_u2 = translation_and_rotation( x1, y1, galaxy_position[0], galaxy_position[1], theta_L, angle='degree'), translation_and_rotation( u1, v1, galaxy_position[0], galaxy_position[1], theta_L, angle='degree'), translation_and_rotation( x2, y2, galaxy_position[0], galaxy_position[1], theta_L, angle='degree'), translation_and_rotation( u2, v2, galaxy_position[0], galaxy_position[1], theta_L, angle='degree')
			x1, y1, u1, v1, x2, y2, u2, v2 = back_rotation_x1[0], back_rotation_x1[1], back_rotation_u1[0], back_rotation_u1[1], back_rotation_x2[0], back_rotation_x2[1], back_rotation_u2[0], back_rotation_u2[1]



	# attempt to get the correct number of CC (iteration over the 'delta' parameter in separate_curves)
	delta = 0.1*np.sqrt( ((max(x1)-min(x1))**2) + ((max(y1)-min(y1))**2)) # the default from separate_curves
	delta_count = 0
	if len(CC_curves) > 2:
		logging.debug("Try to separate the curves correctly by iterating over 'delta'. Current value "
				 "of delta is %s and we found %s curves" % (delta, len(CC_curves)) )

		while len(CC_curves) > 2 and delta_count < max_delta_count:
			delta = delta*delta_increment
			CC_curves, start_idx, end_idx = separate_curves_a(x1, y1, x2, y2, delta=delta)
			delta_count += 1

		logging.debug("Ended loop to separate the curves iterating over 'delta'. Last value "
				 "of delta is %s (after %s iterations) and we found %s curves." % (delta, delta_count, len(CC_curves)) )

		if len(CC_curves) == 1:
			logging.warning('The iteration over delta did not work: it now found only one CC. Maybe lowering delta_increment will work out.')
			radial_curve = [[],[],[],[]]
			tang_curve = CC_curves[0]
		else:
			radial_curve = CC_curves[0] # [ [x1_i], [y1_i], [x2_i], [y2_i] ]
			tang_curve = CC_curves[1] # [ [x1_i], [y1_i], [x2_i], [y2_i] ]



	rad_CC_x,rad_CC_y, tan_CC_x, tan_CC_y = radial_curve[0], radial_curve[1], tang_curve[0], tang_curve[1]

	logging.debug('The radial and tangential CC have %d and %d points each.' % (len(rad_CC_x), len(tan_CC_x) ) )


	# Repeating the 1st element in the end of each array will avoid 'holes' in the plot
	# First, convert the array to a list
	rad_CC_x, rad_CC_y, tan_CC_x, tan_CC_y = list(rad_CC_x), list(rad_CC_y), list(tan_CC_x), list(tan_CC_y)

	if len(rad_CC_x) > 0:
		rad_CC_x.append(rad_CC_x[0])
		rad_CC_y.append(rad_CC_y[0])
	if len(tan_CC_x) > 0:
		tan_CC_x.append(tan_CC_x[0])
		tan_CC_y.append(tan_CC_y[0])

	rad_CC_x, rad_CC_y, tan_CC_x, tan_CC_y = np.array(rad_CC_x), np.array(rad_CC_y), np.array(tan_CC_x), np.array(tan_CC_y)

	# =====================================================================
	# separating the caustics curves
	caustic_curves = []
	for index in range( len(start_idx) ):
		caustic_curves.append([u1[start_idx[index]:end_idx[index]], v1[start_idx[index]:end_idx[index]], 
			       u2[start_idx[index]:end_idx[index]], v2[start_idx[index]:end_idx[index]]])

	logging.debug('Function separate_curves found %d caustics' % len(caustic_curves) )

	if len(caustic_curves) == 1:
		radial_curve = [[],[],[],[]]
		tang_curve = caustic_curves[0]
	else:
		radial_curve = caustic_curves[0] # [ [x1_i], [y1_i], [x2_i], [y2_i] ]
		tang_curve = caustic_curves[1] # [ [x1_i], [y1_i], [x2_i], [y2_i] ]

	rad_caustic_x, rad_caustic_y, tan_caustic_x, tan_caustic_y = radial_curve[0], radial_curve[1], tang_curve[0], tang_curve[1]

	logging.debug('The radial and tangential caustics have %d and %d points each.' % (len(rad_caustic_x), len(tan_caustic_x) ) )

	# Repeating the 1st element in the end of each array will avoid 'holes' in the plot
	# First, convert the array to a list
	rad_caustic_x, rad_caustic_y, tan_caustic_x, tan_caustic_y = list(rad_caustic_x), list(rad_caustic_y), list(tan_caustic_x), list(tan_caustic_y)

	if len(rad_caustic_x) > 0:
		rad_caustic_x.append(rad_caustic_x[0])
		rad_caustic_y.append(rad_caustic_y[0])
	if len(tan_caustic_x) > 0:
		tan_caustic_x.append(tan_caustic_x[0])
		tan_caustic_y.append(tan_caustic_y[0])

	rad_caustic_x, rad_caustic_y, tan_caustic_x, tan_caustic_y = np.array(rad_caustic_x), np.array(rad_caustic_y), np.array(tan_caustic_x), np.array(tan_caustic_y)

	logging.debug('Finished separation of the curves.')

	# plot
	logging.debug('Making the plots... \n')
	if curves_plot != 0: # curves_plot = 0 means you don't want any plots.
		plot_CC(tan_caustic_x, tan_caustic_y, rad_caustic_x, rad_caustic_y, tan_CC_x, tan_CC_y, rad_CC_x, rad_CC_y, curves_plot, show_plot)


	if write_to_file != 0: # wrihte_show =0 means you don't want save the files.
		logging.debug('Writting curves to files %s and %s.' % (tan_curves_file, rad_curves_file) )
		outtan = open(tan_curves_file,"w")			
		outtan.write("# Data file for tangential curves (critical and caustic) \n")
		outtan.write("# x1 x2 y1 y2 \n")
		outtan.write("# where (x1, x2) correspond to the critical curve and (y1,y2) correspond to caustic\n")
		outtan.write("#\n")
		for i in range(len(tan_CC_x)):
			outtan.write("%f %f %f %f\n" %(tan_CC_x[i],tan_CC_y[i],tan_caustic_x[i],tan_caustic_y[i]))
#
		outrad = open(rad_curves_file,"w")			
		outrad.write("# Data file for radial curves (critical and caustic) \n")
		outrad.write("# x1 x2 y1 y2 \n")
		outrad.write("# where (x1, x2) correspond to the critical curve and (y1,y2) correspond to caustic\n")
		outrad.write("#\n")
		for j in range(len(rad_CC_x)):
			outrad.write("%f %f %f %f\n" %(rad_CC_x[j],rad_CC_y[j],rad_caustic_x[j],rad_caustic_y[j]))


	return rad_CC_x, rad_CC_y, tan_CC_x, tan_CC_y, rad_caustic_x, rad_caustic_y, tan_caustic_x, tan_caustic_y





