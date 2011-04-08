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
# This package has functions to i) get the caustic and the critical curves of a given model (see find_CC_new) and ii) 
# plot the caustics and the critical curves (plot_CC).
# The curves are determined using an iterating method that uses gravlens (see find_CC_new).
# The curves are identifyed as radial and tangencial using the function sltools.geometry.separate_curves.

import os
import logging
import numpy as np
import matplotlib.pyplot as pyplot

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
def find_CC_new(lens_model, mass_scale, model_param_8, model_param_9, model_param_10, galaxy_position=[0,0], e_L=0, theta_L=0, shear=0, theta_shear=0, gravlens_params={}, caustic_CC_file='crit.txt', gravlens_input_file='gravlens_CC_input.txt'):
	""" 
	Determines the caustics and critical curves (CC) for a single component lens model. Uses an iterative
	procedure to determine the grid size and the minimum number of points on the curves.

	The CC size may vary considerably depending on the lens-source configuration. Because of this, 
	we developed an iterative method to find the CC, that involves the size of the grid used by 
	gravlens (gridhi1). The initial value of gridhi1 must be an upper limit for the size of the CC, 
	and find_CC_new will try to find the CC. Each time gravlens can't find the CC, we lower	gridhi1 
	by a factor of grid_factor (=5), in order to increase precision. This continues on until 
	gravlens finds the CC or the number of iterations, max_iterations_number (=20), is reached. After 
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
	 - lens_model          <str> : lens name (see gravlens manual table 3.1)
	 - mass_scale        <float> : mass scale of the lens - "parameter 1"
	 - model_param_8     <float> : misc. lens parameter - often scale radio
	 - model_param_9     <float> : misc. lens parameter - often scale radio
	 - model_param_10    <float> : misc. lens parameter - often a power law index
	 - galaxy_position    <list> : [x,y] position of the lens
	 - e_L               <float> : lens ellipticity (default=0)
	 - theta_L           <float> : lens position angle (in degrees) with respect to the vertical 
					  (counterclockwise)
	 - shear             <float> : external shear amplitude
	 - theta_shear       <float> : external shear direction (in degrees)
	 - gravlens_params     <dic> : contains the keys and values of the gravlens configuration 
					  (see default parameters at function set_gravlens_default,inside lens_parameters_new)
	 - caustic_CC_file     <str> : name of the output file with the caustic and CC positions
	 - gravlens_input_file <str> : name of the input file used to run gravlens

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
	logging.debug('Determined the strings that defines the mass model in gravlens through the \'lens_parameters\' function')
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


	x1, y1, u1, v1 = np.loadtxt(caustic_CC_file, usecols = (0,1,2,3), unpack=True) # CC_x, CC_y, caustic_x, caustic_y

	#-----------------------------------------------------------------------------------------------------------
	# redefine gridhi1 according to the CC size

	index_CC = np.argmax(x1**2 + y1**2)
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
	if max( max(np.abs(u1)), max(np.abs(v1)) ) < acceptable_res_limit:
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
	 - <file> : the plot file named plot_filename will be generated only if show_plot=0

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
def run_find_CC(lens_model, mass_scale, model_param_8, model_param_9, model_param_10, galaxy_position=[0.,0.], e_L=0, theta_L=0, shear=0, theta_shear=0, gravlens_params={}, caustic_CC_file='crit.txt', gravlens_input_file='gravlens_CC_input.txt', rad_curves_file='lens_curves_rad.dat', tan_curves_file='lens_curves_tan.dat', curves_plot='crit-caust_curves.png', show_plot=0, write_to_file=0):
	""" 
	This is a pipeline that runs 'find_CC_new', 'separate_CC' and
	'plot_CC'. For details of these functions, see their documentation.

	We separate only the image plane curves (CC) and use the fact that the gravlens output has a 
	one-to-one correspondence, allowing to separate the caustics at the same time.

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
					  (see default parameters at function set_gravlens_default, inside lens_parameters_new)
	 - caustic_CC_file        <str> : name of the output file with the caustic and CC positions
	 - gravlens_input_file    <str> : name of the input file used to run gravlens
	 - rad_curves_file        <str> : 
	 - tan_curves_file        <str> :  (default='crit_tan.txt')
	 - curves_plot            <str> : the name of the generated plot file (including the format, ex. 
					'curves_plot.png'). To see the formats available see matplotlib.pyplot help 
					(default='crit_curves.png'). If curves_plot=0 means no plots.
	 - show_plot              <int> : use 0 for 'no screen display' (default) and 1 for display on 
					  the screen.
	 - write_to_file          <int> : Option to write the curves to a file (see rad_curves_file and 
					  tan_curves_file ). The default is 0, meaning not to write to a file.

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
	 - <file>                : file named 'caustic_CC_file' with the caustic and CC positions
	 - <file>                : separated radial curves (CC + caustics) - rad_curves_file
	 - <file>                : separated tangential curves (CC + caustics) - tan_curves_file
	 - <file>                : curves_plot

	Warning: 
			 show_plot=1 only shows one plot in the screen

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
		logging.warning('Function separate_curves found %d critical curve(s) (expected 2). It is probable that the curves separation function (separate_curves) did not separated them properly. Maybe you are approaching gravlens precision. Try changing units (ex., from arcsec to miliarcsec). Note also that some angles are not very well dealt by separate_curves (ex. 90 degrees).' % len(CC_curves) )

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
		else:
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

	logging.debug('Making the plots... \n')
	# plot
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





