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
# The pipeline run concatenates the functions gravlens.find_cc.find_cc, 
# geometry.separate_curves.separate_curves and gravlens.find_cc.plot_cc, optimizing the curves separation 
# (see run for details). For high-level documentation, see 
# http://twiki.linea.gov.br/bin/view/StrongLensing/SLtoolsGravlens#Module_find_cc_py


import logging
import numpy as np
import time

from sltools.coordinate.translation_and_rotation import translate_and_rotate_coord_system
from sltools.geometry.separate_curves import separate_curves
from sltools.gravlens.find_cc import critcurves, find_cc, plot_cc



#=======================================================================================================
def run(lens_model, mass_scale, model_param_8, model_param_9, model_param_10, galaxy_position=[0.,0.], e_L=0, theta_L=0, shear=0, theta_shear=0, gravlens_params={}, caustic_CC_file='crit.txt', gravlens_input_file='gravlens_CC_input.txt', rad_curves_file='lens_curves_rad.dat', tan_curves_file='lens_curves_tan.dat', curves_plot='crit-caust_curves.png', show_plot=0, write_to_file=0, max_delta_count=20, delta_increment=1.1, grid_factor=5., grid_factor2=3., max_iter_number=20, min_n_lines=200, gridhi1_CC_factor=2., accept_res_limit=2E-4):
	""" 
	This is a pipeline that runs 'find_cc', 'separate_curves' and 'plot_cc'. For details of these 
	functions, see their documentation.

	Besides running the functions above mentioned, run also treats for some cases which the CC
	are not separated correctly. See details below.

	As we expect 2 critical curves, if the function separate_curves finds a different number, we 
	question it. If we find only one CC, then the separation code did not worked properly, and we 
	noticed that for some values of theta_L this occurred quite often. So we find the CC again for
	theta_L=0, separate the CC for this value of theta_L and rotate the separated curves back to the
	original value of theta_L. This proved to be a good solution.

	In the case that more than 2 CC are found, we try an iterating process with the control parameter
	of separate_curves 'delta' (see separate_curves documentation). So we increase the value of delta 
	until there are 2 CC.

	The combination of these 2 implementations above with the default gravlens parameters used in 
	'lens_parameters' allows a good success on separating the curves.
	For more details, see http://twiki.linea.gov.br/bin/view/StrongLensing/SLtoolsGravlens#Module_find_CC_new_py.

	We separate only the image plane curves (CC) and use the fact that the gravlens output has a 
	one-to-one correspondence, allowing to separate the caustics at the same time.

	Warning (known bug): 
	 - show_plot=1 shows only one plot in the screen

	Input:
	 - lens_model               str : Lens name (see gravlens manual table 3.1)
	 - mass_scale             float : Mass scale of the lens - "parameter 1"
	 - model_param_8          float : misc. lens parameter - often scale radio (depends on the lens model)
	 - model_param_9          float : misc. lens parameter - often scale radio (depends on the lens model)
	 - model_param_10         float : misc. lens parameter - often a power law index (depends on the lens model)
	 - galaxy_position  [float,...] : [x,y] position of the lens
	 - e_L                    float : lens ellipticity (ellipticity = 1 - q, where 0<q<=1 is the axial ratio) 
	 - theta_L                float : lens position angle (in degrees) with respect to the vertical 
					                  (counterclockwise)
	 - shear                  float : external shear amplitude
	 - theta_shear            float : external shear direction (in degrees)
	 - gravlens_params          dic : Contains the keys and values of the gravlens configuration 
					  (see default parameters at function set_gravlens_default, inside init_gravlens_parameters)
	 - caustic_CC_file          str : name of the output file with the caustic and CC positions
	 - gravlens_input_file      str : name of the input file used to run gravlens
	 - rad_curves_file          str : name of the file containing the radial curves (optional - see write_to_file) 
	 - tan_curves_file          str : name of the file containing the radial curves (optional - see write_to_file) 
	 - curves_plot              str : the name of the generated plot file (including the extension, ex. 
					'curves_plot.png'). To see the extensions available see matplotlib.pyplot help 
					(default='crit_curves.png'). If curves_plot=0 means no plots.
	 - show_plot                int : use 0 for 'no screen display' (default) and 1 for display on 
					  the screen.
	 - write_to_file            int : Option to write the curves to a file (see rad_curves_file and 
					  tan_curves_file ). The default is 0, meaning not to write to a file.
     - delta_increment        float : 
	 - grid_factor            float : used in the iteration to find the critical curves. Each time no 
				       curves were found, the grid size is divided by grid_factor
	 - min_n_lines              int : minimum number of lines demanded in the CC file to decrease the grid by grid_factor_2 
	 - grid_factor2           float : if the CC file is composed of less than min_n_lines, the code divides
				       the grid size by grid_factor2 
	 - max_iter_number          int : maximum number of time the grid size will be divided by grid_factor
	 - gridhi1_CC_factor      float : the final size of the grid will be gridhi1_CC_factor * (the furthest CC point)
	 - accept_res_limit       float : if the furthest CC point is closer than accept_res_limit from the lens 
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

	logging.info('Starting run for the lens model %s with mass_scale (model_param_1) = %s, model_param_8 = %s, model_param_9 = %s, model_param_10 = %s, galaxy_position = %s, e_L = %s, theta_L = %s, shear = %s, theta_shear = %s' % (lens_model, mass_scale, model_param_8, model_param_9, model_param_10, galaxy_position, e_L, theta_L, shear, theta_shear) )

	x_caustic, y_caustic, x_CC, y_CC, gravlens_params_updated = find_cc(lens_model, mass_scale, model_param_8, model_param_9, model_param_10, galaxy_position, e_L, theta_L, shear, theta_shear, gravlens_params, caustic_CC_file, gravlens_input_file)

	logging.debug('Ran find_cc')

	x1, y1, u1, v1, x2, y2, u2, v2 = np.loadtxt(caustic_CC_file, comments='#', unpack=True)

	logging.debug('Starting separation of the curves...')

	# =====================================================================
	# Separating the critical curves
	#nodes=[]
	#curves = separate_curves(x1, y1, x2, y2,nodes)
	CC_curves, start_idx, end_idx = separate_curves(x1, y1, x2, y2)


	if len(CC_curves) != 2:
		logging.warning('Function separate_curves found %d critical curve(s) (expected 2). It is probable that the curves separation function (separate_curves) did not separate them properly. Maybe you are approaching gravlens precision. Try changing units (ex., from arcsec to miliarcsec). Note also that some angles are not very well dealt by separate_curves (ex. 90 degrees).' % len(CC_curves) )
	else:
		radial_curve = CC_curves[0]
		tang_curve = CC_curves[1]


	if len(CC_curves) == 1 and theta_L != 0:
		logging.debug('Since only one CC was found, I will try to check if setting the lens angle to 0 solves the separation problem.')
		# try to get the CC with theta_L=0 and rotate back to theta_L
		# 1- run find_cc again, but with theta_L=0
		x_caustic, y_caustic, x_CC, y_CC, gravlens_params_if_clause = find_cc(lens_model, mass_scale, model_param_8, model_param_9, model_param_10, galaxy_position, e_L, 0, shear, theta_shear, gravlens_params_updated, caustic_CC_file, gravlens_input_file) # ran with theta_L = 0
		# 2- read data (loadtxt)
		x1, y1, u1, v1, x2, y2, u2, v2 = np.loadtxt(caustic_CC_file, comments='#', unpack=True)
		# 3- separate the curves
		CC_curves, start_idx, end_idx = separate_curves(x1, y1, x2, y2)
		# 4- check again if there are 2 curves 
		if len(CC_curves) == 2:
			logging.debug('Two CC found. It seems the separation was fixed by rotating the CC back and forth.')
			radial_curve = CC_curves[0] # [ [x1_i], [y1_i], [x2_i], [y2_i] ]
			tang_curve = CC_curves[1] # [ [x1_i], [y1_i], [x2_i], [y2_i] ]			
			# 5- rotate back the curves (radial_curve[0], radial_curve[1], tang_curve[0], tang_curve[1])
			out_trans_rot, out_trans_rot2  = translate_and_rotate_coord_system(radial_curve[0], radial_curve[1], galaxy_position[0], galaxy_position[1], theta_L, angle='degree'), translate_and_rotate_coord_system(radial_curve[2], radial_curve[3], galaxy_position[0], galaxy_position[1], theta_L, angle='degree')
			radial_curve = out_trans_rot[0], out_trans_rot[1], out_trans_rot2[0], out_trans_rot2[1]
			out_trans_rot, out_trans_rot2  = translate_and_rotate_coord_system(  tang_curve[0],   tang_curve[1], galaxy_position[0], galaxy_position[1], theta_L, angle='degree'), translate_and_rotate_coord_system(  tang_curve[2],   tang_curve[3], galaxy_position[0], galaxy_position[1], theta_L, angle='degree')
			tang_curve   = out_trans_rot[0], out_trans_rot[1], out_trans_rot2[0], out_trans_rot2[1]
			# 5b- rotate the caustics
			rot_caustic1, rot_caustic2 = translate_and_rotate_coord_system( u1, v1, galaxy_position[0], galaxy_position[1], theta_L, angle='degree'), translate_and_rotate_coord_system( u2, v2, galaxy_position[0], galaxy_position[1], theta_L, angle='degree')
			u1, v1, u2, v2 = rot_caustic1[0], rot_caustic1[1], rot_caustic2[0], rot_caustic2[1]

		# 6- if, again, we haven't the 2 curves, go to 'else'
		if len(CC_curves) < 2:
			radial_curve = [[],[],[],[]]
			tang_curve = CC_curves[0]

		if len(CC_curves) > 2: # in NFW example of kappas = 0.4, rs = 74. , eL=0.5 and theta_L=90 this 'if' is needed (to get the x1,y1, etc positions of the curves)
			back_rotation_x1, back_rotation_u1, back_rotation_x2, back_rotation_u2 = translate_and_rotate_coord_system( x1, y1, galaxy_position[0], galaxy_position[1], theta_L, angle='degree'), translate_and_rotate_coord_system( u1, v1, galaxy_position[0], galaxy_position[1], theta_L, angle='degree'), translate_and_rotate_coord_system( x2, y2, galaxy_position[0], galaxy_position[1], theta_L, angle='degree'), translate_and_rotate_coord_system( u2, v2, galaxy_position[0], galaxy_position[1], theta_L, angle='degree')
			x1, y1, u1, v1, x2, y2, u2, v2 = back_rotation_x1[0], back_rotation_x1[1], back_rotation_u1[0], back_rotation_u1[1], back_rotation_x2[0], back_rotation_x2[1], back_rotation_u2[0], back_rotation_u2[1]



	# attempt to get the correct number of CC (iteration over the 'delta' parameter in separate_curves)
	delta = 0.1*np.sqrt( ((max(x1)-min(x1))**2) + ((max(y1)-min(y1))**2)) # the default from separate_curves
	delta_count = 0
	if len(CC_curves) > 2:
		logging.debug("Try to separate the curves correctly by iterating over 'delta'. Current value "
				 "of delta is %s and we found %s curves" % (delta, len(CC_curves)) )

		while len(CC_curves) > 2 and delta_count < max_delta_count:
			delta = delta*delta_increment
			CC_curves, start_idx, end_idx = separate_curves(x1, y1, x2, y2, delta=delta)
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
	if curves_plot != 0: # curves_plot = 0 means you don't want any plots.
		logging.debug('Making the plots... \n')
		plot_cc(tan_caustic_x, tan_caustic_y, rad_caustic_x, rad_caustic_y, tan_CC_x, tan_CC_y, rad_CC_x, rad_CC_y, curves_plot, show_plot)


	if write_to_file != 0: # write_show =0 means you don't want save the files.
		logging.debug('Writing curves to files %s and %s.' % (tan_curves_file, rad_curves_file) )
		outtan = open(tan_curves_file,"w")			
		outtan.write("# Data file for tangential curves (critical and caustic) \n")
		outtan.write("# File generated by run on %s for the lens model %s with mass_scale (model_param_1) = %s, model_param_8 = %s, model_param_9 = %s, model_param_10 = %s, galaxy_position = %s, e_L = %s, theta_L = %s, shear = %s, theta_shear = %s \n" % (time.ctime(), lens_model, mass_scale, model_param_8, model_param_9, model_param_10, galaxy_position, e_L, theta_L, shear, theta_shear) )
		outtan.write("# x1 x2 y1 y2 \n")
		outtan.write("# where (x1, x2) correspond to the critical curve and (y1,y2) correspond to caustic\n")
		outtan.write("#\n")
		for i in range(len(tan_CC_x)):
			outtan.write("%f %f %f %f\n" %(tan_CC_x[i],tan_CC_y[i],tan_caustic_x[i],tan_caustic_y[i]))
#
		outrad = open(rad_curves_file,"w")			
		outrad.write("# Data file for radial curves (critical and caustic) \n")
		outrad.write("# File generated by run on %s for the lens model %s with mass_scale (model_param_1) = %s, model_param_8 = %s, model_param_9 = %s, model_param_10 = %s, galaxy_position = %s, e_L = %s, theta_L = %s, shear = %s, theta_shear = %s \n" % (time.ctime(), lens_model, mass_scale, model_param_8, model_param_9, model_param_10, galaxy_position, e_L, theta_L, shear, theta_shear) )
		outrad.write("# x1 x2 y1 y2 \n")
		outrad.write("# where (x1, x2) correspond to the critical curve and (y1,y2) correspond to caustic\n")
		outrad.write("#\n")
		for j in range(len(rad_CC_x)):
			outrad.write("%f %f %f %f\n" %(rad_CC_x[j],rad_CC_y[j],rad_caustic_x[j],rad_caustic_y[j]))


	return rad_CC_x, rad_CC_y, tan_CC_x, tan_CC_y, rad_caustic_x, rad_caustic_y, tan_caustic_x, tan_caustic_y





