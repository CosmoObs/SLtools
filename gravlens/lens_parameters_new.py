#!/usr/bin/env python
# ==================================
# Authors:
# Pedro Ferreira - pferreira@dfte.ufrn.br
# ==================================

""" Package to deal with gravlens configuration """

##@package lens_parameters_new
# 
#
# Deals with gravlens configuration. The function set_gravlens_default defines default values for 
# gravlens configuration. Function lens_parameters_new generates the basis gravlens configuration
# in order to the user start his calculations with gravlens. 


#=======================================================================================================
def set_gravlens_default():
	"""
	Defines the configuration of gravlens

	The default configurations is:

	gridmode : 1 # the grid type gravlens will use (1: standard polar grid, 2: use information from 
		       simple convex critical lines)
	gridhi1  : 150 # the grid size in the image plane
	xtol     : 1e-10 # Tolerance on image positions in numerical root finding
	crittol  : 1e-6 # Tolerance for finding critical curves
	inttol   : 1e-6 # Tolerance for numerical integrals
	maxlev   : 4 # Deepest level of subgrid near critical curves
	gallev   : 3 # Deepest level of subgrid near galaxies (other than the primary)
	imglev   : 3 # Deepest level of subgrid near images
	ngrid1   : 30 # dimension of top grid (in radius)
	ngrid2   : 30 # dimension of top grid (in angle)

	Output:
	- gravlens_params_default <dict> : contains the updated keys for gravlens configuration

	"""
	gravlens_params_default = {'gridmode': 1, 'gridhi1' : 150, 'xtol' : 1e-10, 'crittol' : 1e-6, 'inttol' : 1e-6, 'maxlev' : 4, 'gallev' : 3, 'imglev' : 3, 'ngrid1' : 30, 'ngrid2' : 30 }


	return gravlens_params_default


def lens_parameters_new(lens, mass_scale, model_param_8, model_param_9, model_param_10, galaxy_position=[0,0], e_L=0, theta_L=0, shear=0, theta_shear=0, gravlens_params={}):
	"""
	Determine the 'basis string' that defines the lens and gravlens properties.

	The advantage of this function is that the user do not need to know the gravlens syntax (only 
	the lenses available and its parameters).

	Input:
	 - lens              <str> : Lens name (see gravlens manual table 3.1)
	 - mass_scale      <float> : Mass scale of the lens - "parameter 1"
	 - model_param_8   <float> : Vertical size (in pixels) of output image
	 - model_param_9   <float> : Vertical size (in pixels) of output image
	 - model_param_10  <float> : Vertical size (in pixels) of output image
	 - galaxy_position  <list> : [x,y] position of the lens
	 - e_L             <float> : Horizontal central position for output (cut) image
	 - theta_L         <float> : Vertical central position for output (cut) image
	 - shear           <float> : 'pixel' or 'degrees' for size (x_size,y_size) values
	 - theta_shear     <float> : Horizontal size (in pixels) of output image


	Output:
	 - <str> : a string with all the lines gravlens needs to its configuration
	 - <str> : a string line that only defines the lens used (in gravlens format) 
	 - <dic> : all gravlens defaults (defined at 'set_gravlens_default') plus the
		   ones used as input

	"""

	gravlens_params_updated = set_gravlens_default() # sets the gravlens default

	gravlens_params_updated.update( gravlens_params ) # update with the user options
	
	# define the lens model
	setlens = '%s %f %f %f %f %f %f %f %f %f %f'% (lens, mass_scale, galaxy_position[0], galaxy_position[1], e_L, theta_L, shear, theta_shear, model_param_8, model_param_9, model_param_10 )

	inputlens = """gridmode %d
set gridhi1=%0.12f
set ngrid1 = %d
set ngrid2 = %d
set xtol = %0.12f 
set crittol = %0.12f 
set inttol = %0.12f 
set maxlev = %d
set gallev = %d 
set imglev = %d 
startup 1 1
%s
0 0 0 0 0 0 0 0 0 0\n""" % ( int(gravlens_params_updated['gridmode']),
			     float(gravlens_params_updated['gridhi1']),
                             int(gravlens_params_updated['ngrid1']),
                             int(gravlens_params_updated['ngrid2']),
                             float(gravlens_params_updated['xtol']),
                             float(gravlens_params_updated['crittol']),
                             float(gravlens_params_updated['inttol']),
                             int(gravlens_params_updated['maxlev']),
                             int(gravlens_params_updated['gallev']),
                             int(gravlens_params_updated['imglev']), setlens) 

	return inputlens, setlens, gravlens_params_updated
