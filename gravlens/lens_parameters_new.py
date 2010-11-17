##@package lens_parameters_new
# deals with gravlens configuration
#
# 
#

def set_gravlens_default():
	"""
	Defines the configuration of gravlens

	The default configurations is:

	gridhi1 : 150 # the grid size in the image plane
	xtol    : 1e-10 # Tolerance on image positions in numerical root finding
	crittol : 1e-6 # Tolerance for finding critical curves
	inttol  : 1e-6 # Tolerance for numerical integrals
	maxlev  : 4 # Deepest level of subgrid near critical curves
	gallev  : 3 # Deepest level of subgrid near galaxies (other than the primary)
	imglev  : 3 # Deepest level of subgrid near images
	ngrid1  : 30 # dimension of top grid (in radius)
	ngrid2  : 30 # dimension of top grid (in angle)
	plotmode: 1 # Specifies file format for curve plots

	Any of these parameters can be modified through the input.
	
	Output:
	- <dic> : contains the updated keys for gravlens configuration
	"""
	gravlens_params_default = { 'gridhi1' : 150, 'xtol' : 1e-10, 'crittol' : 1e-6, 'inttol' : 1e-6, 'maxlev' : 4, 'gallev' : 3, 'imglev' : 3, 'ngrid1' : 30, 'ngrid2' : 30, 'plotmode' : 1 }


	return gravlens_params_default


def lens_parameters_new(lens, mass_scale, model_param_8, model_param_9, model_param_10, galaxy_position=[0,0], e_L=0, theta_L=0, shear=0, theta_shear=0, gravlens_params={}):
	"""
	Determine the 'basis string' that defines the lens and gravlens properties.

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


	Output:
	 - (ndarray, header) : resultant image array and (updated) header instance

	"""

	gravlens_params_updated = set_gravlens_default() # sets the gravlens default

	gravlens_params_updated.update( gravlens_params ) # update with the user options
	
	# define the lens model
	setlens = '%s %f %f %f %f %f %f %f %f %f %f'% (lens, mass_scale, galaxy_position[0], galaxy_position[1], e_L, theta_L, shear, theta_shear, model_param_8, model_param_9, model_param_10 )

	inputlens = """set gridhi1=%0.12f
set ngrid1 = %d
set ngrid2 = %d
set xtol = %0.12f 
set crittol = %0.12f 
set inttol = %0.12f 
set maxlev = %d
set gallev = %d 
set imglev = %d 
set plotmode = %d
startup 1 1
%s
0 0 0 0 0 0 0 0 0 0\n""" % ( float(gravlens_params_updated['gridhi1']),
                             int(gravlens_params_updated['ngrid1']),
                             int(gravlens_params_updated['ngrid2']),
                             float(gravlens_params_updated['xtol']),
                             float(gravlens_params_updated['crittol']),
                             float(gravlens_params_updated['inttol']),
                             int(gravlens_params_updated['maxlev']),
                             int(gravlens_params_updated['gallev']),
                             int(gravlens_params_updated['imglev']),
                             int(gravlens_params_updated['plotmode']), setlens) 

	return inputlens, setlens, gravlens_params_updated
