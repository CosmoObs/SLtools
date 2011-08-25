#!/usr/bin/env python
# ==================================
# Authors:
# Pedro Ferreira - pferreira@dfte.ufrn.br
# ==================================

""" Package to deal with gravlens configuration """

##@package lens_parameters
# 
#
# Deals with gravlens configuration. The function set_gravlens_default defines default values for 
# gravlens configuration. Function lens_parameters generates the basis gravlens configuration
# so the user can start his calculations with gravlens. 


#=======================================================================================================
def set_gravlens_default():
    """
    Defines the configuration of gravlens

    In order to make broader application of these modules we have set some default values more 
    accurate. For instance, separating the critical curves and caustics (see pipelines.find_CC.run) of 
    some models demand crittol = 1e-7 and inttol = 1e-8 (see description below). The identified cases 
    where this is needed are the NFW with the combination of parameters: kappas = [0.4], 
    rs = [34., 54., 94.], theta = [0, 180] and e_L=0.5 (note that the problem for these values were not 
    in gravlens itself, but in the separation of the curves).

    The default configurations is:

     - gridmode      <int> = 1 : the grid type gravlens will use (1: standard polar grid, 2: use 
    		                 information from simple convex critical lines)
     - gridhi1   <float> = 150 : the grid size in the image plane
     - xtol    <float> = 1e-10 : Tolerance on image positions in numerical root finding
     - crittol  <float> = 1e-7 : Tolerance for finding critical curves
     - inttol   <float> = 1e-8 : Tolerance for numerical integrals
     - maxlev        <int> = 4 : Deepest level of subgrid near critical curves
     - gallev        <int> = 3 : Deepest level of subgrid near galaxies (other than the primary)
     - imglev        <int> = 3 : Deepest level of subgrid near images
     - ngrid1       <int> = 32 : dimension of top grid (in radius)
     - ngrid2       <int> = 32 : dimension of top grid (in angle)

    Output:
     - gravlens_params_default <dict> : contains the sltools default keys for gravlens configuration.
    				        The keys and its values are 'gridmode': 1, 'gridhi1' : 150, 
				        'xtol' : 1e-10, 'crittol' : 1e-7, 'inttol' : 1e-8, 'maxlev' : 4, 
    				        'gallev' : 3, 'imglev' : 3, 'ngrid1' : 32, 'ngrid2' : 32 

    """
    gravlens_params_default = {'gridmode': 1, 'gridhi1' : 150, 'xtol' : 1e-10, 'crittol' : 1e-7, 'inttol' : 1e-8, 'maxlev' : 4, 'gallev' : 3, 'imglev' : 3, 'ngrid1' : 32, 'ngrid2' : 32 }


    return gravlens_params_default


def lens_parameters(lens_model, mass_scale, model_param_8, model_param_9, model_param_10, galaxy_position=(0,0), e_L=0, theta_L=0, shear=0, theta_shear=0, gravlens_params={}):
    """
    Determine the 'base string' that defines the lens and gravlens properties.

    The advantage of this function is that the user do not need to know the gravlens syntax (only 
    the lenses available and its parameters).

    Input:
     - lens_model        <str> : Lens name (see gravlens manual table 3.1)
     - mass_scale      <float> : Mass scale of the lens - "parameter 1"
     - model_param_8   <float> : misc. lens parameter - often scale radio (depends on the lens model)
     - model_param_9   <float> : misc. lens parameter - often scale radio (depends on the lens model)
     - model_param_10  <float> : misc. lens parameter - often a power law index (depends on the lens model)
     - galaxy_position  <list> : [x,y] position of the lens
     - e_L             <float> : Horizontal central position for output (cut) image
     - theta_L         <float> : Vertical central position for output (cut) image
     - shear           <float> : 'pixel' or 'degrees' for size (x_size,y_size) values
     - theta_shear     <float> : Horizontal size (in pixels) of output image
     - gravlens_params    dict : contains the keys and values of the gravlens configuration (see default 
                                 parameters at function set_gravlens_default,inside init_gravlens_parameter)

    Output:
     - <str> : a string with all the lines gravlens needs to its configuration
     - <str> : a string (one line) that only defines the lens used (in gravlens format) 
     - <dic> : all gravlens defaults (defined at 'set_gravlens_default') plus the ones used as input

    """

    galaxy_position = list(galaxy_position)

    gravlens_params_updated = set_gravlens_default() # sets the gravlens default

    gravlens_params_updated.update( gravlens_params ) # update with the user options
    
    # define the lens model
    setlens = '%s %f %f %f %f %f %f %f %f %f %f'% (lens_model  , mass_scale, galaxy_position[0], galaxy_position[1], e_L, theta_L, shear, theta_shear, model_param_8, model_param_9, model_param_10 )

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
