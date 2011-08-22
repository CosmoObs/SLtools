#!/usr/bin/env python
# ==================================
# Authors:
# Pedro Ferreira - pferreira@dfte.ufrn.br
# ==================================

""" Package to lens finite sources with gravlens and measure the images. """

##@package lens_finite_sources_new
# 
#
# This package contains functions to lens finite sources (lensing), extract the image plane images 
# determining mergers, plot sources and images with the caustics and critical curves and apply 
# measurement methods to individual images. Also, the pipeline lens_finite_sources_new concatenates and 
# runs the above functions.

from __future__ import division
import os
import logging
import numpy as np
import pyfits
import random

from sltools.gravlens.lens_parameters_new import lens_parameters_new
from sltools.gravlens.find_CC_new import run_find_CC
from sltools.gravlens.find_CC_new import plot_CC
from sltools.image.imcp import elliminate_disconected as elliminate_disconected
from sltools.image.sextractor import run_segobj



def source_model_sampler(source_type, rs_i, rs_f, es_i, es_f, thetas_i, thetas_f, ns_i=0.5, ns_f=0.5, nover=1, mag=24.0):
    """ 
    Generates dictionary with source properties ("source_model") custom made for lens_finite_sources.

    Given variable values var_i and var_f ('var' being radius, ellipticity, theta or Sersic index 'n'), 
    picks randomly a value between var_i and var_f. If var_i = var_f, var_i is taken.

    Input:
     - source_type <str> : Can be 'sersic' for a Sersic profile, or 'uniform' for a constant surface brightness
     - rs_i      <float> : Inferior limit for the source radius
     - rs_f      <float> : Superior limit for the source radius
     - es_i      <float> : Inferior limit for the source ellipticity
     - es_f      <float> : Superior limit for the source ellipticity
     - thetas_i  <float> : Inferior limit for the source position angle
     - thetas_f  <float> : Superior limit for the source position angle
     - ns_i      <float> : Inferior limit for the source sersic index 'n'
     - ns_f      <float> : Superior limit for the source sersic index 'n'
     - nover     <int>   : Oversampling of the image plane, with nover steps in each direction. 
     - mag       <float> : Magnitude of the source
     
    Output:
     - <dict> : Dictionary with the source properties. The keys are 'source_type', 'rs', 'es', 'thetas', 
                'nover', 'ns', 'mag'.

    """

    # radius --------------------------------------------------
    if rs_i == rs_f:
        rs = rs_i
    else:
        rs = rs_i + (rs_f - rs_i)*random.random()
    # ellipticity ---------------------------------------------
    if es_i == es_f:
        es = es_i
    else:
        es = es_i + (es_f - es_i)*random.random()
    # position angle ------------------------------------------
    if thetas_i == thetas_f:
        thetas = thetas_i
    else:
        thetas = thetas_i + (thetas_f - thetas_i)*random.random()
    # sersic n index
    if ns_i == ns_f:
        ns = ns_i
    else:
        ns = ns_i + (ns_f - ns_i)*random.random()

    return { 'source_type': source_type, 'rs': rs , 'es': es, 'thetas': thetas, 'nover': nover, 'ns': ns, 'mag': mag}


def lensing(lens_model, mass_scale, model_param_8, model_param_9, model_param_10, galaxy_position, e_L, 
            theta_L, shear, theta_shear, gravlens_params, dimpix, source_centers, source_model, 
            ref_magzpt, reference_band, base_name='sbmap', nover_max=3):
    """
    Lenses sources for the given lens model. 
    
    The lensing is performed with the gravlens software, using command SBmap2.

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
     - gravlens_params   <dic> : Contains the keys and values of the gravlens configuration , dimpix, source_centers, source_model, 
     - ref_magzpt      <float> : The zero point magnitude to be used. It only affects the source fluxes.
     - reference_band    <str> :
     - base_name         <str> : 
     - nover_max         <int> :

    Output:
     - image_names <list> :

    """

    inputlens, setlens, same_gravlens_params = lens_parameters_new(lens_model, mass_scale, model_param_8, 
                                                  model_param_9, model_param_10, galaxy_position, e_L, 
                                                  theta_L, shear, theta_shear, gravlens_params)

    half_frame_size = gravlens_params['gridhi1']
    Npix = int( (2*half_frame_size) / dimpix )
    image_names = []
    f199 = open('sersicinput.txt', 'w')
    f199.write(inputlens)
    for i in range (0,len(source_centers)):
        source_type = source_model[i]['source_type'] # type of the source (can be 'sersic' or 'uniform')
        f199.write('setsource 1\n')
        # Nover increase resolution if the source is smaller than the pixel
        nover = source_model[i]['nover']
        while source_model[i]['rs'] < dimpix/nover and nover < nover_max:
            nover += 1
        source_model[i]['nover'] = nover
        # -------------------------------------------------------------------
        if source_type == 'sersic':
            # determine 'totalsourceplaneflux', according to the source magnitude
            totalsourceplaneflux = 10**((2./5)*(ref_magzpt - source_model[i]['mag']))
            f199.write('%s %f %.6f %.6f %f %f %.6f 0 %f macro\n' % (source_type, totalsourceplaneflux, source_centers[i][0], source_centers[i][1], source_model[i]['es'], source_model[i]['thetas'], source_model[i]['rs'], source_model[i]['ns']  ) ) # sersic/uniform F x y e PA halflightr nothing nS macro/micro
        if source_type == 'uniform':
            # fix 'totalsourceplaneflux' to 1
            totalsourceplaneflux = 1
            f199.write('%s %f %.6f %.6f %f %f %.6f 0 0 macro\n' % (source_type, totalsourceplaneflux, source_centers[i][0], source_centers[i][1], source_model[i]['es'], source_model[i]['thetas'], source_model[i]['rs'] ) ) # sersic/uniform F x y e PA halflightr nothing macro/micro
        #------------------------------------------------------------------------
        f199.write('0 0 0 0 0 0 0 0\n')
        image_names.append('%s%05d_%s.fits' % (base_name, i, reference_band) )
        f199.write('SBmap2 %0.9f %0.9f %d %0.9f %0.9f %d %d %s 3\n' % (-half_frame_size, half_frame_size, Npix, -half_frame_size, half_frame_size, Npix, nover, image_names[-1]) ) # <x lo> <hi> <# steps> <y lo> <hi> <# steps> <Nover> <file> <outtype>
        logging.debug( "The source %d is centered at %s and is properties are %s" % (i+1, source_centers[i] ,str(source_model) ) )

    f199.close()



    if len(source_centers) > 0:
        logging.info( "gravlens is lensing %d finite source(s) (this may take several minutes depending on the resolution)..." % len(source_centers) )
        status = os.system('gravlens sersicinput.txt > /dev/null')
        logging.debug('Executed gravlens to lens the finite sources and returned status %s' % str(status) )
    else:
        logging.info( "There were no sources to be lensed" )

    return image_names, half_frame_size


def identify_images(frame_name, params=[], args={}, preset='sims'):
    """
    Identify the objects on an image.
    
    """

    dict_run_segobj = run_segobj(frame_name, params=[], args={}, preset='') # 'params' are the SExtractor 
    # parameters to output (strings, see SE's default.param). 'args' ({str:str,}) are SE command-line arguments

    objimgname, segimgname, catfilename = dict_run_segobj['OBJECTS'], dict_run_segobj['SEGMENTATION'], dict_run_segobj['CATALOG'] # catfilename is the SE output catalog with the 'params' collumns

    frame_data = pyfits.getdata(frame_name) # to get the header add 'header=True'
    nonzero_frame_data = np.nonzero(frame_data)

    logging.debug('Running elliminate_disconected')
    img_1 = elliminate_disconected(nonzero_frame_data) # get all arguments of this single image
    logging.debug('Finished running elliminate_disconected')
    #img[img_1] = 0

    # make a loop over all images on the frame

    return objimgname, segimgname


def get_src_img_coords(source_name, image_name):
    """
    Get the nonzero points of the source and image plane.
    
    Input:
     - source_name <str> : source plane FITS image
     - image_name  <str> : image plane FITS image

    Output:
     - x_src <ndarray> : 1d array (int) of the source x coordinates (in pixels)
     - y_src <ndarray> : 1d array (int) of the source y coordinates (in pixels)
     - x_img <ndarray> : 1d array (int) of the images x coordinates (in pixels)
     - y_img <ndarray> : 1d array (int) of the images y coordinates (in pixels)
         
    """


    src_data = pyfits.getdata(source_name)
    img_data = pyfits.getdata(image_name)
    src_data = np.nonzero(src_data)
    img_data = np.nonzero(img_data)

    x_src = src_data[1] # o pyfits inverte os eixos
    y_src = src_data[0] # o pyfits inverte os eixos

    x_img = img_data[1] # o pyfits inverte os eixos
    y_img = img_data[0] # o pyfits inverte os eixos        

    return x_src, y_src, x_img, y_img


def pix_2_arcsec(x_src, y_src, x_img, y_img, half_frame_size, dimpix):
    """
    Convert from pixel to arcsec coordinates.
    
    Input:
     - x_src          <ndarray> : 1d array (int) of the source x coordinates (in pixels)
     - y_src          <ndarray> : 1d array (int) of the source y coordinates (in pixels)
     - x_img          <ndarray> : 1d array (int) of the images x coordinates (in pixels)
     - y_img          <ndarray> : 1d array (int) of the images y coordinates (in pixels)
     - half_frame_size  <float> : half the size of the FITS file (in arcsec)
     - dimpix           <float> : size of the pixel, in arcsec

    Output:
     - x_src_arcsec  <ndarray> : 1d array (int) of the source x coordinates (in arcsec)
     - y_src_arcsec  <ndarray> : 1d array (int) of the source y coordinates (in arcsec)
     - x_img_arcsec  <ndarray> : 1d array (int) of the images x coordinates (in arcsec)
     - y_img_arcsec  <ndarray> : 1d array (int) of the images y coordinates (in arcsec)

    """
    
    Npix = int( (2*half_frame_size) / dimpix )
    x_src_arcsec = (2*half_frame_size/(Npix-1))*x_src - half_frame_size -2*half_frame_size/(Npix - 1) # arcsec x coordinate
    y_src_arcsec = (2*half_frame_size/(Npix-1))*y_src - half_frame_size -2*half_frame_size/(Npix - 1) # arcsec y coordinate

    x_img_arcsec = (2*half_frame_size/(Npix-1))*x_img - half_frame_size -2*half_frame_size/(Npix - 1) # arcsec x coordinate
    y_img_arcsec = (2*half_frame_size/(Npix-1))*y_img - half_frame_size -2*half_frame_size/(Npix - 1) # arcsec y coordinate

    return x_src_arcsec, y_src_arcsec, x_img_arcsec, y_img_arcsec, 

#=================================================================================================================
def lens_finite_sources_new(lens_model, mass_scale, model_param_8, model_param_9, model_param_10, dimpix, source_centers, ref_magzpt, reference_band, source_model, galaxy_position=[0,0], e_L=0, theta_L=0, shear=0, theta_shear=0, gravlens_params={}, caustic_CC_file='crit.txt',  gravlens_input_file='gravlens_CC_input.txt', rad_curves_file='lens_curves_rad.dat', tan_curves_file='lens_curves_tan.dat', curves_plot='caustic_CC_curves.png', show_plot=0, write_to_file=0, max_delta_count=20, delta_increment=1.1, grid_factor=5., grid_factor2=3., max_iter_number=20, min_n_lines=200, gridhi1_CC_factor=1.2, accept_res_limit=2E-4, nover_max=3, SE_params={}, SE_args={}, preset='sims', base_name='sbmap'):
    """
    This is a pipeline that ...

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
     - gravlens_params   <dic> : Contains the keys and values of the gravlens configuration , dimpix, source_centers, source_model, 
     - ref_magzpt      <float> :
     - reference_band    <str> :
     - nover_max         <int> :
     - 

    Output:
     - 

    """

    rad_CC_x, rad_CC_y, tan_CC_x, tan_CC_y, rad_caustic_x, rad_caustic_y, tan_caustic_x, tan_caustic_y = run_find_CC(lens_model, mass_scale, model_param_8, model_param_9, model_param_10, galaxy_position, e_L, theta_L, shear, theta_shear, gravlens_params, caustic_CC_file, gravlens_input_file, rad_curves_file, tan_curves_file, curves_plot, show_plot, write_to_file, max_delta_count, delta_increment, grid_factor, grid_factor2, max_iter_number, min_n_lines, gridhi1_CC_factor, accept_res_limit)

    index_CC = np.argmax(tan_CC_x**2 + tan_CC_y**2) # it is not necessary to use the lens center since we want the furthest point anyway
    # gridhi1_CC_factor must be defined/iterated internaly to encompass all images
    gravlens_params['gridhi1'] =  gridhi1_CC_factor * ( (tan_CC_x[index_CC]**2 + tan_CC_y[index_CC]**2)**0.5 )

    # obtain the sources
    source_names = lensing('none', mass_scale, model_param_8, model_param_9, model_param_10, 
                          galaxy_position, e_L, theta_L, shear, theta_shear, gravlens_params, float(dimpix), 
                          source_centers, source_model, ref_magzpt, reference_band, 'src_sbmap', nover_max)[0]
    # obtain the images
    image_names, half_frame_size = lensing(lens_model, mass_scale, model_param_8, model_param_9, model_param_10, 
                          galaxy_position, e_L, theta_L, shear, theta_shear, gravlens_params, float(dimpix), 
                          source_centers, source_model, ref_magzpt, reference_band, base_name, nover_max)

    # loop nos frames (each image plane)
    for n in range( len(image_names) ):


#        # zerar a contagem dos pixels fracos (<0.1max(count))? (attention to uniform sources)
#        counts = frame_data[nonzero_frame_data]; arg_max = np.argmax(counts); max_count_pix = [ nonzero_frame_data[1][arg],nonzero_frame_data[0][arg] ]
#        args_low_count = np.where(frame_data < np.max(counts)/1000.)
#        frame_data[args_low_count] = 0
#        nonzero_frame_data_new = np.nonzero(frame_data)
        objimgname, segimgname = identify_images(image_names[n], SE_params, SE_args, preset='sims')

        # plot sources and images with the caustics and critical curves 
        x_src, y_src, x_img, y_img = get_src_img_coords(source_names[n], image_names[n])
        x_src_arcsec, y_src_arcsec, x_img_arcsec, y_img_arcsec, = pix_2_arcsec(x_src, y_src, x_img, y_img, half_frame_size, dimpix)

        plot_CC(tan_caustic_x, tan_caustic_y, rad_caustic_x, rad_caustic_y, tan_CC_x, tan_CC_y, rad_CC_x, rad_CC_y, 'src_imgs.png', show_plot=0, src_x=x_src_arcsec, src_y=y_src_arcsec, img_x=x_img_arcsec, img_y=y_img_arcsec )
        # ==============================================================

    # extract the image plane images 

    # determine mergers

    # apply measurement methods to individual images







