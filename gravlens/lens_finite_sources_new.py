#!/usr/bin/env python
# ==================================
# Authors:
# Pedro Ferreira - pferreira@dfte.ufrn.br
# ==================================

""" Lenses a finite source with gravlens. """


from __future__ import division
import os
import logging


from sltools.gravlens.lens_parameters_new import lens_parameters_new

##@package lens_finite_sources
# 
# lenses sources (centered in source_centers) with gravlens
# Treats both sersic and uniform sources (uniform sources do not have sersic parameter 'n').
# In the case of uniform sources, the flux is set to unity. Always!
# We increase the resolution of the calculus by increasing Nover (maximum Nover allowed = 3)
# 


#=================================================================================================================
def lens_finite_sources_new(lens_model, mass_scale, model_param_8, model_param_9, model_param_10, dimpix, source_centers, ref_magzpt, reference_band, source_model, galaxy_position=[0,0], e_L=0, theta_L=0, shear=0, theta_shear=0, gravlens_params={}):

    inputlens, setlens, gravlens_params_updated = lens_parameters_new(lens_model, mass_scale, model_param_8, 
                                                  model_param_9, model_param_10, galaxy_position, e_L, theta_L, 
                                                  shear, theta_shear, gravlens_params)

    gridhi1 = gravlens_params_updated['gridhi1']
    half_frame_size = gridhi1
    Npix = int( (2*half_frame_size) / dimpix)
    image_names = []
    f199 = open('sersicinput.txt', 'w')
    f199.write(inputlens)
    for i in range (0,len(source_centers)):
        source_type = source_model[i]['source_type'] # type of the source (can be 'sersic' or 'uniform')
        f199.write('setsource 1\n')
        # Nover increase resolution if the source is smaller than the pixel
        nover = source_model[i]['nover']
        while source_model[i]['rs'] < dimpix/nover and nover < 3:
            nover += 1
        source_model[i]['nover'] = nover
        # -------------------------------------------------------------------
        if source_type == 'sersic':
            # determine 'totalsourceplaneflux', according to the source magnitude
            totalsourceplaneflux = 10**((2/5)*(ref_magzpt - source_model[i][reference_band]))
            f199.write('%s %f %.6f %.6f %f %f %.6f 0 %f macro\n' % (source_type, totalsourceplaneflux, source_centers[i][0], source_centers[i][1], source_model[i]['es'], source_model[i]['thetas'], source_model[i]['rs'], source_model[i]['ns']  ) ) # sersic/uniform F x y e PA halflightr nothing nS macro/micro
        if source_type == 'uniform':
            # fix 'totalsourceplaneflux' to 1
            totalsourceplaneflux = 1
            f199.write('%s %f %.6f %.6f %f %f %.6f 0 0 macro\n' % (source_type, totalsourceplaneflux, source_centers[i][0], source_centers[i][1], source_model[i]['es'], source_model[i]['thetas'], source_model[i]['rs'] ) ) # sersic/uniform F x y e PA halflightr nothing macro/micro
        #------------------------------------------------------------------------
        f199.write('0 0 0 0 0 0 0 0\n')
        f199.write('SBmap2 %0.9f %0.9f %d %0.9f %0.9f %d %d sbmap%05d_%s.fits 3\n' % (-1*half_frame_size, half_frame_size, Npix, -1*half_frame_size, half_frame_size, Npix, nover, i, reference_band) ) # <x lo> <hi> <# steps> <y lo> <hi> <# steps> <Nover> <file> <outtype>	
        image_names.append('sbmap%05d_%s.fits' % (i, reference_band) )
        logging.debug( "The source %d is centered at %s and is properties are %s" % (i+1, source_centers[i] ,str(source_model) ) )

    f199.close()



    if len(source_centers) > 0:
        logging.info( "gravlens is lensing %d finite source(s) (this may take several minutes depending on the resolution)..." % len(source_centers) )
        status = os.system('gravlens sersicinput.txt > /dev/null')
        logging.debug('Executed gravlens to lens the finite sources and returned status %s' % str(status) )
    else:
        logging.info( "There were no sources to be lensed" )
	return image_names








