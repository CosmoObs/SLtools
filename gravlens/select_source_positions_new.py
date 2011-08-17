#!/usr/bin/env python
# ==================================
# Authors:
# Pedro Ferreira - pferreira@dfte.ufrn.br
# ==================================

""" Select the position of sources near the tangential caustic """


##@package select_source_positions
#
#
#  Select the position of sources near the tangential caustic (documentation in the old style)


from __future__ import division
import os
import pyfits
import math
import random
import numpy as np
import logging


from sltools.gravlens.find_CC_new import run_find_CC
from sltools.gravlens.lens_parameters_new import lens_parameters_new




#=======================================================================================================
def compute_deformation_rectangle(thetal, tan_caustic_x, tan_caustic_y, control_rectangle):
    """
    Determines the 4 cusps of the tangential caustic and determines a rectangle that encloses it. 

    The rectangle vertices are on distances proportional to the cusp distance.

    Input:
     - thetal             <float> : inclination of the lens with respect to the vertical
                                   (counterclockwise)
     - tan_caustic_x       <list> : list of x coordinates of the tangential caustic
     - tan_caustic_y       <list> : list of y coordinates of the tangential caustic
     - control_rectangle  <float> : determines how many times the vertices of the losangle 
                                    (inscribed in the rectangle) will be away from the cusps.

    Output:
     - deformation_rectangle_pt1  <list> : pair (x,y) of upper right corner of the rectangle
     - deformation_rectangle_pt2  <list> : pair (x,y) of lower left corner of the rectangle

    """

    # converts tan_caustic_x and tan_caustic_y to the system of coordinates in which thetal = 0 (to find the 4 cusps)
    # xcrit0 = cos(theta)*x + sin(theta)*y
    # ycrit0 = -sin(theta)*x + cos(theta)*y
    theta = thetal*(math.pi)/180. # convert thetal to radians
    tan_caustic_x = np.array(tan_caustic_x)
    tan_caustic_y = np.array(tan_caustic_y)
    xcrit0 =  math.cos(theta)*tan_caustic_x + math.sin(theta)*tan_caustic_y
    ycrit0 = -math.sin(theta)*tan_caustic_x + math.cos(theta)*tan_caustic_y
    #global cusp1, cusp2, cusp3, cusp4, xlosango, ylosango
    cusp1 = np.argmax(ycrit0)
    cusp2 = np.argmax(xcrit0)
    cusp3 = np.argmin(ycrit0)
    cusp4 = np.argmin(xcrit0)
    xlosango1 = tan_caustic_x[cusp1]*control_rectangle # x coord of point 1 of the losangle
    ylosango1 = tan_caustic_y[cusp1]*control_rectangle # y coord of point 1 of the losangle
    xlosango2 = tan_caustic_x[cusp2]*control_rectangle # x coord of point 2 of the losangle
    ylosango2 = tan_caustic_y[cusp2]*control_rectangle # y coord of point 2 of the losangle
    xlosango3 = tan_caustic_x[cusp3]*control_rectangle # x coord of point 3 of the losangle
    ylosango3 = tan_caustic_y[cusp3]*control_rectangle # y coord of point 3 of the losangle
    xlosango4 = tan_caustic_x[cusp4]*control_rectangle # x coord of point 4 of the losangle
    ylosango4 = tan_caustic_y[cusp4]*control_rectangle # y coord of point 4 of the losangle
    xlosango = [xlosango1, xlosango2, xlosango3, xlosango4]
    ylosango = [ylosango1, ylosango2, ylosango3, ylosango4]

    return xlosango, ylosango

#---------------------------------------------------------------------------------------------------------------------

## nsource
# determines the number of point sources to be generated using a poisson of rectangle_area*source_surface_density.
#
#@param rectangle_area
#@param source_surface_density
#@return number of point sources to be generated
def nsource(rectangle_area, source_surface_density):

    nsources = np.random.poisson(rectangle_area*source_surface_density )
    return nsources

#---------------------------------------------------------------------------------------------------------------------    

## source_positions
# determines the source centers that are candidates to generate an arc
#
#@param source_selector_control_params
#@param deformation_rectangle
#@param nsources
#@param inputlens
#@return source_centers [former potarcxy], image_centers, image_distortions    
def source_positions(minimum_distortion, deformation_rectangle, nsources, inputlens):

    x_vert_1, y_vert_1 = deformation_rectangle[0]
    x_vert_2, y_vert_2 = deformation_rectangle[1]
    #-----------------------------------------------------------------------------------------------------------
    #GENERATE RANDOM POINTS NEAR THE CAUSTIC -------------------------------------------------------------------
    #global usq, vsq, usq2, vsq2
    usq = []
    vsq = []
    for i in range(0,nsources):
        xrand =  x_vert_2 + ( x_vert_1 - x_vert_2 ) * random.random() # x coordinate inside the square
        yrand =  y_vert_2 + ( y_vert_1 - y_vert_2 ) * random.random() # y coordenada inside the square
        usq.append(xrand) # x coord of a point inside the square, but outside the losangle
        vsq.append(yrand) # y coord of a point inside the square, but outside the losangle
    #----------------------------------------------------------------------------------------------------------------------
    ########### DETERMINE WHICH CENTERS HAS mu_t/mu_r > minimum_distortion IN ORDER TO GENERATE SERSIC SOURCES IN THESE CENTERS ############
    #----------------------------------------------------------------------------------------------------------------------
    f = open('findimginput.txt', 'w')
    f.write(inputlens)
    for i in range (0,len(usq)): # determine the positions of the images of the source
        #f = open('findimginput.txt', 'a')
        f.write('findimg %0.6f %0.6f findimgoutput%05d.txt\n'  % ( usq[i],vsq[i], i ))
    f.close()
    # runs gravlens to obtain the images of the positions inside the rectangle (usq,vsq)
    os.system('gravlens findimginput.txt > /dev/null')
    f1 = open('gravlensmagtensorcorte.txt', 'w')
    f1.write(inputlens)
    imagens = [] # contains the number of images of the source j
    image_positions = [] # image_positions[i][j] is the (x,y) position of the image j of the source i
    for i in range (0,len(usq)):
        filemag = open('findimgoutput%05d.txt' % i , 'r').readlines()
        imagens.append(len(filemag)-2)
        image_positions.append([])
        for j in range (2,len(filemag)):    
            f1.write('magtensor %0.6f %0.6f\n' % ( float(filemag[j].split()[0]), float(filemag[j].split()[1])  ) )
            image_positions[i].append( [float(filemag[j].split()[0]), float(filemag[j].split()[1])] )
    # if gravlens does not find any images, the grid is not big enough
    if 0 in imagens and len(imagens) != 0:
        logging.debug( 'There are %d point images outside the grid' % imagens.count(0) )
        return False
    f1.close()
    # obtain with the magtensor all the magnifications
    os.system('gravlens gravlensmagtensorcorte.txt > saidamagtensorcorte.txt')
    f = open('saidamagtensorcorte.txt' , 'r').readlines()
    nlinha = len(f) - 6*(np.sum(imagens)) # I am using this to get the matrix elements without worrying about the number of lines the file has (which depends on the number of parameters of gravlens we modify). 'nlinha' is the number of lines that DON'T contain the output data of magtensor
    source_centers = []
    image_centers = [] # image_positions[i][j] is the (x,y) position of the image j of the source i
    image_distortions = [] # image_distortions[i][j] is the (mu_t,mu_r) distortion of the image j of the source i
    usq2 = []
    vsq2 = []
    contj = 0
    for i in range(0,len(imagens) ): # calculates the magnifications at each point
        mutemp = 0    
        distortion_temp = []
        for j in range (0,imagens[i]): # note that the index 'i' is the same as of usq
            a11 = float(f[nlinha + (contj)*6].split()[0]) # element 11 of the mag tensor
            a12 = float(f[nlinha + (contj)*6].split()[1]) # element 12 (=21) of the mag tensor
            a22 = float(f[nlinha + (contj)*6 + 1].split()[1]) # element 22 of the mag tensor
            mu_t = ( ((a11 + a22)/2.)/(a11*a22 - a12**2) - ((( (a22 - a11)/2. )**2 + a12**2 )**0.5 )/(abs(a11*a22 - a12**2)) )**(-1)
            mu_r = ( ((a11 + a22)/2.)/(a11*a22 - a12**2) + ((( (a22 - a11)/2. )**2 + a12**2 )**0.5 )/(abs(a11*a22 - a12**2)) )**(-1)
            distortion_temp.append( [mu_t,mu_r] )
            contj = contj + 1
            if abs(mu_t/mu_r) > mutemp:
                mutemp = abs(mu_t/mu_r)
        if mutemp > minimum_distortion:
            source_centers.append([usq[i], vsq[i]])
            image_centers.append( image_positions[i] )
            image_distortions.append( distortion_temp )
        else:
            usq2.append(usq[i])
            vsq2.append(vsq[i])
    os.system('rm -f findimgoutput*.txt')
    return source_centers, image_centers, image_distortions  #, usq2, vsq2, cusp1, cusp2, cusp3, cusp4, xlosango, ylosango, usq, vsq

#---------------------------------------------------------------------------------------------------------------------

##@package select_source_positions [formerly SourceSelector]
# Identifies source_positions close to caustics (i.e. that could produce arcs)
#
# Given the lens parameters and source surface density, identifies the source_positions that have a local distortion above minimum_distortion 
# 
#@param lens_model (kappas, xsl,el,thetal) currently a dictionary with the lens model parameters 
#@param gravlens_params (set of gravlens parameters: gridhi1, maxlev, etc)
#@param source_selector_control_params (minimum_distortion, control_rectangle) 
#@param minimum_distortion : minimum ratio of the tangential and radial magnifications to an image be considered an arc 
#@param control_rectangle : control parameter to define distance to cusps 
#@param src_density_or_number
#@return source_centers [former potarcxy] or "False" (in cases where no critical curves are found), and image_centers, image_distortions
def select_source_positions(lens_model, mass_scale, model_param_8, model_param_9, model_param_10, galaxy_position=[0.,0.], e_L=0, theta_L=0, shear=0, theta_shear=0, gravlens_params={}, src_density_or_number=1, minimum_distortion=0., control_rectangle=2., caustic_CC_file='crit.txt', gravlens_input_file='gravlens_CC_input.txt', rad_curves_file='lens_curves_rad.dat', tan_curves_file='lens_curves_tan.dat', curves_plot='crit-caust_curves.png', show_plot=0, write_to_file=0, max_delta_count=20, delta_increment=1.1, grid_factor=5., grid_factor2=3., max_iter_number=20, min_n_lines=200, gridhi1_CC_factor=2., accept_res_limit=2E-4):
    """
    
    """

    # find and separate critical curves
    out_run_find_CC = run_find_CC(lens_model, mass_scale, model_param_8, model_param_9, model_param_10, galaxy_position, e_L, theta_L, shear, theta_shear, gravlens_params, caustic_CC_file, gravlens_input_file, rad_curves_file, tan_curves_file, curves_plot, show_plot, write_to_file, max_delta_count, delta_increment, grid_factor, grid_factor2, max_iter_number, min_n_lines, gridhi1_CC_factor, accept_res_limit)

    logging.info('Ran run_find_CC and obtained and separated the caustic and critical curves')

    rad_CC_x, rad_CC_y, tan_CC_x, tan_CC_y, rad_caustic_x, rad_caustic_y, tan_caustic_x, tan_caustic_y = out_run_find_CC
    #----------------------------------------------------------

    #------------ call compute_deformation_rectangle ---------------------------
    xlosango, ylosango = compute_deformation_rectangle(float(theta_L), tan_caustic_x, tan_caustic_y, control_rectangle)
    #---------------------------------------------------------------------------

    #-----------------------------------------------------------------------------------------------------------

    # vert_1 : upper right corner
    # vert_2 : lower left corner
    x_vert_1, y_vert_1 = max(xlosango), max(ylosango)
    x_vert_2, y_vert_2 = min(xlosango), min(ylosango)
    deformation_rectangle = [x_vert_1,y_vert_1],[x_vert_2, y_vert_2]
    logging.debug('Obtained the deformation rectangle: %s' % str(deformation_rectangle) )

    rectangle_area = ( x_vert_1 - x_vert_2 )*( y_vert_1 - y_vert_2 ) # deltaX*deltaY
    # if src_density_or_number is a float, it is a density, if it is a integer it is the number of sources
    if type(src_density_or_number) == int:
        nsources = src_density_or_number
    if type(src_density_or_number) == float:    
        nsources = nsource(rectangle_area, src_density_or_number) # generates 'nsource' sources per arcsec^2
    logging.debug('Number of point sources generated = %d' % nsources)
    #-----------------------------------------------------------------------------------------------------------
    # redefine gridhi1
    tan_CC_x = np.array(tan_CC_x)
    tan_CC_y = np.array(tan_CC_y)
    index = np.argmax(tan_CC_x**2 + tan_CC_y**2)
    image_plane_factor = float(gridhi1_CC_factor);
    gravlens_params['gridhi1'] =  image_plane_factor * ( (tan_CC_x[index]**2 + tan_CC_y[index]**2)**0.5 )
    inputlens, setlens, gravlens_params = lens_parameters_new(lens_model, mass_scale, model_param_8, model_param_9, model_param_10, galaxy_position, e_L, theta_L, shear, theta_shear, gravlens_params)
    #-----------------------------------------------------------------------------------------------------------
    source_positions_output = source_positions(minimum_distortion, deformation_rectangle, nsources, inputlens)
    while source_positions_output == False:
        gravlens_params['gridhi1'] = float( gravlens_params['gridhi1'] ) * 1.15;
        logging.debug( 'loop in source_positions: gridhi1 = %f' % float( gravlens_params['gridhi1'] ) )
        inputlens, setlens, gravlens_params = lens_parameters_new(lens_model, mass_scale, model_param_8, model_param_9, model_param_10, galaxy_position, e_L, theta_L, shear, theta_shear, gravlens_params)       
        source_positions_output = source_positions(minimum_distortion, deformation_rectangle, nsources, inputlens)
    logging.debug( 'gridhi1 = %s' % float( gravlens_params['gridhi1'] ) )
    source_centers = source_positions_output[0]
    image_centers = source_positions_output[1]
    image_distortions = source_positions_output[2]

    logging.debug('Selected positions for finite sources: %s' % str(source_centers) )
    logging.debug('Images of the point sources: %s' % str(image_centers) )
    logging.debug('Corresponding mu_t and mu_r for each image: %s' % str(image_distortions) )

    #from arcgeneratormodules import grafico
    #grafico(tan_caustic_x,tan_caustic_y, rad_caustic_x,rad_caustic_y, usq2, vsq2, source_centers, cusp1, cusp2, cusp3, cusp4, xlosango, ylosango, usq, vsq, halo_id+'1')
    #grafico(tan_caustic_x,tan_caustic_y, rad_caustic_x,rad_caustic_y, [0], [0], source_centers, cusp1, cusp2, cusp3, cusp4, xlosango, ylosango, [0], [0],halo_id+'2')
    return source_centers, image_centers, image_distortions 
