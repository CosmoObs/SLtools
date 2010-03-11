from __future__ import division
import os
import pyfits
import math
import random
import numpy as np
import logging

#import tools;
#from tools.gravlens import lens_parameters
from lens_parameters import lens_parameters # file with the functions

#---------------------------------------------------------------------------------------------------------------------

## get_source_model
# this function executes gravlens to get the file with the critical curves ("crit.txt")
#
#@param inputlens (the base text a input file of gravlens must have, "set gridhi1 = #", "setlens 1 1", etc)
#@return File "crit.txt" with the critical curves positions
def critcurves(gravinput):
	fmag = open('gravlensmagtensor.txt', 'w')
	fmag.write(gravinput)
	fmag.write('plotcrit crit.txt\n' )
	fmag.close()
	os.system('gravlens gravlensmagtensor.txt > /dev/null')

#---------------------------------------------------------------------------------------------------------------------

## find_CC
# this function iteratively runs "critcurves" (iteranting on gridhi1) to find the critical curves
#
#@param lens_model
#@param gravlens_params
#@return "inputlens" and "setlens": two strins that are used to compose gravlens' input files
def find_CC(lens_model, gravlens_params):
	inputlens, setlens = lens_parameters(lens_model, gravlens_params) # inputlens is the gravlens input (sets general parameters and the lens parameters)	# setlens is the gravlens input line that concerns the lens (ex: nfw 1 0 ...)
	logging.debug('Determined the strings that defines the mass model in gravlens through the \'lensparameter\' function')
	#--------------------------------------

	#-----------------------------
	critcurves(inputlens) # gets the critical curves (crit.txt file)
	logging.debug('Got the critical curves (crit.txt file) with function \"critcurves\"')
	#-----------------------------

	# ITERATION ON gridhi1 (TO GET BETTER RESOLUTION ON CC)
	counter = 0
	while os.path.isfile('./crit.txt') == False and counter < 20: # looks for the C.C. file (crit.txt) - if this file does not exist, add a note to arccatalog and quit
		gravlens_params['gridhi1'] /= 5. # gridhi1 /= 5. 
		inputlens, setlens = lens_parameters(lens_model, gravlens_params) # inputlens is the gravlens input (sets general parameters and the lens parameters)
		critcurves(inputlens) # gets the critical curves (crit.txt file)
		counter += 1
	logging.debug( '(Number of iterations on gridhi1 = %d)' % (counter) )

	if os.path.isfile('./crit.txt') == False: # looks for the C.C. file (crit.txt) - if this file does not exist, returns 'False'
		return False

 	if len(open('./crit.txt').readlines() ) < 200: # these correspond to critical cases when the CC have too few points.
		gravlens_params['gridhi1'] /= 3. # gridhi1 /= 3. 
		inputlens, setlens = lens_parameters(lens_model, gravlens_params) # inputlens is the gravlens input (sets general parameters and the lens parameters)
		critcurves(inputlens) # gets the critical curves (crit.txt file)
	logging.debug( 'gridhi1 = %f and number of iterations on gridhi1 = %d' % (gravlens_params['gridhi1'],  counter) )
	return inputlens, setlens

#----------------------------------------------------------------------------------- 


## tang_caustic
# this function uses crit.txt and determines the tangencial caustic
#
#@param inputlens (string that caracterizes the gravlens parameters and the lens parameters)
#@return dictionary {'tan_caustic_x', 'tan_caustic_y', 'tan_CC_x', 'tan_CC_y', 'rad_caustic_x', 'rad_caustic_y'}
def tang_caustic(inputlens):
	fcrit = open('crit.txt', 'r').readlines()
	xcritsrc = []
	ycritsrc = []
	xcritimg = []
	ycritimg = []
	f1 = open('gravlensmagtensor2.txt', 'w')
	f1.write(inputlens) # record the lines relatives to the lens on the new file
	f1.close()
	for i in range (10,len(fcrit)): # because the first 10 lines of this file are only comments
		f1 = open('gravlensmagtensor2.txt', 'a')
		f1.write('magtensor %0.12f %0.12f\n' % ( float(fcrit[i].split()[0]) , float(fcrit[i].split()[1])  ) )
		xcritimg.append(float(fcrit[i].split()[0]))
		ycritimg.append(float(fcrit[i].split()[1]))
		xcritsrc.append(float(fcrit[i].split()[2]))
		ycritsrc.append(float(fcrit[i].split()[3]))
		f1.close()
	os.system('gravlens gravlensmagtensor2.txt > saidamagtensor.txt') # runnins gravlens n times make this function MUCH slower.
	# up to now, we have a file containing all the magnification matrices, besides the lists with all the coodinates of the caustics and C.C.
	f = open('saidamagtensor.txt' , 'r').readlines()
	nlinha = len(f) - len(xcritimg)*6 # I am using this to get the matrix elements without worrying about the number of lines the file has (which depends on
				# the number of parameters of gravlens we modify). 'nlinha' is the number of lines that DON'T contain the output data of magtensor 
	tan_caustic_x = [] # tangential caustic x coordinate
	tan_caustic_y = [] # tangential caustic y coordinate
	rad_caustic_x = [] # radial caustic x coordinate
	rad_caustic_y = [] # radial caustic y coordinate
	tan_CC_x = [] # tangential CC x coordinate
	tan_CC_y = [] # tangential CC y coordinate
	rad_CC_x = [] # radial CC x coordinate
	rad_CC_y = [] # radial CC y coordinate
	for i in range(0,len(xcritimg) ): # computes the magnification of each point of the caustic
		a11 = float(f[nlinha + i*6].split()[0]) # element 11 of the mag tensor
		a12 = float(f[nlinha + i*6].split()[1]) # element 12 (=21) of the mag tensor
		a22 = float(f[nlinha + i*6 + 1].split()[1]) # element 22 of the mag tensor
		lambda_t = ((a11 + a22)/2)/(a11*a22 - a12**2) - ((( (a22 - a11)/2 )**2 + a12**2 )**0.5 )/(abs(a11*a22 - a12**2)) 
		lambda_r = ((a11 + a22)/2)/(a11*a22 - a12**2) + ((( (a22 - a11)/2 )**2 + a12**2 )**0.5 )/(abs(a11*a22 - a12**2)) 
		
		if abs(lambda_t) < 1e-6:
			tan_caustic_x.append(xcritsrc[i])
			tan_caustic_y.append(ycritsrc[i])
			tan_CC_x.append(xcritimg[i])
			tan_CC_y.append(ycritimg[i])
		if abs(lambda_r) < 1e-6:
			rad_caustic_x.append(xcritsrc[i])
			rad_caustic_y.append(ycritsrc[i])
			rad_CC_x.append(xcritimg[i])
			rad_CC_y.append(ycritimg[i])
#	tan_CC_x = np.array(tan_CC_x)
#	tan_CC_y = np.array(tan_CC_y)
#	##############################################################################
#	indice = np.argmax(tan_CC_x**2 + tan_CC_y**2)
#	a = 3.0*( (tan_CC_x[indice]**2 + tan_CC_y[indice]**2)**0.5 ) # define 'a', which will be the region, in arcsec, to be converted to fits
	return {'tan_caustic_x': tan_caustic_x, 'tan_caustic_y': tan_caustic_y, 'tan_CC_x': tan_CC_x, 'tan_CC_y': tan_CC_y, 'rad_caustic_x': rad_caustic_x, 'rad_caustic_y': rad_caustic_y}

#---------------------------------------------------------------------------------------------------------------------


## compute_deformation_rectangle
# determines the 4 cusps of the tangential caustic and determines a rectangle that encloses it. The rectangle vertices are on distances proportional to the cusp distance
#
#@param thetal (inclination of the lens with respect to the vertical, counterclockwise)
#@param tan_caustic_x
#@param tan_caustic_y
#@param source_selector_control_params
#@return deformation_rectangle: [x_vert_1,y_vert_1],[x_vert_2, y_vert_2] (vert_1 refers to the upper right corner, while vert_2 refers to the botton left corner)
def compute_deformation_rectangle(thetal, tan_caustic_x, tan_caustic_y, source_selector_control_params):
	control_rectangle = source_selector_control_params['control_rectangle']
	# converts tan_caustic_x and tan_caustic_y to the system of coordinates in which thetal = 0 (to find the 4 cusps)
	# xcrit0 = cos(theta)*x + sin(theta)*y
	# ycrit0 = -sin(theta)*x + cos(theta)*y
	theta = thetal*(math.pi)/180 # convert thetal to radians
	tan_caustic_x = np.array(tan_caustic_x)
	tan_caustic_y = np.array(tan_caustic_y)
	xcrit0 = math.cos(theta)*tan_caustic_x + math.sin(theta)*tan_caustic_y
	ycrit0 = -math.sin(theta)*tan_caustic_x + math.cos(theta)*tan_caustic_y
	#global cusp1, cusp2, cusp3, cusp4, xlosango, ylosango
	cusp1 = np.argmax(ycrit0)
	cusp2 = np.argmax(xcrit0)
	cusp3 = np.argmin(ycrit0)
	cusp4 = np.argmin(xcrit0)
	xlosango1 = tan_caustic_x[cusp1]*control_rectangle # x coord of point 1 of the losangle
	ylosango1 = tan_caustic_y[cusp1]*control_rectangle # y coord of point 1 of the losangle
	xlosango2 = tan_caustic_x[cusp2]*control_rectangle  # x coord of point 2 of the losangle
	ylosango2 = tan_caustic_y[cusp2]*control_rectangle # y coord of point 2 of the losangle
	xlosango3 = tan_caustic_x[cusp3]*control_rectangle  # x coord of point 3 of the losangle
	ylosango3 = tan_caustic_y[cusp3]*control_rectangle # y coord of point 3 of the losangle
	xlosango4 = tan_caustic_x[cusp4]*control_rectangle # x coord of point 4 of the losangle
	ylosango4 = tan_caustic_y[cusp4]*control_rectangle # y coord of point 4 of the losangle
	xlosango = [xlosango1, xlosango2, xlosango3, xlosango4]
	ylosango = [ylosango1, ylosango2, ylosango3, ylosango4]

	# vert_1 : upper right corner
	# vert_2 : lower left corner
	x_vert_1 = max(xlosango)
	y_vert_1 = max(ylosango)
	x_vert_2 = min(xlosango)
	y_vert_2 = min(ylosango)
	deformation_rectangle = [x_vert_1,y_vert_1],[x_vert_2, y_vert_2]
	return deformation_rectangle

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
def source_positions(source_selector_control_params, deformation_rectangle, nsources, inputlens):
	minimum_distortion = source_selector_control_params['minimum_distortion']
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
		filemag	= open('findimgoutput%05d.txt' % i , 'r').readlines()
		imagens.append(len(filemag)-2)
		image_positions.append([])
		for j in range (2,len(filemag)):	
			f1.write('magtensor %0.6f %0.6f\n' % ( float(filemag[j].split()[0]), float(filemag[j].split()[1])  ) )
			image_positions[i].append( [float(filemag[j].split()[0]), float(filemag[j].split()[1])] )
	# if gravlens does not finf any images, the grid is not big enough
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
			mu_t = ( ((a11 + a22)/2)/(a11*a22 - a12**2) - ((( (a22 - a11)/2 )**2 + a12**2 )**0.5 )/(abs(a11*a22 - a12**2)) )**(-1)
			mu_r = ( ((a11 + a22)/2)/(a11*a22 - a12**2) + ((( (a22 - a11)/2 )**2 + a12**2 )**0.5 )/(abs(a11*a22 - a12**2)) )**(-1)
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
def select_source_positions(lens_model, gravlens_params, source_selector_control_params, src_density_or_number): # source_selector_control_params = [minimum_distortion, control_rectangle]
	# find critical curves (returns 'False' if can't find them)
	out_find_CC = find_CC(lens_model, gravlens_params)
	if out_find_CC == False:
		return False
	logging.debug('Obtained critical curves')
	inputlens, setlens = out_find_CC
	#----------------------------------------------------------

	#----------------------------------------------------------------------------------------------
	# we now get the outputs of the function tang_caustic (radial and tangential caustic and CC)
	outtang_caustic = tang_caustic(inputlens)
	tan_caustic_x = outtang_caustic['tan_caustic_x'] # tangencial caustic x coordinate
	tan_caustic_y = outtang_caustic['tan_caustic_y'] # tangencial caustic y coordinate
	tan_CC_x = outtang_caustic['tan_CC_x'] # tangential CC x coordinate
	tan_CC_y = outtang_caustic['tan_CC_y'] # tangential CC y coordinate
	rad_caustic_x = outtang_caustic['rad_caustic_x'] # radial caustic x coordinate
	rad_caustic_y = outtang_caustic['rad_caustic_y'] # radial caustic y coordinate
	logging.debug('Obtained the radial and tangential critical curves through the function \"tang_caustic\"')
	#----------------------------------------------------------------------------------------------

	#------------ call compute_deformation_rectangle ---------------------------
	deformation_rectangle = compute_deformation_rectangle(lens_model['thetal'], tan_caustic_x, tan_caustic_y, source_selector_control_params)
	logging.debug('Obtained the deformation rectangle: %s' % str(deformation_rectangle) )
	#---------------------------------------------------------------------------

	#-----------------------------------------------------------------------------------------------------------
	x_vert_1, y_vert_1 = deformation_rectangle[0]
	x_vert_2, y_vert_2 = deformation_rectangle[1]
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
	image_plane_factor = source_selector_control_params['image_plane_factor']
	gravlens_params['gridhi1'] =  image_plane_factor * ( (tan_CC_x[index]**2 + tan_CC_y[index]**2)**0.5 )
	logging.debug( 'gridhi1 = %f' % gravlens_params['gridhi1'] )
	inputlens, setlens = lens_parameters(lens_model, gravlens_params)
	#-----------------------------------------------------------------------------------------------------------
	source_positions_output = source_positions(source_selector_control_params, deformation_rectangle, nsources, inputlens)
	while source_positions_output == False:
		gravlens_params['gridhi1'] *=  1.15
		logging.debug( 'loop in source_positions: gridhi1 = %f' % gravlens_params['gridhi1'] )
		inputlens, setlens = lens_parameters(lens_model, gravlens_params)		
		source_positions_output = source_positions(source_selector_control_params, deformation_rectangle, nsources, inputlens)
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

















