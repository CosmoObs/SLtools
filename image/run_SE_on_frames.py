import os
import logging

logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(levelname)s %(message)s',
                    filename='AddArcs.log')
                    #filemode='w')


##@package run_SE_on_frames
# runs SExtractor on a list of images
#
#@param input_image_names
#@param sextractor_params
#@return the indices that indicates which images were sextracted succesfully
def run_SE_on_frames(input_image_names, sextractor_params):
	##############################  SExtractor ######################################
	# creates the SExtractor default configuration file
	os.system('sex -dd > default.sex')
	os.system('echo "NUMBER" > default.param')
	os.system('echo "X_IMAGE" >> default.param')
	os.system('echo "Y_IMAGE" >> default.param')
	os.system('echo "X2_IMAGE" >> default.param')
	os.system('echo "Y2_IMAGE" >> default.param')
	os.system('echo "XY_IMAGE" >> default.param')
	#---------------------------------------------------------------------------------
	# gets SE variables from the config file
	DEBLEND_MINCONT_merger = sextractor_params['DEBLEND_MINCONT_merger']
	DETECT_MINAREA = sextractor_params['DETECT_MINAREA']
	THRESH_TYPE = sextractor_params['THRESH_TYPE']
	DETECT_THRESH = sextractor_params['DETECT_THRESH']
	ANALYSIS_THRESH = sextractor_params['ANALYSIS_THRESH']
	PIXEL_SCALE = sextractor_params['PIXEL_SCALE']
	SEEING_FWHM = sextractor_params['SEEING_FWHM']
	BACK_TYPE = sextractor_params['BACK_TYPE']
	BACK_VALUE = sextractor_params['BACK_VALUE']
	MEMORY_OBJSTACK = sextractor_params['MEMORY_OBJSTACK']
	MEMORY_PIXSTACK = sextractor_params['MEMORY_PIXSTACK']
	MEMORY_BUFSIZE = sextractor_params['MEMORY_BUFSIZE']
	#----------------------------------------------
	# now we run SE on the images. The list 'SEstatus' will tell if the run was OK
	img_list = [] # list that contains the images properly generated and SExtracted 
	for i in range (len(input_image_names)): # runs SE on all images (objects e segmentation outputs)
		SEstatus = []
		SEstatus.append( os.system('sex %s -c default.sex -DEBLEND_MINCONT %f -CHECKIMAGE_TYPE OBJECTS -CATALOG_NAME catalogoF%05d.cat -CHECKIMAGE_NAME obj%s -FILTER N -DETECT_MINAREA %d -THRESH_TYPE %s -DETECT_THRESH %f -ANALYSIS_THRESH %f -PIXEL_SCALE %f -SEEING_FWHM %f -CATALOG_TYPE ASCII -BACK_TYPE %s -BACK_VALUE %s -MEMORY_OBJSTACK %d -MEMORY_PIXSTACK %d -MEMORY_BUFSIZE %d > /dev/null 2> /dev/null ' %(input_image_names[i], DEBLEND_MINCONT_merger, i, input_image_names[i], DETECT_MINAREA, THRESH_TYPE, DETECT_THRESH, ANALYSIS_THRESH, PIXEL_SCALE,SEEING_FWHM, BACK_TYPE, BACK_VALUE, MEMORY_OBJSTACK, MEMORY_PIXSTACK, MEMORY_BUFSIZE)) )

		SEstatus.append( os.system('sex %s -c default.sex -DEBLEND_MINCONT %f -CHECKIMAGE_TYPE SEGMENTATION -CATALOG_NAME catalogoF%05d.cat -CHECKIMAGE_NAME segF%s -FILTER N -DETECT_MINAREA %d -THRESH_TYPE %s -DETECT_THRESH %f -ANALYSIS_THRESH %f -PIXEL_SCALE %f -SEEING_FWHM %f -CATALOG_TYPE ASCII -BACK_TYPE %s -BACK_VALUE %s -MEMORY_OBJSTACK %d -MEMORY_PIXSTACK %d -MEMORY_BUFSIZE %d > /dev/null 2> /dev/null ' %(input_image_names[i], DEBLEND_MINCONT_merger, i, input_image_names[i], DETECT_MINAREA, THRESH_TYPE, DETECT_THRESH, ANALYSIS_THRESH, PIXEL_SCALE,SEEING_FWHM, BACK_TYPE, BACK_VALUE, MEMORY_OBJSTACK, MEMORY_PIXSTACK, MEMORY_BUFSIZE)) ) # merger arcs

		if max(SEstatus) != 0:
			print "BEWARE: SE error output number %d in image %d." % (max(SEstatus), i)
			logging.debug('BEWARE: SE returned status: %d, %d, %d in image %d\n' % (SEstatus[0], SEstatus[1], SEstatus[2], i))
		else:
			img_list.append(i)
	return img_list	
