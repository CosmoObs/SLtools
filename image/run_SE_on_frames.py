
##@package run_SE_on_frames
#
# Runs SExtractor on a list of images
#
#@param input_image_names
#@param sextractor_params
#@return the indices that indicates which images were sextracted succesfully

# ---
import os;
import logging;
# ---
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
	# get SE variables from the config file
	try:
		DEBLEND_MINCONT_merger = float(sextractor_params['deblend_mincont_merger']);
	except KeyError:
		DEBLEND_MINCONT_merger = float(sextractor_params['DEBLEND_MINCONT_merger']);
	try:
		DETECT_MINAREA = float(sextractor_params['detect_minarea']);
	except KeyError:
		DETECT_MINAREA = float(sextractor_params['DETECT_MINAREA']);
	try:
		THRESH_TYPE = sextractor_params['thresh_type'];
	except KeyError:
		THRESH_TYPE = sextractor_params['THRESH_TYPE'];
	try:
		DETECT_THRESH = float(sextractor_params['detect_thresh']);
	except KeyError:
		DETECT_THRESH = float(sextractor_params['DETECT_THRESH']);
	try:
		ANALYSIS_THRESH = float(sextractor_params['analysis_thresh']);
	except KeyError:
		ANALYSIS_THRESH = float(sextractor_params['ANALYSIS_THRESH']);
	try:
		PIXEL_SCALE = float(sextractor_params['pixel_scale']);
	except KeyError:
		PIXEL_SCALE = float(sextractor_params['PIXEL_SCALE']);
	try:
		SEEING_FWHM = float(sextractor_params['seeing_fwhm']);
	except KeyError:
		SEEING_FWHM = float(sextractor_params['SEEING_FWHM']);
	try:
		BACK_TYPE = sextractor_params['back_type'];
	except KeyError:
		BACK_TYPE = sextractor_params['BACK_TYPE'];
	try:
		BACK_VALUE = sextractor_params['back_value'];
	except KeyError:
		BACK_VALUE = sextractor_params['BACK_VALUE'];
	try:
		MEMORY_OBJSTACK = float(sextractor_params['memory_objstack']);
	except KeyError:
		MEMORY_OBJSTACK = float(sextractor_params['MEMORY_OBJSTACK']);
	try:
		MEMORY_PIXSTACK = float(sextractor_params['memory_pixstack']);
	except KeyError:
		MEMORY_PIXSTACK = float(sextractor_params['MEMORY_PIXSTACK']);
	try:
		MEMORY_BUFSIZE = float(sextractor_params['memory_bufsize']);
	except KeyError:
		MEMORY_BUFSIZE = float(sextractor_params['MEMORY_BUFSIZE']);

	#----------------------------------------------
	# now we run SE on the images. The list 'SEstatus' will tell if the run was OK
	img_list = [] # list that contains the images properly generated and SExtracted 
	for i in range (len(input_image_names)): # runs SE on all images (objects e segmentation outputs)
		SEstatus = []
		SEstatus.append( os.system('sex %s -c default.sex -DEBLEND_MINCONT %f -CHECKIMAGE_TYPE OBJECTS -CATALOG_NAME catalogoF%05d.cat -CHECKIMAGE_NAME obj%s -FILTER N -DETECT_MINAREA %d -THRESH_TYPE %s -DETECT_THRESH %f -ANALYSIS_THRESH %f -PIXEL_SCALE %f -SEEING_FWHM %f -CATALOG_TYPE ASCII -BACK_TYPE %s -BACK_VALUE %s -MEMORY_OBJSTACK %d -MEMORY_PIXSTACK %d -MEMORY_BUFSIZE %d > /dev/null 2> /dev/null ' %(input_image_names[i], DEBLEND_MINCONT_merger, i, input_image_names[i], DETECT_MINAREA, THRESH_TYPE, DETECT_THRESH, ANALYSIS_THRESH, PIXEL_SCALE,SEEING_FWHM, BACK_TYPE, BACK_VALUE, MEMORY_OBJSTACK, MEMORY_PIXSTACK, MEMORY_BUFSIZE)) )

		SEstatus.append( os.system('sex %s -c default.sex -DEBLEND_MINCONT %f -CHECKIMAGE_TYPE SEGMENTATION -CATALOG_NAME catalogoF%05d.cat -CHECKIMAGE_NAME segF%s -FILTER N -DETECT_MINAREA %d -THRESH_TYPE %s -DETECT_THRESH %f -ANALYSIS_THRESH %f -PIXEL_SCALE %f -SEEING_FWHM %f -CATALOG_TYPE ASCII -BACK_TYPE %s -BACK_VALUE %s -MEMORY_OBJSTACK %d -MEMORY_PIXSTACK %d -MEMORY_BUFSIZE %d > /dev/null 2> /dev/null ' %(input_image_names[i], DEBLEND_MINCONT_merger, i, input_image_names[i], DETECT_MINAREA, THRESH_TYPE, DETECT_THRESH, ANALYSIS_THRESH, PIXEL_SCALE,SEEING_FWHM, BACK_TYPE, BACK_VALUE, MEMORY_OBJSTACK, MEMORY_PIXSTACK, MEMORY_BUFSIZE)) ) # merger arcs

		if max(SEstatus) != 0:
			logging.warning('BEWARE: SE returned status: %d, %d in image %d\n' % (SEstatus[0], SEstatus[1], i));
		else:
			img_list.append(i)

	return img_list;
# -
