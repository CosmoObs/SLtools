#!/usr/bin/env python
# ==================================
# Authors:
# Pedro Ferreira - pferreira@dfte.ufrn.br
# ==================================

""" Reescalonates image to simulate other band """

##@package create_color_img
#
#
# Re-scale counts of an image using a magzpt and magnitude of both reference image and desired image


from __future__ import division
import pyfits



# =================================================================================================
def create_color_img( ref_img, SEXMGZPT_ref, SEXMGZPT_color, out_color, mag_ref, mag_color ):
	"""
	Rescales a given image to simulate other band of an object in the image.

	This function uses the magnitude definition \f$ m = m_0 - 2.5log(F/F_0) \f$ to calculate the
	factor that multiplies every pixel count of the original image:
	\f$ factor = 10^{ (2/5)( SEXMGZPT_{color} - SEXMGZPT_{ref} - mag_{color} + mag_{ref} ) }\f$,
	where SEXMGZPT_color and SEXMGZPT_ref are the zero point magnitudes of the output and input 
	images, while mag_color and mag_ref are the output and input magnitudes of an object of the 
	image. 
	Multiplying an image by this factor simulates an image in the desired band.

	Input:
	 - ref_img    <str> : name of the input fits image
	 - SEXMGZPT_ref <float> : the zero point magnitude of the input image
	 - SEXMGZPT_color <float> : the zero point magnitude of the image in the desired band
	 - out_color <str> : a one string caracter, designating the output band (ex. 'r' or 'g'). 
			     The output image will have this sufix in its name
	 - mag_ref <float> : magnitude of an object in the input image
	 - mag_color <float> : magnitude of the same object in the output image/band

	Output:
	 - <file> : the new image rescaled to another band
	 - <str>  : name of the new band image

	"""

	dados = pyfits.getdata(ref_img)  

	factor = 10**( (2/5)*( SEXMGZPT_color - SEXMGZPT_ref - mag_color + mag_ref ) )

	imgname_out = ''
	for i in range( len(ref_img) - 6): # the last 6 characters corresponds to 'color.fits'
		imgname_out = imgname_out + ref_img[i]
	imgname_out = imgname_out + out_color + '.fits'

	pyfits.writeto('%s' % (imgname_out), dados*factor )

	return imgname_out
# -
