from __future__ import division
import pyfits

##@package create_color_img
# reescalonate counts of an image using a magzpt and magnitude of both reference image and desired image
#
#@param ref_imgl 
#@param SEXMGZPT_ref
#@param SEXMGZPT_color
#@param out_color
#@param mag_ref
#@param mag_color
#@return name of the color image (ad the image itself)
def create_color_img(ref_img, SEXMGZPT_ref, SEXMGZPT_color, out_color, mag_ref, mag_color):
	dados = pyfits.getdata(ref_img)  
	factor = 10**( (2/5)*( SEXMGZPT_color - SEXMGZPT_ref - mag_color + mag_ref ) )
	imgname_out = ''
	for i in range( len(ref_img) - 6): # the last 6 characters corresponds to 'color.fits'
		imgname_out = imgname_out + ref_img[i]
	imgname_out = imgname_out + out_color + '.fits'
	pyfits.writeto('%s' % (imgname_out), dados*factor )
	return imgname_out


