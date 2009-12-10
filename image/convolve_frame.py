from math import sqrt
import os
import string

##@package convolve_frame
# convolves the specified frames with the specified PSF.
#
#@param arc_files
#@param PSF
#@param dimpix
#@return name of the convolved frame and the convolved image file
def convolve_frame(arc_files, PSF, dimpix):

	# Instrumental parameters conversion (")|(px):
	# The parameters below remain fixed yet
	Arcsec_per_px = dimpix
	FWHM_arcsec = PSF
	##########################################
	# FWHM value in arcsec:
	#    FWHM(px) = FWHM(")/resolucao-espacial("/px)
	#    FWHM_arcsec = input("FWHM(\"): ")
	##########################################
	#FWHM_px = FWHM_arcsec/Arcsec_per_px
	##########################################
	# FWHM = 2*sqrt(2*ln(2))*Sigma
	# Sigma(px) = FWHM(px)/2.35482
	##########################################
	#Sigma_px = FWHM_px/2.35482
	##########################################
	# Half Light Radius:
	#    HLR(px) = FWHM_arcsec/(2*sqrt(Arcsec_per_px^2))
	##########################################
	HLR_px = FWHM_arcsec/(2*sqrt(Arcsec_per_px))
	# Pixel sides dimension (considering square pixels)
	pixel_xy = sqrt(Arcsec_per_px)

	Arcs_convolved = []
	for Arc_filename in arc_files:

		if ( string.find(Arc_filename,'.fits') == -1 ) :
			Arcs_convolved.append(Arc_filename)
 			continue
		# Using IRAF's gauss routine to convolve arc image with a gaussian PSF
		#    cl> gauss image.in image.out sigma=Sigma(px)
		# Usage: gauss('arc.fits','arc_PSF.fits',sigma=Sigma_px)
		Arc_psf_file = string.replace(Arc_filename,'.fits','_psf.fits')
		#gauss(Arc_filename,Arc_psf_file,sigma=Sigma_px)
		Filename_OUT = ''
		for i in range( len(Arc_filename) - 6): # the last 6 characters corresponds to 'color.fits'
			Filename_OUT = Filename_OUT + Arc_filename[i]
		Filename_OUT = Filename_OUT + 'psf_' + Arc_filename[-6] + '.fits' # Arc_filename[-6] is the 6th str from the end (the color)
		os.system('echo "setsource 1" > gravlens.in')
		os.system('echo "sersic 1 0 0 0 0 %s 0 0.5 macro" >> gravlens.in' % (HLR_px))
		os.system('echo "0 0 0 0 0 0 0 0" >> gravlens.in')
		os.system('echo "convolve2 %s 3 %s %s %s 3" >> gravlens.in' % (Arc_filename,pixel_xy,pixel_xy,Filename_OUT))
		os.system('gravlens gravlens.in > /dev/null' )
		#os.remove('gravlens.in')
		Arcs_convolved.append(Filename_OUT)
	return (Arcs_convolved)

