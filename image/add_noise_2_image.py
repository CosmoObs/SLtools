import pyfits
from numpy import zeros, where
from pylab import poisson
# poissonian noise addition
def add_noise_2_image(Arc_psf_file):
	Arcs_noise = []
	for Arc_filename in Arc_psf_file:
		# Unfortunately the convolution process outputs to a file..
		# Lets open it to an array..
		ArcImg_array,ArcImg_header = pyfits.getdata(Arc_filename,header=True)
		###########################################################################
		# Here I am treating a strange "behauevoiubfd" of the Arc images... !-P (negative counts)
		Arc_tmp = zeros(ArcImg_array.shape)
		Arc_tmp[where(ArcImg_array>0)] = ArcImg_array[where(ArcImg_array>0)]
		ArcImg_array = poisson(Arc_tmp)
		###########################################################################

		# Re-sort pixel intensities of arc image in a poissonian fashion
		#    >>> arco = poisson(arco)
		ArcImg_array = poisson(ArcImg_array)

		#ArcImg_header.header.add_history('PSF convolved and noise added.')

		# Outputs Image to fits files
		Filename_OUT = ''
		for i in range( len(Arc_filename) - 6): # the last 6 characters corresponds to 'color.fits'
			Filename_OUT = Filename_OUT + Arc_filename[i]
		Filename_OUT = Filename_OUT + 'ns_' + Arc_filename[-6] + '.fits' # Arc_filename[-6] is the 6th str from the end (the colour)

		#Filename_OUT = string.replace(Arc_filename,'.fits','_psf_nz.fits')
		pyfits.writeto(Filename_OUT,ArcImg_array,ArcImg_header)
		
		Arcs_noise.append(Filename_OUT)

	return (Arcs_noise)
