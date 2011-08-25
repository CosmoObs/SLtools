def convert_arcsec_pix(x,pix_scale):
	'''
	Convert the units from arcsec to pixels.
	
	Input:
	- x <float>: value in arcsec
	- pix_scale <float>: pixel scale (arcsec/pix)
	
	Output:
	- <float>: value in pixels
	'''
	return x / pix_scale

def convert_pix_arcsec(x,pix_scale):
	'''
	Convert the units from pixels to arcsec.
	
	Input:
	- x  <float>: value in pixel
	- pix_scale <float>: pixel scale (arcsec/pix)
	
	Output:
	- <float>: value in arcsec
	'''
	return x * pix_scale
