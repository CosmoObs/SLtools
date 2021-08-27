#!/usr/bin/python
# =====================================================
# Authors:
# Cristina Furlanetto - furlanetto.cristina@gmail.com
# Eduardo Valadão - eduardovaladao98@gmail.com
# =====================================================
 
import math
import convolve 
import add_noise
import numpy as np 
import paint_arcs_v1_1 as pa
import astropy.io.fits as pyfits
from convert_arcsec_pix import convert_arcsec_pix

def run(params, seeing, background, pix_scale, suffix_name=''):

    hdu = pa.paint_arcs(params, pix_scale)
    
    arc_array = hdu.data
    arc_header = hdu.header
    
    pyfits.writeto("arc"+suffix_name+".fits", arc_array, arc_header, overwrite=True)
     
    dim_x = arc_header["NAXIS1"]
    dim_y = arc_header["NAXIS2"]
    
    # Convolving the arc
    seeing = convert_arcsec_pix(seeing, pix_scale)	
    sigma = seeing/(2.0*math.sqrt(2.0*math.log(2.0)))
    arc_header['APSF'] = (seeing, 'Arcs PSF (arcsec)')	
    
    arc_conv_array = convolve.gauss_convolution_fft(arc_array, 3.0, sigma)
    pyfits.writeto("arc_conv"+suffix_name+".fits", arc_conv_array, arc_header, overwrite=True)
        
    # Adding Poisson noise
    arc_conv_noise_array = add_noise.add_poisson_noise(arc_conv_array)
    pyfits.writeto("arc_conv_noise"+suffix_name+".fits", arc_conv_noise_array, arc_header, overwrite=True)	
    
    # Creating background image with noise
    bkg_array = background*np.ones((dim_y,dim_x)) 
    
    bkg_noise_array = add_noise.add_poisson_noise(bkg_array)
    
    pyfits.writeto("background"+suffix_name+".fits", bkg_noise_array, overwrite=True)
    
    # Adding arc to the background image
    arc_conv_noise_bkg = arc_conv_noise_array + bkg_noise_array
    arc_header['BKG'] = (background, 'Background Mean Value')	
    
    pyfits.writeto("arc_conv_noise_bkg"+suffix_name+".fits", arc_conv_noise_bkg, arc_header, overwrite=True)	


params = [0.6, 4.75, 4.0, 90, 1, 0.12, 0.8, 20, 31.83, 0, 0, 0, 4, 4, 4, 2]
run(params, 0.5, 75, 0.04, "einstein3")

