#!/usr/bin/env python
#-*- coding:utf-8 -*-
"""@file

gauss_convolution module

"""


"""@package gauss_convolution

Module to convolve images with a gaussian kernel of size defined as a 
multiple of sigma.

See <A HREF="http://twiki.on.br/bin/view/StrongLensing/CompareConvolutionCodes"> CompareConvolutionCodes</A> wiki page for a quantitative comparison with the results of IRAF and Gravlens convolution tasks. 

"""


import numpy
from scipy import mgrid,signal


def gauss_kernel(n_sig,sigma,sigmay=None):
   	""" Creates a normalized 2D gaussian kernel array for convolution.
	The optional keyword argument sigmay allows for a different
        sigma in the y direction.
	The size of the 2D gaussian kernel array is defined as a multiple (n_sig)
	of sigma.
	
	@param n_sig  multiple of sigma (integer) that defines the size of the 2D gaussian kernel
	@param sigma  sigma of the gaussian kernel
	@param sigmay  sigma in the y direction of the gaussian kernel
	
	@return normalized 2D gaussian kernel array for convolution
    	"""

    x_length = int(n_sig * sigma + 0.5) #Add 0.5 to approximate to nearest integer
    
    if not sigmay:
        sigmay = sigma
	y_length = x_length
    else:
        y_length = int(n_sig * sigmay + 0.5)        
       
    x, y = mgrid[-x_length:x_length+1, -y_length:y_length+1]
    g = numpy.exp(-(x**2/(2*(float(sigma)**2))+y**2/(2*(float(sigmay)**2))))
    return g / g.sum()


def fftconvolve_image(im_array, n_sig, sig, sigy=None) :
    	""" Convolve the image array with a gaussian kernel of size defined as a 
        multiple (n_sig) of sigma (sig) using FFT. 
	The optional keyword argument sigy allows for a different
        sigma in the y direction.
        
        
        @param im_array  image array to be convolved
        @param n_sig  multiple of sigma (integer) that defines the size of the 2D gaussian kernel
        @param sig sigma of the gaussian kernel
        @param sigy  sigma in the y direction of the gaussian kernel
        
        @return image (fft) convolved array 
    	"""

    im_kernel_array = gauss_kernel(n_sig, sig, sigmay=sigy)
    fftconv_image = signal.fftconvolve(im_array,im_kernel_array,mode = 'same')

    return (fftconv_image) 




def convolve_image(im_array, n_sig, sig, sigy=None) :
    	""" Convolve the image array with a gaussian kernel of size defined as a 
        multiple (n_sig) of sigma (sig). 
	The optional keyword argument sigy allows for a different
        sigma in the y direction.
        
        
        @param im_array  image array to be convolved
        @param n_sig  multiple of sigma (integer) that defines the size of the 2D gaussian kernel
        @param sig sigma of the gaussian kernel
        @param sigy  sigma in the y direction of the gaussian kernel
        
        @return image convolved array 
        
    	"""

    im_kernel_array = gauss_kernel(n_sig, sig, sigmay=sigy)
    conv_image = signal.convolve(im_array,im_kernel_array,mode = 'same')

    return (conv_image)   

