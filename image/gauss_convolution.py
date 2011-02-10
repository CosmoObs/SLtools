#!/usr/bin/env python
#-*- coding:utf-8 -*-

"""This module convolves images with a gaussian kernel"""

##@package gauss_convolution
#
# Module to convolve images with a gaussian kernel of size defined as a 
# multiple of the standard deviation sigma. 
#
# The standard deviation sigma can be obtained from the FWHM using the expression 
# FWHM = 2 \sqrt{2 \ln 2} \sigma 
#
# The convolution can be done using the function convolve_image (normal convolution) 
# or the function fftconvolve_image (convolution using Fast Fourier Tranform). Both 
# functions have three mandatory arguments: the image array to be convolved, the 
# standard deviation sigma of the gaussian kernel and multiple of sigma that defines
# the size of the 2D gaussian kernel. There is also an optional argument that allows 
# convolution with a different sigma in the y direction.
#
# See <A HREF="http://twiki.on.br/bin/view/StrongLensing/CompareConvolutionCodes"> CompareConvolutionCodes</A> 
# wiki page for a quantitative comparison with the results of IRAF and Gravlens convolution tasks. 



import numpy
from scipy import mgrid,signal


def gauss_kernel(n_sig,sigma,sigmay=None):
    """ 
    Creates a normalized 2D gaussian kernel array for convolution.
    The optional keyword argument sigmay allows for a different
    sigma in the y direction.
    the size of the 2D gaussian kernel array is defined as a multiple (n_sig)
    of sigma.
    
    Input:
    - n_sig: multiple of sigma that defines the size of the 2D gaussian kernel
    - sigma: standard deviation sigma of the gaussian kernel
    - sigmay: standard deviation sigma in the y direction of the gaussian kernel

	Output:
    - (ndarray): normalized 2D gaussian kernel array for convolution
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
    """ 
    Convolve the image array with a gaussian kernel of size defined as a 
    multiple (n_sig) of sigma (sig) using FFT. 
    The optional keyword argument sigy allows for a different
    sigma in the y direction.

    Input:
    - im_array: image array to be convolved
    - n_sig: multiple of sigma that defines the size of the 2D gaussian kernel
    - sig: standard deviation sigma of the gaussian kernel
    - sigy: standard deviation sigma in the y direction of the gaussian kernel
    
    Output:
    - (ndarray): image (fft) convolved array 
    
    """

    im_kernel_array = gauss_kernel(n_sig, sig, sigmay=sigy)
    fftconv_image = signal.fftconvolve(im_array,im_kernel_array,mode = 'same')

    return (fftconv_image) 

def convolve_image(im_array, n_sig, sig, sigy=None) :
    """ 
    Convolve the image array with a gaussian kernel of size defined as a 
	multiple (n_sig) of sigma (sig). 
	The optional keyword argument sigy allows for a different
	sigma in the y direction.
        
	Input:    
	- im_array: image array to be convolved
	- n_sig: multiple of sigma that defines the size of the 2D gaussian kernel
	- sig: standard deviation sigma of the gaussian kernel
	- sigy: standard deviation sigma in the y direction of the gaussian kernel
        
	Output:
	- (ndarray): image convolved array 
        
    """

    im_kernel_array = gauss_kernel(n_sig, sig, sigmay=sigy)
    conv_image = signal.convolve(im_array,im_kernel_array,mode = 'same')

    return (conv_image)   

