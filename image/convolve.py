#!/usr/bin/env python
#-*- coding:utf-8 -*-
# ===================================================
# Authors:
# Cristina Furlanetto - furlanetto.cristina@gmail.com
# ===================================================



"""This module convolves images with a gaussian kernel"""

##@package convolve
#
# Module to convolve images with an isotropic gaussian kernel of size defined as a 
# multiple of the standard deviation sigma. 
#
# The standard deviation sigma can be obtained from the FWHM using the expression 
# \f$FWHM = 2 \sqrt{2 \ln 2} \sigma \f$
#
# The convolution can be done using the function gauss_convolution (normal convolution) 
# or the function gauss_convolution_fft (convolution using Fast Fourier Tranform). Both 
# functions have three mandatory arguments: the image array to be convolved, the 
# standard deviation sigma of the gaussian kernel and multiple of sigma that defines
# the size of the 2D gaussian kernel.
#
# See <A HREF="http://twiki.linea.gov.br/bin/view/StrongLensing/CompareConvolutionCodes"> CompareConvolutionCodes</A> 
# wiki page for a quantitative comparison with the results of IRAF and Gravlens convolution tasks. 



import numpy
from scipy import mgrid,signal


def gauss_kernel(n_sig,sigma):
    """ 
    Creates a normalized 2D isotropic gaussian kernel array for convolution.
    The size of the 2D gaussian kernel array is defined as a multiple (n_sig)
    of sigma.
    
    Input:
    - n_sig <float> : multiple of sigma that defines the size of the 2D gaussian kernel (usually n_sig=4)
    - sigma <float>: standard deviation sigma of the gaussian kernel

	Output:
    - <ndarray>: normalized 2D gaussian kernel array for convolution
    """

    x_length = int(n_sig * sigma + 0.5) #Add 0.5 to approximate to nearest integer
    y_length = x_length
        
       
    x, y = mgrid[-x_length:x_length+1, -y_length:y_length+1]
    g = numpy.exp(-(x**2/(2*(float(sigma)**2))+y**2/(2*(float(sigma)**2))))
    return g / g.sum()


def gauss_convolution_fft(im_array, n_sig, sig) :
    """ 
    Convolves the image array with an isotropic gaussian kernel of size defined as a 
    multiple (n_sig) of sigma (sig) using FFT. 

    Input:
    - im_array <ndarray>: image array to be convolved
    - n_sig <float>: multiple of sigma that defines the size of the 2D gaussian kernel (usually n_sig=4)
    - sig <float>: standard deviation sigma of the gaussian kernel
    
    Output:
    - <ndarray>: image (fft) convolved array 
    
    """

    im_kernel_array = gauss_kernel(n_sig, sig)
    fftconv_image = signal.fftconvolve(im_array,im_kernel_array,mode = 'same')

    return (fftconv_image) 

def gauss_convolution(im_array, n_sig, sig) :
    """ 
    Convolves the image array with an isotropic gaussian kernel of size defined as a 
	multiple (n_sig) of sigma (sig). 
	
	Input:    
	- im_array <ndarray>: image array to be convolved
	- n_sig <float>: multiple of sigma that defines the size of the 2D gaussian kernel (usually n_sig=4)
	- sig <float>: standard deviation sigma of the gaussian kernel
        
	Output:
	- <ndarray>: image convolved array 
        
    """

    im_kernel_array = gauss_kernel(n_sig, sig)
    conv_image = signal.convolve(im_array,im_kernel_array,mode = 'same')

    return (conv_image)   


def moffat_kernel(n_rs,beta,r_s):
    """ 
    Creates a normalized 2D Moffat kernel array for convolution.
    The size of the 2D kernel array is defined as a multiple (n_rs)
    of the scale radius (r_s).
	    
    Input:
    - n_rs <float>: multiple of r_s that defines the size of the 2D Moffat kernel (usually n_rs=4)
    - beta <float>: profile slope
    - r_s <float>: scale radius
	
    Output:
    - <ndarray>: normalized 2D kernel array for convolution
    """

    x_length = int(n_rs * r_s + 0.5) #Add 0.5 to approximate to nearest integer
    y_length = x_length
        

    x, y = mgrid[-x_length:x_length+1, -y_length:y_length+1]
	
    m = 1. /((1+(x**2+y**2)/r_s**2)**beta)
		

    return m / m.sum()	
	
	
	
def moffat_convolution_fft(im_array,n_rs,beta,r_s) :
    """ 
    Convolves the image array with an Moffat kernel of size defined as a 
    multiple (n_rs) of r_s using FFT. 

    Input:
    - im_array <ndarray>: image array to be convolved
    - n_rs <float>: multiple of r_s that defines the size of the 2D Moffat kernel (usually n_rs=4)
    - beta <float>: profile slope
    - r_s <float>: scale radius
    
    Output:
    - <ndarray>: image (fft) convolved array 
    
    """

    im_kernel_array = moffat_kernel(n_rs,beta,r_s)
    fftconv_image = signal.fftconvolve(im_array,im_kernel_array,mode = 'same')

    return (fftconv_image) 



def moffat_convolution(im_array,n_rs,beta,r_s) :
    """ 
    Convolves the image array with an Moffat kernel of size defined as a 
    multiple (n_rs) of r_s. 
	
    Input:    
    - im_array <ndarray>: image array to be convolved
    - n_rs <float>: multiple of r_s that defines the size of the 2D Moffat kernel (usually n_rs=4)
    - beta <float>: profile slope
    - r_s <float>: scale radius	
	
    Output:
    - <ndarray>: image convolved array 
        
    """

    im_kernel_array = gauss_kernel(n_rs,beta,r_s)
    conv_image = signal.convolve(im_array,im_kernel_array,mode = 'same')

    return (conv_image)   


