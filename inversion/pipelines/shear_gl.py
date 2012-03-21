#!/usr/bin/env python

from math import sqrt, fabs

#from functions import generate_fits_image

import numpy as np

import pyfits

import Image

def generate_fits_image(xy_matrx, file_name = 'out.fits'):
    """
    Function to generate a fits image from a matrix

    Input:
     - xy_matrx [[x1y1,x1y2...],[[x2y1,x2y2...]],... ] : file name
    Output:
     - NONE
    """


    hdu = pyfits.PrimaryHDU()
    hdu.data = xy_matrx

    hdu.writeto(file_name, clobber=True)

def creat_sehar_map():
    pot_der_path = 'pot_gl.txt'
    pos_x, pos_y, pot, defx, defy, pot_xx, pot_yy, pot_xy = np.loadtxt( \
                 pot_der_path, unpack=True, usecols={0, 1, 2, 3, 4, 5, 6, 7}, \
                 skiprows=2)


    param1 = 0.25*(pot_xx-pot_yy)**2
    pot_xy_2 = pot_xy**2
    kappa = 0.5*(pot_xx + pot_yy)
    gamma = np.sqrt (param1 + pot_xy_2)#/(1.0-kappa)
    pixn = 201*100
    print pos_x[pixn], pos_y[pixn]
    print pot_xx[pixn], pot_yy[pixn], pot_xy[pixn]
    print kappa[pixn], gamma[pixn]
     

    
    #image = Image.open('exif.jpg')
    #image.show()
    #image.load()
    #xsize, ysize = image.size
    #print xsize, ysize
    
    #r, g, b = image.split()
    #rdata = r.getdata()
    #gdata = g.getdata()
    #bdata = b.getdata()
    pixels = sqrt(len(pot_xx))
    #print pixels
    
    gamma_dat = np.reshape( gamma, (pixels, pixels) )
    gamma_dat = np.transpose(gamma_dat)
    
    gamma1 = np.reshape( 0.5*(pot_xx-pot_yy), (pixels, pixels) )
    gamma1 = np.transpose(gamma1)
    #gamma1 = np.sqrt(gamma1)
    
    pot_xy_2 = np.reshape( pot_xy, (pixels, pixels) )
    pot_xy_2 = np.transpose(pot_xy_2)

   
    tot_mag = np.reshape( ( ((1.0 - pot_xx)*(1.0 - pot_yy) - pot_xy**2)**(-1) ), (pixels, pixels) )
    tot_mag = np.transpose(tot_mag)

    kappa = np.reshape( kappa, (pixels, pixels) )
    kappa = np.transpose(kappa)
    
    generate_fits_image(gamma_dat, file_name='shear_gl.fits')
    generate_fits_image(gamma1, file_name='gamma1_gl.fits')
    generate_fits_image(pot_xy_2, file_name='gamma2_gl.fits')
    generate_fits_image(tot_mag, file_name='ampli_gl.fits')
    generate_fits_image(kappa, file_name='kappa.fits')
    
def compute_diff(image1, image2, file_name):
    file1 = pyfits.open(image1)
    data1 = file1[0].data


    file2 = pyfits.open(image2)
    data2 = file2[0].data
    diff = ( data1 - data2)#*data1[0][0]/data2[0][0])

    generate_fits_image(diff, file_name)

    print data1[0][0], data2[0][0], data1[0][0]/data2[0][0]
    print data1[0][0], data2[0][0]* data1[0][0]/data2[0][0]
    
if __name__ == '__main__':
    creat_sehar_map()
    #compute_diff('ampli.fits', 'ampli_gl.fits', 'dif_ampli.fits')
    #compute_diff('ampli_small.fits', 'ampli_gl_small.fits', 'dif_ampli_small.fits')

    #compute_diff('gamma1.fits', 'gamma1_gl.fits', 'dif_gamma1.fits')
    #compute_diff('gamma2.fits', 'gamma2_gl.fits', 'dif_gamma2.fits')
    
    #compute_diff('mass.fits', 'mass_gl.fits', 'dif_mass.fits')
    #compute_diff('kappa.fits', 'mass.fits', 'dif_mass.fits')



































