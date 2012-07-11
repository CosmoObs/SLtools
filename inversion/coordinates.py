#!/usr/bin/env python
# ==================================
# Author: Gabriel Bartosch Caminha - gbcaminha@gmail.com
# ==================================

"""
FIXME - Package to deal with astronomiacal coordinates
"""

import pywcs
import pyfits
import numpy as np

import os
import glob

import matplotlib.pylab as plt

def image_limits(hdr):
    """
    Function to get the image limites in wcs
    Input
     - hdr pyfits_hdr : image hdr
    Output
     -
    #FIXME - finish documentation
    """
    im_wcs = pywcs.WCS(hdr)

    pixcrd = np.array([[0,0]], np.float_)

    coor1 = im_wcs.wcs_pix2sky( [[0,0]], 1)
    coor2 = im_wcs.wcs_pix2sky( [[hdr["NAXIS1"], hdr["NAXIS2"]]], 1)

    coor3 = im_wcs.wcs_pix2sky( [[0, hdr["NAXIS2"]]], 1)
    coor4 = im_wcs.wcs_pix2sky( [[hdr["NAXIS1"], 0]], 1)

    c0min = min([coor1[0][0], coor2[0][0], coor3[0][0], coor4[0][0]])
    c0max = max([coor1[0][0], coor2[0][0], coor3[0][0], coor4[0][0]])

    c1min = min([coor1[0][1], coor2[0][1], coor3[0][1], coor4[0][1]])
    c1max = max([coor1[0][1], coor2[0][1], coor3[0][1], coor4[0][1]])

    return [c0min, c0max], [c1min, c1max]

def cs82_limits():
    """
    Function to get the limits from cs82 images
    Input
     -
    Output
     -
    #FIXME - finish documentation
    """
    cs82_dir = "/data1/CS82/CS82_data_products/coadd_V2.7A/images_V2.7A"

    cs82_imgs = glob.glob( cs82_dir + "/*.fit*" )

    limits = []

    for i in range(len(cs82_imgs)):
        #print cs82_imgs[i]
        im_hdu =  pyfits.open(cs82_imgs[i])
        im_hdr = im_hdu[0].header

        coor = image_limits(im_hdr)
        limits.append(coor)
#        ver_x = [coor[0][0], coor[0][0], coor[0][1], coor[0][1], coor[0][0]]
#        ver_y = [coor[1][0], coor[1][1], coor[1][1], coor[1][0], coor[1][0]]
#        plt.plot(ver_x, ver_y)

#    plt.show()
    return limits, cs82_imgs




def main_coordinates():
    """
    Main function
    Input
     -
    Output
     -
    #FIXME - finish documentation
    """
    print "Hello world! from main_coordinates()"
    """
    im_hdu = im_hdu = pyfits.open("/data2/sdss_img_dir/fpC-200006-g5-0629.fit")
    im_hdr = im_hdu[0].header
    print image_limits(im_hdr)

    cs_hdu = pyfits.open("/data1/CS82/CS82_data_products/coadd_V2.7A/images_V2.7A/S82p11p_y.V2.7A.swarp.cut.fits")
    cs_hdr = cs_hdu[0].header
    print image_limits(cs_hdr)
    """

    limits = cs82_limits()
    print len(limits)

if __name__ == '__main__':

    print "Hello world!"
    main_coordinates()
