#!/usr/bin/env python
# ==================================
# Author: Gabriel Bartosch Caminha - gbcaminha@gmail.com
# ==================================

"""
Package to compute to projected mass distrubution
"""

import composite_model
import argparse
import astropy.io.fits as pyfits
import math
import functions
import lenstool

import numpy as np
import matplotlib.pylab as plt


def main_compute_projmass(args):
    """
    Main function
    Input
     -
    Output
     -
    #FIXME - finish documentation
    """
    args.radius = args.radius/3600.0

    k_map = pyfits.open(args.kappa_map)
    k_data = k_map[0].data
    k_data_tmp = k_data

    pix_dim = math.fabs(k_map[0].header["CDELT1"])
    pix_unit = k_map[0].header["CUNIT1"]
    if pix_unit != "deg":
        print "Error, pixel unit not in deg"
    shape = k_map[0].data.shape

    x_axis = np.linspace(-(shape[0] - 1.0)/2.0*pix_dim , \
                          (shape[0] - 1.0)/2.0*pix_dim, shape[0])
    y_axis = np.linspace(-(shape[1] - 1.0)/2.0*pix_dim , \
                          (shape[1] - 1.0)/2.0*pix_dim, shape[1])
    proj_mass = 0.0
    for i_x in range(shape[0]):
        for i_y in range(shape[1]):
            if x_axis[i_x]**2.0 + y_axis[i_y]**2.0 < args.radius**2.0:
                #k_data_tmp[i_x][i_y] = 0.0
                proj_mass += k_data_tmp[i_x][i_y]
    print "%fe14 M_sol" % proj_mass

    if args.plot_cont:
        circ = functions.make_circunference(args.radius*3600, 0, 0)
        plt.plot(circ[0], circ[1], "k--", linewidth = 2)
        plt.contour(x_axis*3600.0, y_axis*3600.0, k_data)
        plt.show()

    return proj_mass

if __name__ == '__main__':

    parser = argparse.ArgumentParser(\
        description="Program to compute the projected mass.", \
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--kappa_map", help = "Fits file containing the kappa \
                         map")
    parser.add_argument("--radius", help = "Radius where the projected mass \
                         will be computed in arcseconds", default = 1.0, \
                         type = float)
    parser.add_argument("--plot_cont", help = "If true, plot the mass \
                        contourns", default = False, type = bool)
    args = parser.parse_args()

    if args.kappa_map == None:
        parser.error("Error: you should define kappa_map")


    #print "Hello world!"
#    main_compute_projmass(args)
    lenstool.compute_projmass(args)
