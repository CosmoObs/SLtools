#!/usr/bin/env python
# ==================================
# Author: Gabriel Bartosch Caminha - gbcaminha@gmail.com
# ==================================
"""
Package to plot caustics and critica curves created by  gravlens or lenstool
"""

from lenstool import plot_curves_lenstool, plot_points, read_images_lenstool
from gravlenspipeline import GravlensPipeline as GP

import matplotlib.pyplot as plt
import numpy as np
import os

import argparse

def main_plot_curves_py(plot_show = True):
    """
    Main function

    Input:
     - NONE
    Output:
     - NONE
    """

    parser = argparse.ArgumentParser(description="Plot critical curves and \
             caustics, from lenstool and/or gravlens outputs.", \
             formatter_class = argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--img", help = "Plot or not the images", \
                        default = False, action = "store_true")

    parser.add_argument("--src", help = "Plot or not the source", \
                        default = False, action = "store_true")

    parser.add_argument("--img_file", help = "File with image positions", \
                        default = ["image.all"], type = str, nargs = "+")

    parser.add_argument("--lt_file", help = "File with curves from Lenstool", \
                        default = "ce.dat", type = str)

    parser.add_argument("--src_file", help = "File with source position", \
                        default = "source.mul", type = str)





    args = parser.parse_args()

#    lenstool_file = 'ce.dat'
    lenstool_file = args.lt_file
    gravlens_file = 'crit.dat'
    lenstool_img_file = args.img_file
    lenstool_src_file = args.src_file

    is_lenstool_file = os.path.isfile(lenstool_file)
    is_gravlens_file = os.path.isfile(gravlens_file)

    is_lenstool_img = True
    for i in args.img_file:
        if (not os.path.isfile(i)):
            is_lenstool_img = False
            print "File \"" + i + "\" not found."

    is_lenstool_src = os.path.isfile(lenstool_src_file)

    if is_lenstool_src and args.src:
        lenstool_src_pos = [np.loadtxt(lenstool_src_file, usecols = (1, 2))]
        plot_points(source = lenstool_src_pos, plt_show = False, \
                    color = 'b')
    if is_lenstool_img and args.img:
        for i in args.img_file:
            lenstool_img_pos = read_images_lenstool(i)
            plot_points(images = lenstool_img_pos, plt_show = False, \
                    color= 'm', label = i, marker = "o")
    if is_gravlens_file:
        gp_obj = GP()
        gp_obj.plot_curves(courves_file = gravlens_file, plt_show = False)
    if is_lenstool_file:
        plot_curves_lenstool(label = 'Lenstool', courves_file = lenstool_file, \
                             plt_show=False)
    if not (is_gravlens_file) and not (is_lenstool_file):
        print 'error(main_plot_curves_py): no critical line file found'
    if plot_show:
        plt.show()
if __name__ == '__main__':
    main_plot_curves_py(plot_show = True)
