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

def main_plot_curves_py(plot_show = True):
    """
    Main function

    Input:
     - NONE 
    Output:
     - NONE
    """
    lenstool_file = 'ce.dat'
    gravlens_file = 'crit.dat'
    lenstool_img_file = 'image.all'
    lenstool_src_file = 'source.mul'
    
    is_lenstool_file = os.path.isfile(lenstool_file)
    is_gravlens_file = os.path.isfile(gravlens_file)
    is_lenstool_img = os.path.isfile(lenstool_img_file)
    is_lenstool_src = os.path.isfile(lenstool_src_file)
    
    if is_lenstool_src:
        lenstool_src_pos = [np.loadtxt(lenstool_src_file, usecols = (1, 2))]
        plot_points(source = lenstool_src_pos, plt_show = False, marker='*', \
                    color = 'b')
    if is_lenstool_img:
        lenstool_img_pos = read_images_lenstool(lenstool_img_file)
        plot_points(images = lenstool_img_pos, plt_show = False, marker='o', \
                    color='m')
    if is_gravlens_file:
        gp_obj = GP()
        gp_obj.plot_curves(courves_file = gravlens_file, plt_show = False)
    if is_lenstool_file:
        plot_curves_lenstool(label = 'Lenstool', courves_file = lenstool_file, \
                             plt_show=False)
    if not (is_gravlens_file) and not (is_lenstool_file):
        print 'error: no critical line file found'
    if plot_show:
        plt.show()
if __name__ == '__main__':
    main_plot_curves_py(plot_show = True)
