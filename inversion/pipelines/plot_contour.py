#!/usr/bin/env python
# ==================================
# Author: Gabriel Bartosch Caminha - gbcaminha@gmail.com
# ==================================
"""
Package to plot contour plots
"""
import plot_utils as pu

import argparse
import numpy as np
import matplotlib.pyplot as plt


def main_plot_contour():
    """
    docstring for
    Input
     -
    Output
     -
    #FIXME - finish documentation
    """

    parser = argparse.ArgumentParser(description="Make a contour.", \
                         formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("in_file", \
                    help = "file with the data to be ploted (x_i, y_i, z_i)")

    parser.add_argument("--levels", help = "levels values to be ploted", \
                        type = float, default = None, nargs = "+")

    parser.add_argument("--cols", help = "Columns to be ploted", type = int, \
                        default = [0, 1, 2], nargs = 3)

    parser.add_argument("--labels", help = "Plot or not the labels on the \
                        lines", type = bool, default = False)

    args = parser.parse_args()

    print "plot_contour.py: Plotting file ", args.in_file

    data = np.loadtxt(args.in_file, unpack = True)

    pu.plot_contours(data[args.cols[0]], data[args.cols[1]], \
       data[args.cols[2]], levels = args.levels, plt_labels = args.labels)

    plt.show()

if __name__ == '__main__':
    main_plot_contour()
