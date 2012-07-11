#!/usr/bin/env python
# ==================================
# Author: Gabriel Bartosch Caminha - gbcaminha@gmail.com
# ==================================

"""
FIXME - package description
"""

import numpy as np
import parser as ps
import functions as fc

import argparse


if __name__ == '__main__':

#    print "Hello world!"

    parser = argparse.ArgumentParser(description = \
                           "Add a new column to a table.", \
                         formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("in_file", help = "File name", type = str)

    parser.add_argument("--position", help = "Position for the new column.", \
                        type = int, default = 0)

    parser.add_argument("--precision", help = "Number of decimal cases.", \
                        type = int, default = 6)

    parser.add_argument("--value", help = "Value to be added.", \
                        type = float, default = 0.0)


    args = parser.parse_args()

    data = np.loadtxt(args.in_file, unpack = False)

    fmt=""
    for i in range(len(data[0]) + 1):
        fmt = fmt + "%." + str(args.precision)  + "E  "
    fmt = fmt[:len(fmt) - 2]


    for i in range(len(data)):
        print fmt % tuple(fc.flatten([data[i][0:args.position], args.value, \
                          data[i][args.position:len(data[i])]]))

