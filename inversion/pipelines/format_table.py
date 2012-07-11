#!/usr/bin/env python
# ==================================
# Author: Gabriel Bartosch Caminha - gbcaminha@gmail.com
# ==================================

"""
FIXME - package description
"""

import numpy as np
import parser as ps

import argparse

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description = \
                           "Re-write a table file using scientific notation.", \
                         formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("in_file", help = "File name", type = str)

    parser.add_argument("--precision", help = "Number of decimal cases.", \
                        type = int, default = 6)

    args = parser.parse_args()

    data = np.loadtxt(args.in_file, unpack = False)

    fmt = ""   
    for i in range(len(data[0])):
        fmt = fmt + "%." + str(args.precision)  + "e  "
    fmt = fmt[:len(fmt) - 2]

    for i in range(len(data)):
        print fmt % tuple(data[i])

