#!/usr/bin/env python

"""Module to plot caustics and critical curves from gravlens"""

##@package plot_gravlens_crit
#
# Executable package: YES
#
# This module uses sltools.geometry.separate_curves to separate radial and tangential curves
#
# To see the package help message, type:
#
# > python plot_gravlens_crit --help
#
# The input of 'plot_gravlens_crit'  is a gravlens file created by the 'plotcrit' command and plotmode 1
#
# The result is a png file with the caustics and critical curves.
#
# Example:
# > python plot_gravlens_crit.py -i testing/curves.txt


from matplotlib.pyplot import *
from numpy import * 
from random import *
import sys

from sltools.geometry import separate_curves


def plot_gravlens_crit(file):
    '''Function to plot gravlens caustics and critical curves
    @param input_file: gravlens gravlens file created by the 'plotcrit' command and plotmode 1
    '''

    # The gravlens plotmode 1 output are consecutive line segments defining a set of 2D curves 

    x1, y1, u1, v1, x2, y2, u2, v2 = loadtxt(file, unpack=True)

    colors='rygbrygbrygb'

    clf() # Clear previous plot
    axes().set_aspect('equal')

    # Plot critical curves

    subplot(121)
    title("Critical curves")

    curves = separate_curves(x1, y1, x2, y2)
    i=0

    for curve in curves:
        x1, y1, x2, y2 = curve
        plot([x1, x2], [y1, y2], color = colors[i])  # plot line segments from (x1,y1) to (x2,y2)
        i=i+1

    # Plot caustics

    subplot(122)
    title("Caustics")

    curves = separate_curves(u1, v1, u2, v2)
    i=0

    for curve in curves:
        u1, v1, u2, v2 = curve
        plot([u1, u2], [v1, v2], color = colors[i])
        i=i+1

    savefig(file + '.png')            

    return 0

# \cond
if __name__ == "__main__":

    from optparse import OptionParser;
    parser = OptionParser();

    parser.add_option('-i','--input_file',
            dest='file', default=None,
            help='Gravlens plotcrit output with plotmode 1');

    (opts,args) = parser.parse_args();
    file = opts.file;

    if ( opts=={} ):
      print"";
      parser.print_help();
      print ""; 
      sys.exit(1);
    if ( file==None ):
      print ""; 
      parser.print_help();
      print ""; 
      sys.exit(1);

    plot_gravlens_crit(file)

    print "Created %s file." % (file + '.png')
# \endcond
