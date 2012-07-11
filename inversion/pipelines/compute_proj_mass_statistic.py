#!/usr/bin/env python
# ==================================
# Author: Gabriel Bartosch Caminha - gbcaminha@gmail.com
# ==================================

"""
FIXME - package description
"""

import compute_projmass
import lenstool
import functions

import os
import numpy as np
import matplotlib.pylab as plt

class Namespace(object):
    pass

def main_compute_proj_mass_statistic():
    """
    Function to compute the projected mass for several lens inversions, using
    the output from inversion_stat.py
    Input
     -
    Output
     -
    #FIXME - finish documentation
    """
    print "Hello world! from main_compute_proj_mass_statistic()"

    directory = "/home/gcaminha/Programs/gbclib/runs/test/proj_mass/1err0.0_ell0.2_vdisp500.0/"

    args = Namespace()
    args.kappa_map = "kappa_map.fits"
    args.radius = 1.0
    args.plot_cont = False


    mass = []
    mass_im = []

    best_file = []
    im_file = []

    for i in sorted(os.listdir(directory)):
        if "best.par" in i:
            best_file.append(i)
        if "im_generate_images" in i:
            im_file.append(i)

    for i in range(len(im_file)):
        print best_file[i], im_file[i]
        lens = functions.extract_second_identifiers(directory + best_file[i], \
                                                    "potentiel")
        lens_im = functions.extract_second_identifiers(directory + im_file[i], \
                                                    "potential")
        print lens
#        print dir(args)
        lenstool.compute_kappa_map(lens, 50, 1)
        os.system("lenstool kappa_map.par -n")
        mass.append(lenstool.compute_projmass(args))

        lenstool.compute_kappa_map(lens_im, 50, 1)
        os.system("lenstool kappa_map.par -n")
        mass_im.append(lenstool.compute_projmass(args))



    mass = np.array(mass)
    mass_im = np.array(mass_im)

    plt.hist((mass - mass_im)/mass*100, normed = True)
    plt.show()

if __name__ == '__main__':

    print "Hello world!"
    main_compute_proj_mass_statistic()
