#!/usr/bin/env python

import sys
from math import exp
import pylab as pl 
import numpy as np
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
#import myhistogram
import minuit
import matplotlib #import step
import glob
import os

filenames = glob.glob('*.cat')

filename_param= 'fits_binned_objects_mag_aper_star.dat' #'fits_binned_objects_mag_auto_gal.dat'
ID_test, A_value, b_value, mag_inf, mag_lim = np.loadtxt(filename_param, usecols=(0, 1, 2, 3, 4), unpack=True)


for i in range(len(filenames)):
    os.system('python plot_max_vero_final2.py '+filenames[i]+' '+str(mag_inf[i])+' '+str(mag_lim[i])+' '+str(A_value[i])+' '+str(b_value[i]))
