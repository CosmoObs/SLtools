from __future__ import division
import sys
import os
import numpy as np
import math
import cmath
from scipy import optimize
import collections



def find_mag_lim(cut_data,fit_params):
    
    a = fit_params[0]
    b = fit_params[1]

    rel_diff = np.absolute(1-cut_data[1]/mag_function(a,b,cut_data[0]))

    print rel_diff
    
    mag_mask = rel_diff > 0.1

    mag_lim = cut_data[0][mag_mask].min()

    return mag_lim


'''# comment this to use the code as module
X = open(sys.argv[1]).readlines()

m, N, dN = [], [], []

for i in range(0,len(X)):
    m.append(float(X[i].split()[0]))
    N.append(float(X[i].split()[1]))
    dN.append(math.sqrt(float(X[i].split()[1])))

data = np.array([m, N, dN])


#fit_params = [0.00000761, 0.22]
fit_params = [0.001864, 0.297]
print find_mag_lim(data,fit_params)
'''
