from __future__ import division
import sys
import os
import numpy as np
import math
import cmath
from scipy import optimize
import collections


def mag_function(m,a,b):

    if type(m) != 'numpy.ndarray':
        m = np.array(m)
    
    return a*10**(b*m)

#def cut_mag_data(data, m_inf, m_sup):

    

    


def fit_mag_data(data,a,b):

    x = data[0]
    y = data[1]
    yerror = data[2]
    
    p0=[a,b]

    p, cov = optimize.curve_fit(mag_function, x, y, p0, yerror)
    
    return p


X = open(sys.argv[1]).readlines()

m, N, dN = [], [], []

for i in range(12,len(X)-6):
    m.append(float(X[i].split()[0]))
    N.append(float(X[i].split()[1]))
    dN.append(math.sqrt(float(X[i].split()[1])))

data = [m, N, dN]

p = fit_mag_data(data,0,0.3)

print p, type(p)
