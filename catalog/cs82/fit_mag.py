from __future__ import division
import sys
import os
import numpy as np
import math
import cmath
from scipy import optimize
import collections


def bin_mag_data(mag_data, bin_size):

    bin_inf = math.floor(mag_data.min()/bin_size)*bin_size

    bin_sup = math.ceil(mag_data.max()/bin_size)*bin_size

    bin_vals = np.arange(bin_inf-bin_size/2, bin_sup + bin_size/2, bin_size)

    N, m = np.histogram(mag_data, bins = bin_vals)

    m = np.delete(m,len(m)-1)+bin_size/2

    binned_data = np.array([m, N, np.sqrt(N)])

    return binned_data



def mag_function(m,a,b):

    if type(m) != 'numpy.ndarray':
        m = np.array(m)
    
    return a*10**(b*m)


def find_mag_sup(data):

    if type(data) != 'numpy.ndarray':
        data = np.array(data)

    mag_sup = data[0][data[1].argmax()] - 1

    return mag_sup


def cut_mag_data(data, m_inf, m_sup):
    
    if type(data) != 'numpy.ndarray':
        data = np.array(data)

    mag_inf_mask = data[0] >= m_inf
    mag_sup_mask = data[0] <= m_sup

    cut_data = data[:, mag_sup_mask & mag_inf_mask]

    return cut_data


def fit_mag_data(data,a,b):

    if type(data) != 'numpy.ndarray':
        data = np.array(data)

    x = data[0]
    y = data[1]
    yerror = data[2]
    
    p0=[a,b]

    p, cov = optimize.curve_fit(mag_function, x, y, p0, yerror)

    return p


def find_mag_lim(cut_inf_data,fit_params):
    
    a = fit_params[0]
    b = fit_params[1]

    rel_diff = np.absolute(1-cut_inf_data[1]/mag_function(cut_inf_data[0],a,b))
    
    mag_mask = rel_diff > 0.1

    mag_lim = cut_inf_data[0][mag_mask].min()

    return mag_lim


############################################## Pseudo Main Function

'''
X = open(sys.argv[1]).readlines()

# Use later in the main function: data = np.transpose(np.loadtxt("/home/brunomor/Documents/Trabalho/CS82/CS82_SExtractor/Maria_Python_Code_Temporary/test_data.txt"))


m, N, dN = [], [], []

for i in range(0,len(X)):
    m.append(float(X[i].split()[0]))
    N.append(float(X[i].split()[1]))
    dN.append(math.sqrt(float(X[i].split()[1])))

data = [m, N, dN]

data = np.array(data)


mag_sup = find_mag_sup(data)

cut_data = cut_mag_data(data,20,mag_sup)

cut_inf_data = cut_mag_data(data,20,99)

p = fit_mag_data(cut_data,0,0.3)

print p

# Use later in the main function: np.savetxt("savetxt_test",np.transpose(cut_data),fmt='%f')

fit_params = list(p)

print find_mag_lim(cut_inf_data,fit_params)
'''
