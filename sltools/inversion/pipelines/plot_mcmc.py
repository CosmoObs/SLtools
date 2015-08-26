#!/usr/bin/env python
# ==================================
# Author: Gabriel Bartosch Caminha - gbcaminha@gmail.com
# ==================================
from plot_utils import plot_scatter_hist, plot_histogram
import numpy as np


# MAIN
#-------------------------------------------------------------------------------
catalog_path = 'bayes.dat'


mcmc_sample = np.loadtxt(catalog_path, unpack=True)

f = open(catalog_path, 'r')
line_f = f.readline()

comments = []
while line_f[0] == '#':
    comments.append(line_f.split())
    line_f = f.readline()
f.close()

#print comments

params = []
for i in range(2, len(comments)-1):
    params.append(comments[i][2])

#if only one parameter
if len(params) == 1:
    name_out = params[0] + '.eps'
    print 'making ', name_out
    plot_histogram(mcmc_sample[2], xlabel=params[0], filename=name_out)

#compute all parameters combinations
count1 = 0
while count1 < (len(params)-1):
    for i in range(len(params)-1):
        if count1 < i+1:
            name_out = params[count1] + '-' + params[i+1] + '.eps'
            print 'making ', name_out
            plot_scatter_hist(mcmc_sample[count1+2], mcmc_sample[i+3], \
            xlabel = params[count1], ylabel=params[i+1], \
            filename = name_out)
    count1 = count1 + 1

f.close()
#-------------------------------------------------------------------------------
