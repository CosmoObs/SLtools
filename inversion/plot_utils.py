# ==================================
# Author: Gabriel Bartosch Caminha - gbcaminha@gmail.com
# ==================================
import numpy as np

import matplotlib.pyplot as plt

from matplotlib.ticker import NullFormatter
from matplotlib import rc

rc('text', usetex=True)
font = {'family' : 'serif', \
        'size' : 17}

rc('font', **font)


################################################################################
#function adapted from:
#http://matplotlib.sourceforge.net/examples/pylab_examples/scatter_hist.html
def plot_scatter_hist(x_in, y_in, xlabel='x', ylabel='y', filename=''):
    # the random data
    #x = np.random.randn(1000)
    #y = np.random.randn(1000)

    nullfmt   = NullFormatter()         # no labels

    # definitions for the axes
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    bottom_h = left_h = left+width+0.02

    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom_h, width, 0.2]
    rect_histy = [left_h, bottom, 0.2, height]

    # start with a rectangular Figure
    plt.figure(1, figsize=(8, 8))

    axScatter = plt.axes(rect_scatter)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    axHistx = plt.axes(rect_histx)
    axHisty = plt.axes(rect_histy)

    # no labels
    axHistx.xaxis.set_major_formatter(nullfmt)
    axHisty.yaxis.set_major_formatter(nullfmt)

    # the scatter plot:
    axScatter.scatter(x_in, y_in, s=2)
    axScatter.set_xlim( (min(x_in), max(x_in)) )
    axScatter.set_ylim( (min(y_in), max(y_in)) )

    #bins = np.arange(-lim, lim + binwidthx, binwidthx)
    axHistx.hist(x_in, bins=20)
    axHisty.hist(y_in, bins=20, orientation='horizontal')

    axHistx.set_xlim( axScatter.get_xlim() )
    axHisty.set_ylim( axScatter.get_ylim() )

    if filename:
        plt.savefig(filename)
    else:
        plt.show()

    #plt.clf()
################################################################################
def plot_histogram(x_in, xlabel='x', bins=20, filename='', legend = ' ', \
                   plot = True, plt_range = None):

    plt.figure(1, figsize=(9, 6))
    plt.hist(x_in, bins, normed=True, label=legend, range = plt_range)
    plt.xlabel(xlabel)

    plt.legend(prop = {'size' : 14} )
    if filename:
        plt.savefig(filename)
        plt.clf()
    elif plot:
        plt.show()


################################################################################
def plot_mult_hist(n_lin, n_col, data, x_axes):
    plt.figure(1, figsize=(n_col*3, n_lin*3))
    sub_num = 1
    
    for lin in range(n_lin):
        for col in range(n_col):
            plt.subplot(n_lin, n_col, sub_num)
            plt.hist(data[lin][col])
            #print n_lin, n_col, sub_num
            #print data[lin][col]
            if lin == n_lin-1: plt.xlabel(x_axes[col])
            sub_num = sub_num + 1
    plt.show()
################################################################################
