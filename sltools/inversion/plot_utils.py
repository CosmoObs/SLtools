# ==================================
# Author: Gabriel Bartosch Caminha - gbcaminha@gmail.com
# ==================================
import numpy as np

import functions as fc

import constants as cs

#change the matplotlib display to allow to run in a terminal with no X window
#FIXME - not tested yet
if not fc.is_X_running():
    import matplotlib
    matplotlib.use('agg')

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
def plot_contours(x_vec, y_vec, z_vec, linestyles = '-', plt_labels = False, \
                  manual = False, levels = None):
    """
    Function to plot contours givem a grid data

    Input:
     - x_vec         [] : x vector
     - y_vec         [] : y vector
     - z_vec         [] : z vector
     - linestyles   str : matplotlib linestyle
     - plt_labels  bool : plot or not the contours values
     - manual      bool : set (or not) the labels manualy
     - levels        [] : values for the iso-contours levels
    Output:
     - NONE
    """

    plt.figure(1, figsize=(10, 10))

    x0_dim = len( ( y_vec == min(y_vec) ).nonzero()[0] )
    x1_dim = len( ( x_vec == min(x_vec) ).nonzero()[0] )

    print "plot_contours: Axis 1 points = ", x0_dim
    print "plot_contours: Axis 1 min and max = ", min(x_vec), " ", max(x_vec)

    print "plot_contours: Axis 2 points = ", x1_dim
    print "plot_contours: Axis min and max = ", min(y_vec), " ", max(y_vec)

    x0_vec = np.linspace(min(x_vec), max(x_vec), x0_dim)
    x1_vec = np.linspace(min(y_vec), max(y_vec), x1_dim)

    x2_vec =  np.reshape(z_vec, (x0_dim, x1_dim)).transpose()

    CS = plt.contour ( x0_vec, x1_vec, x2_vec, levels = levels, colors = "k", \
                       linewidths = 2.5, linestyles = linestyles)

    if plt_labels:
        plt.clabel(CS, fontsize=13, fmt = "%.2f", manual = manual, \
                   rightside_up = True)
################################################################################
def plot_contour_func(func, x_min, x_max, x_npt, y_min, y_max, y_npt, *args):
    """
    Function to add plot contours for a function.
    The contours will be ploted considering the two first arguments of 'func'

    Input:
     - func  function : function used to plot the contours
     - x_min    float : x minimum value
     - x_max    float : x maximum value
     - x_npt    float : grid size along the x axies
     - y_min    float : y minimum value
     - y_max    float : y maximum value
     - y_npt    float : grid size along the y axies
     - *args    *args : 'func' arguments
    Output:
     - NONE
    """
    vec_a = np.linspace(x_min, x_max, num = x_npt)
    vec_b = np.linspace(y_min, y_max, num = y_npt)

    r_vec = np.zeros([len(vec_a), len(vec_b)])

    x = []
    y = []
    for i in range(len(vec_a)):
        for j in range(len(vec_b)):
            r_vec[i][j] = func(vec_a[i], vec_b[j], *args)
            x.append( vec_a[i] )
            y.append( vec_b[j] )

#    levels = [2, 3, 5, 10, 50]
    levels = [0.01, 0.03, 0.05, 0.07]
    plot_contours(x, y, r_vec, levels = levels, plt_labels = True)
################################################################################
def plot_circ(radius):
    """
    Function to plot a simple circumference

    Input:
     - radius  float : circumference radius
    Output:
     - NONE
    """
    theta = np.linspace(0.0, 2.0*cs.pi, num = 200)
    x = radius*np.cos(theta)
    y = radius*np.sin(theta)
    plt.plot(x, y, linewidth = 2)
################################################################################
def plot_func(func, x_min, x_max, *args):
    """
    Function to plot a function

    Input:
     - func  function : functions to be ploted
     - *args    *args : 'func' arguments
    Output:
     - NONE
    """

    n_pts = 200
    #x_min = 1.0E6
    #x_max = 1.0E11

    loglog = False

    x_vec = np.linspace(x_min, x_max, num = n_pts)
    y_vec = []
    for i in x_vec:
        y_vec.append( func(i, *args) )

    if loglog:
        plt.loglog(x_vec, y_vec)
    else:
        plt.plot( x_vec, y_vec )
