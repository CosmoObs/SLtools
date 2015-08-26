from __future__ import division
import numpy as np
import math
import cmath
from time import strftime
import matplotlib.font_manager as fm
from matplotlib import pyplot as plt

def plot_data(multiple_data_array_1, multiple_data_array_2, formats, colors, legends, legends_loc, x_label, y_label, x_lims, y_lims, plot_title):

    for i in range(len(multiple_data_array_1)):
        
        plt.plot(multiple_data_array_1[i], multiple_data_array_2[i], formats[i], color = colors[i], label = legends[i])

    prop = fm.FontProperties(size=10)
    plt.legend(loc = legends_loc, numpoints=1, prop=prop)
    # bbox_to_anchor = (0., 1.) could be useful to chose the place of the legend

    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(plot_title + ' - ' + strftime("%d %b %Y"), fontsize=14)
    plt.xlim(x_lims[0], x_lims[1])
    plt.ylim(y_lims[0], y_lims[1])

    return

