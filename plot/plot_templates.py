from __future__ import division
from matplotlib import pyplot as plt
import numpy as np
import math
import cmath
from time import strftime


def plot_data(multiple_data_arrays, formats, colors, legends, x_label, y_label, x_lims, y_lims, plot_title):

    for i in range(len(multiple_data_arrays)):
        
        plt.plot(multiple_data_arrays[i][0], multiple_data_arrays[i][1], formats[i], color = colors[i], label = legends[i])

    plt.legend(loc = 'best',numpoints=1)

    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(plot_title + ' - ' + strftime("%d %b %Y"))
    plt.xlim(x_lims[0], x_lims[1])
    plt.ylim(y_lims[0], y_lims[1])

    return

data = []

x = np.arange(10)
data.append([x,x])
data.append([x,x**2])
data.append([x,np.sqrt(x)])
data = np.array(data)


plot_data(data, ['.','.','.'], ['b','g','r'], ['mag1','mag2','mag3'], 'mag', 'counts', [x.min(), x.max()], [data.min(), data.max()], 'Test Plot')

plt.savefig('test.png')
plt.clf()
