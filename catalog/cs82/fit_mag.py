from __future__ import division
import sys
import os
import numpy as np
import matplotlib.pyplot as pl 
import numpy as np
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import math
import cmath
from scipy import optimize
import collections


def get_image_number_tests(input_file):

    string_inf = input_file.index('cut_') + 4
    string_sup = input_file.index('.fit') - 1

    if string_inf == string_sup:
        image_number = float(input_file[input_file.find('.fit')-1])

    else:
        image_number = float(input_file[string_inf:string_sup+1])

    return image_number


def get_image_number_catalogs(input_file):

    string_inf = input_file.index('cut_') + 4
    string_sup = input_file.index('.fit') - 1

    if string_inf == string_sup:
        image_number = float(input_file[input_file.find('.fit')-1])

    else:
        image_number = float(input_file[string_inf:string_sup+1])

    return image_number



def bin_mag_data(mag_data, bin_size):

    bin_inf = math.floor(mag_data.min()/bin_size)*bin_size

    bin_sup = math.ceil(mag_data.max()/bin_size)*bin_size

    bin_vals = np.arange(bin_inf-bin_size/2, bin_sup + bin_size/2, bin_size)

    N, m = np.histogram(mag_data, bins = bin_vals)

    m = np.delete(m,len(m)-1)+bin_size/2

    binned_data = np.array([m, N, np.sqrt(N)])

    return binned_data



def cut_unpopulated_bins(binned_data, N_min):

    mask = binned_data[1] >= N_min

    return binned_data[:,mask]



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



def find_mag_lim(cut_inf_data,fit_params,bin_size):
    
    a = fit_params[0]
    b = fit_params[1]

    rel_diff = 1-cut_inf_data[1]/mag_function(cut_inf_data[0],a,b)
    
    mag_mask = rel_diff > 0.1

    mag_lim = cut_inf_data[0][mag_mask].min()-bin_size

    return mag_lim



def make_plot(binned_data, field_names, folder_path, fit_params, m_inf, m_sup, mag_lim, bin_size, gal_cut, S2N_cut, image_number, plot_inf, plot_sup):

    pl.clf()
    
    pl.figure(figsize=(15,10))

    ax=pl.subplot(1,1,1)
    ax.set_yscale("log", nonposy='mask')

    plot_data = pl.errorbar(binned_data[0], binned_data[1], binned_data[2], color='g', ecolor='k', fmt='-', elinewidth=0.5)

    fit_x = np.arange(m_inf, binned_data[0].max(), 0.01)
    fit_y = mag_function(fit_x, fit_params[0], fit_params[1])

    plot_fit = pl.plot(fit_x, fit_y, 'r', lw=1, label='A='+str('%.3e' % float(fit_params[0]))+', b='+str('%.3f' % float(fit_params[1]))+' - ($mag_{inf}$'+'= '+str(m_inf)+', $mag_{sup}$'+'= '+str(m_sup)+') - '+'$m_{comp}$'+'= '+str(mag_lim)+' , bin='+str(bin_size)+' , S/N>='+str(S2N_cut))
    
    pl.axvline(x=mag_lim, color='k',ls='--')
    
    pl.xlim(plot_inf, plot_sup)
    pl.ylim(8, 18000)

    pl.legend(loc='best')

    pl.xlabel('mag', fontsize=25)
    pl.ylabel('N', fontsize=25)

    if gal_cut:
        pl.title('S82p28m_'+ field_names[0]+ '_gal_' + str(int(image_number)), fontsize=25)
        pl.savefig(folder_path + 'S82p28m_'+ field_names[0]+ '_gal_bins_'+ str(bin_size) + '_S2N_'+ str(int(S2N_cut)) + '_' + str(int(image_number)) + '.png')

    else:
        pl.title('S82p28m_'+ field_names[0]+ '_star_' + str(int(image_number)), fontsize=25)
        pl.savefig(folder_path + 'S82p28m_'+ field_names[0]+ '_star_bins_'+ str(bin_size) + '_S2N_'+ str(int(S2N_cut)) + '_' + str(int(image_number)) + '.png')



