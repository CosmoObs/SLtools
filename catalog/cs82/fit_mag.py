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


'''

# This is not needed, except if we can be sure that astrometry is always performed
# with swarp. If this is the case, then we can cut out all the suffixes.

def get_tile_name(input_file):

    input_name = os.path.basename(input_file)

    string_sup = input_name.index('.V2.') - 2

    tile_name = input_name[0:string_sup]

    return tile_name
'''

def get_entry_tile_cs82_ldac(hdulist_string, entry_name):

    ra = hdulist_string[hdulist_string.find(entry_name):]

    ra = ra[0:ra.find('/')]

    ra = ''.join(ra.split())

    ra = str(float(ra[ra.find('=')+1:]))

    return ra




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



def make_plot_mag_pipe(binned_data, tile_name, folder_path, fit_params, m_inf, m_sup, mag_lim, bin_size, gal_cut, S2N_cut, stell_idx, plot_inf, plot_sup):

    pl.clf()
    
    pl.figure(figsize=(15,10))

    ax=pl.subplot(1,1,1)
    ax.set_yscale("log", nonposy='mask')

    plot_data = pl.errorbar(binned_data[0], binned_data[1], binned_data[2], color='g', ecolor='k', fmt='-', elinewidth=0.5)

    fit_x = np.arange(m_inf, binned_data[0].max(), 0.01)
    fit_y = mag_function(fit_x, fit_params[0], fit_params[1])

    plot_fit = pl.plot(fit_x, fit_y, 'r', lw=1, label='A='+str('%.3e' % float(fit_params[0]))+', b='+str('%.3f' % float(fit_params[1]))+' - ($mag_{inf}$'+'= '+str(m_inf)+', $mag_{sup}$'+'= '+str(m_sup)+') - '+'$m_{comp}$'+'= '+str(mag_lim)+' , bin='+str(bin_size)+' , S/N>='+str(S2N_cut)+' , CLASS_STAR >(<)'+str(stell_idx))
    
    pl.axvline(x=mag_lim, color='k',ls='--')
    
    pl.xlim(plot_inf, plot_sup)
    pl.ylim(8, 18000)

    pl.legend(loc='best')

    pl.xlabel('mag', fontsize=25)
    pl.ylabel('N', fontsize=25)

    if gal_cut:
        pl.title(tile_name + '_gal', fontsize=25)
        pl.savefig(folder_path + '/Plots/'+ tile_name + '_gal_bins_'+ str(bin_size) + '_S2N_'+ str(int(S2N_cut)) + '_stell_'+ str(stell_idx) + '.png')

    else:
        pl.title(tile_name + '_star', fontsize=25)
        pl.savefig(folder_path + '/Plots/'+ tile_name + '_star_bins_'+ str(bin_size) + '_S2N_'+ str(int(S2N_cut)) + '_stell_'+ str(stell_idx) + '.png')



