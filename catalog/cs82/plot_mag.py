#!/usr/bin/env python

import sys
from math import exp
import pylab as pl 
import numpy as np
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
#import myhistogram
import minuit
import matplotlib #import step
import glob
import mag_completeness


filename = sys.argv[1] #'S82p28m_y.V2.3A.swarp.cut_0.cat'

# reading the data
NUMBER, ALPHA_J2000, DELTA_J2000, X_IMAGE, Y_IMAGE, MAG_ISO, MAGERR_ISO, MAG_ISOCOR, MAGERR_ISOCOR, MAG_APER, MAGERR_APER, MAG_AUTO, MAGERR_AUTO, MAG_BEST, MAGERR_BEST, FLUX_RADIUS, CLASS_STAR, MU_MAX, FLAGS = np.loadtxt(filename, usecols=(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14 ,15, 16, 17, 18), unpack=True)

# names for data 
image_name = filename.replace('_y.V2.3A.swarp.cut', '')
image_name = image_name.replace('.cat', '')
print '------------', image_name

m_inf=sys.argv[2]
m_lim=sys.argv[3]
A_value=sys.argv[4]
b_value=sys.argv[5]

'''
#----------------------------
print 'm_inf = ', m_inf
print 'm_sup = ', m_lim
print 'A=', A_value
print 'b=', b_value
'''

#================================= For histograms ==================================================

# function to do counting for histogram
def mag_count(m, binm):
    histmag = pl.hist(m, bins=binm, rwidth=1, visible=False)
    histm = histmag[0][:]
    binsm = histmag[1][:-1]
    idx=histm>10           
    return [binsm[idx], histm[idx]]

def mag_count_cut(m, binm):
    histmag = pl.hist(m, bins=binm, rwidth=1, visible=False)
    histm = histmag[0][:]
    binsm = histmag[1][:-1]
    idx=histm>=0           
    return [binsm[idx], histm[idx]]


# bin for histogram
dn_mag= 0.5
#-----------MUDAR----------
#-- for galaxy
#binm_hist= [i*0.5 for i in range(24, 61, 1)]

#-- for star 
binm_hist= [i*0.5 for i in range(14, 61, 1)]

# functions to plot pdf in histogram
#def f(mag,b): return  (b*((mag**(b-1))/(m_lim**b))*dn_mag*N_tot) # fit para potencia
#def f(mag,b): return  (b*np.log(10.0)*10**((mag-m_lim)*b)*dn_mag*N_tot) # fit para exponencial 
def f(mag,A, b): return  (A*10**(b*mag)) # fit para exponencial com 2 paremtros livres


#-----------MUDAR----------
# filtering the data for histogram
'''
# ---- for galaxy cut:
cut='CLASS_STAR < 0.8'
def mag_filter_hist(mag, classtar, flags):
    mag_filter=(mag!= 99)*(classtar<0.8)                         #*(flags==0)    #*(0.186/mag_err>10)    #sigma=0.186/mag
    return mag_filter

'''
# ---- for star cut
cut='CLASS_STAR > 0.8'
def mag_filter_hist(mag, classtar, flags):
    mag_filter=(mag!= 99)*(classtar>0.8)                         #*(flags==0)    #*(0.186/mag_err>10)    #sigma=0.186/mag
    return mag_filter


x=mag_count_cut(MAG_APER[mag_filter_hist(MAG_APER, CLASS_STAR, FLAGS)], [i*0.5 for i in range(44, 61, 1)])[0]
y=mag_count_cut(MAG_APER[mag_filter_hist(MAG_APER, CLASS_STAR, FLAGS)], [i*0.5 for i in range(44, 61, 1)])[1]
data=np.array([x, y, np.sqrt(y)])

fit_params=[float(A_value), float(b_value)]
mag_comp = mag_completeness.find_mag_lim(data,fit_params)


#----------------------------

# function to plot ML fit curve in histogram
m=np.arange(min(mag_count(MAG_APER[mag_filter_hist(MAG_APER, CLASS_STAR, FLAGS)], binm_hist)[0])-0.5, max(mag_count(MAG_APER[mag_filter_hist(MAG_APER, CLASS_STAR, FLAGS)], binm_hist)[0])+0.5,0.01)

t2=f(m, float(A_value), float(b_value))

# function to plot histogram
def hist_mag(m1, binm_hist, image_name, cut, m_inf, m_lim, A_value, b_value, mag_comp):
    print '------------ mag_comp:', mag_comp  #x and y:'                     #max(mag_count(m1, binm_hist)[0])+0.5
    #a =[int(i) for i in mag_count(m1, binm_hist)[0] if i%2!=1.5 and i%2!=0.5]
    a=[i*0.5 for i in range(26,52,2) ]
    pl.clf()
    pl.figure(figsize=(15,10))
    ax=pl.subplot(1,1,1)
    ax.set_yscale("log", nonposy='mask')
    p1=pl.errorbar(mag_count(m1, binm_hist)[0], mag_count(m1, binm_hist)[1], np.sqrt(mag_count(m1, binm_hist)[1]), color='g', ecolor='k', fmt='-', elinewidth=0.5)
    p3=pl.plot(m, t2, 'r', lw=1, label='A='+str('%.3e' % float(A_value))+', b='+str('%.3f' % float(b_value))+' - ($mag_{min}$'+'= '+str(m_inf)+', $mag_{max}$'+'= '+str(m_lim)+') - '+'$m_{comp}$'+'= '+str(mag_comp))
    #pl.xlim(min(mag_count(m1, binm_hist)[0])-0.5, max(mag_count(m1, binm_hist)[0])+0.5)
    pl.axvline(x=mag_comp, color='k',ls='--')
    #pl.xlim(18.5, 26.0)
    #pl.ylim(10, 100000)

    pl.xlim(16.5, 25)#pl.xlim(12.5, 25)
    pl.ylim(10, 10000)

    pl.legend(loc='best')
    #pl.legend(p3, 'A='+str(A_value)+', b='+str(b_value)+' - ($mag_{min}$'+'= '+str(m_inf)+', $mag_{max}$'+'= '+str(m_lim)+')')
    pl.title(image_name+' MAG_APER ('+cut+')' , fontsize=25)
    pl.xlabel('mag', fontsize=25)
    pl.ylabel('N', fontsize=25)
    #pl.xticks([int(i) for i in mag_count(m1, binm_hist)[0] if i%2!=1.5 and i%2!=0.5])
    #pl.xticks(a)
    pl.savefig(image_name+'_completeness_MAG_APER_star.png')
    #pl.show()

# call the function to plot histogram
hist_mag(MAG_APER[mag_filter_hist(MAG_APER, CLASS_STAR, FLAGS)], binm_hist, image_name, cut, m_inf, m_lim, A_value, b_value, mag_comp)




