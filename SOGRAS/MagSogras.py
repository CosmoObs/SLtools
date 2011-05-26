#!/usr/bin/env python
# ====================================================================
# Authors:
# Maria Elidaiana da S. Pereira - mariaeli@cbpf.br
# ====================================================================

import matplotlib
matplotlib.use("Agg")
import sys
import numpy as np
import pylab as pl 

filename = sys.argv[1]

# path where are the cluster's files
cluster=filename.replace('/home/maria/Area_de_Trabalho/SOAR/clusters_sogras/', 'SOGRAS/')
cluster=cluster.replace('/','_')

size_bin=sys.argv[2]

# reading the data
MAG_AUTO_G, MAGERR_AUTO_G, MAG_AUTO_R, MAGERR_AUTO_R, MAG_AUTO_I, MAGERR_AUTO_I, ELLIPTICITY, THETA_IMAGE, CLASS_STAR  = np.loadtxt(filename, usecols=(5, 6, 7, 8, 9, 10, 13, 14 ,15), unpack=True)

# filtering the data
def mag_filter(mag, mag_err, classtar):
    mag_filter=(mag!= 99)*(classtar<0.8)#*(0.186/mag_err>10)    #sigma=0.186/mag
    #print '|', cluster, '|', min(mag[mag_filter]), '|' , max(mag[mag_filter]), '|',  len(mag[mag_filter]), '|',  len(mag), '|'
    return mag_filter

# Function to get the values of histogram : for 1 band
def mag_hist(mag, mag_err,classtar):
    if size_bin == '1':
      binm= [i for i in range(12, 35, 1)]
      histmag = pl.hist(mag[mag_filter(mag, mag_err, classtar)], bins=binm, rwidth=1)
      histm = histmag[0][:]
      binsm = histmag[1][:-1]
    elif size_bin == '0.5':
      binm=[i*0.5 for i in range(26, 51, 1)]
      histmag = pl.hist(mag[mag_filter(mag, mag_err, classtar)], bins=binm, rwidth=0.5)
      histm = histmag[0][:]
      binsm = histmag[1][:-1]
      print histm, binsm
    return [binsm,histm]

# function to plot in log-linear: for 1 band
def mag_plot(magg,magr,magi, classtar):
    if size_bin == '1':
      x=[i for i in range(12, 35, 1)]
    elif size_bin == '0.5':
      x=[i*0.5 for i in range(26, 62, 1)]
    mg=mag_hist(magg, classtar)
    mr=mag_hist(magr, classtar)
    mi=mag_hist(magi, classtar)
    pl.clf()
    pl.figure(figsize=(15,10))
    ax=pl.subplot(111)
    ax.set_yscale("log", nonposy='mask')
    pg=pl.errorbar(mg[0], mg[1], np.sqrt(mg[1]), color='c', ecolor='k')
    pr=pl.errorbar(mr[0], mr[1], np.sqrt(mr[1]), color='b',ecolor='k')
    pi=pl.errorbar(mi[0], mi[1], np.sqrt(mi[1]), color='r',ecolor='k')
    pl.legend((pg[0],pr[0],pi[0]),('g','r','i'), loc = 'upper left')
    pl.xlabel('MAG_AUTO')
    pl.ylabel('N')
    pl.xlim(16,26)
    pl.xticks(x)
    pl.title('Magnitudes: '+cluster)
    pl.savefig('/home/maria/Area_de_Trabalho/SOAR/maglim_resultados_'+'bin'+(size_bin.replace('.',''))+'/'+cluster.replace('.cat','.png'))

# function to plot in log-linear: for 3 bands
def mag_plot_all(magg, magr,magi,magg_err, magr_err,magi_err, classtar):
    if size_bin == '1':
      x=[i for i in range(12, 35, 1)]
    elif size_bin == '0.5':
      x=[i*0.5 for i in range(26, 51, 1)]
    mg=mag_hist(magg, magg_err, classtar)
    mr=mag_hist(magr, magr_err, classtar)
    mi=mag_hist(magi, magi_err, classtar)
    pl.clf()
    pl.figure(figsize=(15,10))
    ax=pl.subplot(111)
    ax.set_yscale("log", nonposy='mask')
    pg=pl.errorbar(mg[0], mg[1], np.sqrt(mg[1]), color='c', ecolor='k')
    pr=pl.errorbar(mr[0], mr[1], np.sqrt(mr[1]), color='b',ecolor='k')
    pi=pl.errorbar(mi[0], mi[1], np.sqrt(mi[1]), color='r',ecolor='k')
    pl.legend((pg[0],pr[0],pi[0]),('g','r','i'), loc = 'upper left')
    pl.xlabel('MAG_AUTO')
    pl.ylabel('N')
    pl.xlim(16,26)
    pl.xticks([i for i in xrange(16, 26, 1)])
    pl.title('Mag Sogras High z')
    pl.savefig('/home/maria/Area_de_Trabalho/SOAR/maglim_resultados_SograsLow_z_bin05/Sogras_Low_z.png')

# call the function
mag_plot_all(MAG_AUTO_G, MAG_AUTO_R, MAG_AUTO_I, MAGERR_AUTO_G, MAGERR_AUTO_R, MAGERR_AUTO_I, CLASS_STAR)
#mag_plot(MAG_AUTO_G, MAG_AUTO_R, MAG_AUTO_I, CLASS_STAR)

