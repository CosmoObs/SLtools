#!/usr/bin/env python
# ====================================================================
# Authors:
# Maria Elidaiana da S. Pereira - mariaeli@cbpf.br
# ====================================================================

import matplotlib
matplotlib.use("Agg")
import sys
import numpy as np
import matplotlib
import pylab as pl 


f =  'Low_and_High_z' 
#filename = sys.argv[1]
#f= filename.replace('.cat', '')
#size_bin=sys.argv[2]

'''
MAG_AUTO_G, MAGERR_AUTO_G, MAG_AUTO_R, MAGERR_AUTO_R, MAG_AUTO_I, MAGERR_AUTO_I, ELLIPTICITY, THETA_IMAGE, CLASS_STAR, dr8_phot_z, dr8_phot_z_err  = np.loadtxt(filename, usecols=(5, 6, 7, 8, 9, 10, 13, 14 ,15, 19, 20), unpack=True)

'''

MAG_AUTO_G, MAGERR_AUTO_G, MAG_AUTO_R, MAGERR_AUTO_R, MAG_AUTO_I, MAGERR_AUTO_I, ELLIPTICITY, THETA_IMAGE, CLASS_STAR, dr8_phot_z, dr8_phot_z_err  = np.loadtxt('Sogras_GalcleanPhotoz_Low_z.cat', usecols=(5, 6, 7, 8, 9, 10, 13, 14 ,15, 19, 20), unpack=True)

MAG_AUTO_GH, MAGERR_AUTO_GH, MAG_AUTO_RH, MAGERR_AUTO_RH, MAG_AUTO_IH, MAGERR_AUTO_IH, ELLIPTICITYH, THETA_IMAGEH, CLASS_STARH, dr8_phot_zH, dr8_phot_z_errH  = np.loadtxt('Sogras_GalcleanPhotoz_High_z.cat', usecols=(5, 6, 7, 8, 9, 10, 13, 14 ,15, 19, 20), unpack=True)

# function for filtering
def mag_filter(mag, mag_err, classtar):
    mag_filter=(mag!= 99)*(classtar<0.8)*(0.186/mag_err>10)    #sigma=0.186/mag
    #print '|', cluster, '|', min(mag[mag_filter]), '|' , max(mag[mag_filter]), '|',  len(mag[mag_filter]), '|',  len(mag), '|'
    return mag_filter

def photoz_filter(photoz, photoz_err, classtar, mag, mag_err):
    photoz_filter=(photoz>0)*(photoz<10)*(classtar<0.8)       #*(0.186/mag_err>10)    #sigma=0.186/mag
    #print '|', cluster, '|', min(mag[mag_filter]), '|' , max(mag[mag_filter]), '|',  len(mag[mag_filter]), '|',  len(mag), '|'
    return photoz_filter

# fuction to plot histogram of redshift: individual
def photoz_plot(pz, e_pz, classtar, mag, mag_err, f):
    binz= [i*0.1 for i in range(0, 101, 1)]
    ly=max(pl.hist(pz[photoz_filter(pz, e_pz, classtar, mag, mag_err)], bins=binz)[0]) #, max(pl.hist(tm2, bins=binm)[0]))
    #
    pl.clf()
    pl.figure(figsize=(15,10))
    pl.subplot(1,1,1)
    pl.hist(pz[photoz_filter(pz, e_pz, classtar, mag, mag_err)], binz, ec='k', histtype='step',lw=3) # label='ns_1',
    #pl.legend(loc='best')
    pl.title('Sogras Redshift '+f, fontsize=25)
    pl.xlabel('z', fontsize=25)
    pl.ylabel('N', fontsize=25)
    pl.ylim(0,ly+10)
    pl.xlim(0,1)
    pl.xticks([i*0.2 for i in xrange(0, 5, 1)])
    #pl.savefig(save_path+'Hist_'+f+'_photoz.png')
    pl.savefig('/home/maria/Area_de_Trabalho/SOAR/Photoz_Sogras/'+f+'.png')

# fuction to plot histogram of redshift: low and high
def photoz_plot_all(pz, e_pz, classtar, mag, mag_err, pzH, e_pzH, classtarH, magH, mag_errH, f):
    binz= [i*0.1 for i in range(0, 101, 1)]
    ly=max(pl.hist(pzH[photoz_filter(pzH, e_pzH, classtarH, magH, mag_errH)], bins=binz)[0]) #, max(pl.hist(tm2, bins=binm)[0]))
    #
    pl.clf()
    pl.figure(figsize=(15,10))
    pl.subplot(1,1,1)
    pl.hist(pz[photoz_filter(pz, e_pz, classtar, mag, mag_err)], binz, ec='k', label='Low z', histtype='step',lw=3, normed=1) # normed is n/(len(x)*dbin)
    pl.hist(pzH[photoz_filter(pzH, e_pzH, classtarH, magH, mag_errH)], binz, ec='r', label='High z', histtype='step',lw=3, normed=1)
    pl.legend(loc='best')
    pl.title('Sogras Redshift '+f, fontsize=25)
    pl.xlabel('z', fontsize=25)
    pl.ylabel('N', fontsize=25)
    #pl.ylim(0,ly+100)
    #pl.ylim(0,1)
    pl.xlim(0,1)
    pl.xticks([i*0.2 for i in xrange(0, 5, 1)])
    pl.show()
    pl.savefig('Hist_'+f+'_photoz_Normed.png')
    

# call the fuction
photoz_plot_all(dr8_phot_z, dr8_phot_z_err, CLASS_STAR, MAG_AUTO_I, MAGERR_AUTO_I, dr8_phot_zH, dr8_phot_z_errH, CLASS_STARH, MAG_AUTO_IH, MAGERR_AUTO_IH,f) 
#photoz_plot(dr8_phot_z, dr8_phot_z_err, CLASS_STAR, MAG_AUTO_I, MAGERR_AUTO_I, f)  




