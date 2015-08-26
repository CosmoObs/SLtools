#!/usr/bin/env python
# ====================================================================
# Authors:
# Maria Elidaiana da S. Pereira - mariaeli@cbpf.br
# ====================================================================

import sys
import glob
import numpy as np
import pylab as pl
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

def scatter_ns(tot_mag, e1, e2, e3, e4, save_path,f):
    k1=(tot_mag>=18)*(tot_mag<=20)  #|(tot_mag2>=18)*(tot_mag2<=20)
    k2=(tot_mag>20)*(tot_mag<=22)   #|(tot_mag2>=20)*(tot_mag2<=22)
    k3=(tot_mag>22)*(tot_mag<=24)   #|(tot_mag2>=22)*(tot_mag2<=24)


    # scatter plot dos ns nao ordenados
    pl.clf()
    pl.figure(figsize=(15,10))
    #pl.plot(e2[k1], e1[k1], 'ko')
    #pl.plot(e2[k2], e1[k2], 'ro')
    #pl.plot(e2[k3], e1[k3], 'yo')
    #pl.legend(('mag: 18 - 20','mag: 20 - 22 ','mag: 22 - 24'),loc='best')
    pl.plot(e2, e1, 'ko', label='Low z')
    pl.plot(e4, e3, 'ro', label='High z')
    #pl.plot(e5, e6, 'yo', label='1 + 2 sersics')
    pl.plot([0, 1,2,23], [0, 1,2,23], 'k-')
    pl.legend(loc='best')   
    pl.title('Scatter Sersic index '+f, fontsize=25)
    pl.xlabel('$ns_A$', fontsize=25)
    pl.ylabel('$ns_B$', fontsize=25)
    pl.xlim(0,2) #1.5
    pl.ylim(0,10)  #5
    pl.savefig(save_path+'Scatter_'+f+'_1e2_sersics_order.png')
    #pl.show()

def hist_disk_bulge(tm1, tm2, tm3, tm4, save_path,f):
    #k=[filter_galclean(tm1)*filter_galclean(tm2)]
    disk_bulge = 10**(((tm2 - tm1)/2.5))
    disk_bulge2 = 10**(((tm4 - tm3)/2.5))
    #disk_bulge3 = 10**(((tm6 - tm5)/2.5))

    # histograma razao disco-bojo
    bindb= [i*0.2 for i in range(-25, 25)]
    ly=max(pl.hist(disk_bulge2, bins=bindb)[0])
    pl.clf()
    pl.figure(figsize=(15,10))
    #pl.hist(disk_bulge, bins=bindb, ec='k', label='disk-bulge', histtype='step',lw=3)# hist(pos_ang, bins=40, width=10, align = 'center')
    pl.hist(disk_bulge, bins=bindb, ec='k', label='Low z', histtype='step',lw=3)
    pl.hist(disk_bulge2, bins=bindb, ec='r', label='High z', histtype='step',lw=3)
    #pl.hist(disk_bulge3, bins=bindb, ec='y', label='1 + 2 sersics', histtype='step',lw=3)
    pl.title('Disk-bulge ratio '+f, fontsize=25)
    pl.xlabel('F_A/F_B',fontsize=25)
    pl.ylabel('N',fontsize=25)
    pl.legend(loc='best')
    pl.xlim(-1,5)
    pl.ylim(0,ly+10)
    pl.savefig(save_path+'Hist_'+f+'_1e2_bulge_to_disk.png')
    ####pl.show()


def hist_disk_boss(db1, db2, tm1, tm2, save_path,f):
    # histograma disk_box
    bind= [i*0.5 for i in range(-6, 11, 1)]
    ly=max(pl.hist(db1, bins=bind)[0])
    #ly=max(max(pl.hist(db1, bins=bind)[0]), max(pl.hist(db2, bins=bind)[0]))
    pl.clf()
    pl.figure(figsize=(15,10))  
    pl.hist(db1, bins=bind, ec='k', histtype='step',lw=3)  
    #pl.hist(db1, bins=bind, ec='k', label='$ns_B$', histtype='step',lw=3)
    #pl.hist(db2, bins=bind, ec='r', label='$ns_A$', histtype='step',lw=3)
    pl.title('Diskyness(-)/ boxyness(+) '+f, fontsize=25)
    pl.legend(loc='best')
    pl.xlabel('c', fontsize=25)
    pl.ylabel('N', fontsize=25)
    pl.xlim(-3,3)
    pl.ylim(0,ly+10)
    pl.savefig(save_path+'Hist_'+f+'_disk_box_2_sersics.png')
    ###pl.show()


def hist_mag_soma(tm1, tm2, tm3, tm4, save_path,f): #(tm1, tm2, save_path,f):
    mtot = np.log10(10**((tm1/2.5)+(tm2/2.5)))
    mtot2 = np.log10(10**((tm3/2.5)+(tm4/2.5)))
    #mtot3 = np.log10(10**((tm5/2.5)+(tm6/2.5)))

    binm= [i*0.5 for i in range(24, 91, 1)]
    ly=max(pl.hist(mtot2, bins=binm)[0])#, max(pl.hist(tm2, bins=binm)[0]))
    pl.clf()
    pl.figure(figsize=(15,10))
    pl.subplot(1,1,1)
    pl.hist(mtot, binm, ec='k', label='Low z', histtype='step',lw=3)
    pl.hist(mtot2, binm, ec='r', label='High z', histtype='step',lw=3)
    #pl.hist(mtot3, binm, ec='y', label='1 + 2 sersics', histtype='step',lw=3)
    #pl.hist(tm1, binm, ec='k', label='$ns_B$', histtype='step',lw=3)
    #pl.hist(tm2, binm, ec='r', label='$ns_A$', histtype='step',lw=3)
    pl.legend(loc='best')
    pl.title('Total magnitude '+f, fontsize=25)
    pl.xlabel('mag', fontsize=25)
    pl.ylabel('N', fontsize=25)
    pl.ylim(0,ly+10)
    pl.savefig(save_path+'Hist_'+f+'_mag_sum_1e2_sersics.png')


def hist_mag(tm1, tm2, save_path,f):
    binm= [i*0.5 for i in range(24, 91, 1)]
    ly=max(max(pl.hist(tm1, bins=binm)[0]), max(pl.hist(tm2, bins=binm)[0]))
    pl.clf()
    pl.figure(figsize=(15,10))
    pl.subplot(1,1,1)
    pl.hist(tm1, binm, ec='k', histtype='step',lw=3)
    #pl.hist(tm1, binm, ec='k', label='$ns_B$', histtype='step',lw=3)
    #pl.hist(tm2, binm, ec='r', label='$ns_A$', histtype='step',lw=3)
    pl.legend(loc='best')
    pl.title('Total magnitude '+f, fontsize=25)
    pl.xlabel('mag', fontsize=25)
    pl.ylabel('N', fontsize=25)
    pl.ylim(0,ly+10)
    pl.savefig(save_path+'Hist_'+f+'_mag_2_sersics.png')
    ##pl.show()


def hist_mag_soma_tot_mag(tm1, tm2, tm3, tm4, save_path,f): #(tm1, tm2, save_path,f):
    mtot = np.log10(10**((tm1/2.5)+(tm2/2.5)))
    #mtot2 = np.log10(10**((tm3/2.5)+(tm4/2.5)))
    #mtot3 = np.log10(10**((tm5/2.5)+(tm6/2.5)))
    #mtot3 = 
    binm= [i*0.5 for i in range(24, 91, 1)]
    ly=max(pl.hist(tm4, bins=binm)[0])#, max(pl.hist(tm2, bins=binm)[0]))
    pl.clf()
    pl.figure(figsize=(15,10))
    pl.subplot(1,1,1)
    #pl.hist(mtot, binm, ec='k', label='2 sersics', histtype='step',lw=3)
    pl.hist(tm1, binm, ec='r', label='1 sersic - col.1', histtype='step',lw=3)
    pl.hist(tm2, binm, ec='b', label='1 sersic - col.2', histtype='step',lw=3)
    pl.hist(tm3, binm, ec='g', label='1 sersic', histtype='step',lw=3)
    pl.hist(tm4, binm, ec='y', label='all sersics', histtype='step',lw=3)
    #pl.hist(tm1, binm, ec='k', label='$ns_B$', histtype='step',lw=3)
    #pl.hist(tm2, binm, ec='r', label='$ns_A$', histtype='step',lw=3)
    pl.legend(loc='best')
    pl.title('Total magnitude '+f, fontsize=25)
    pl.xlabel('mag', fontsize=25)
    pl.ylabel('N', fontsize=25)
    pl.ylim(0,ly+10)
    pl.savefig(save_path+'Hist_'+f+'_mag_sum_all_sersics.png')


def  hist_ax_ratio(ar1, ar2, tm1, tm2, save_path,f):
    # histograma ax_rat 
    binax= [i*0.05 for i in range(1, 21, 1)]
    ly=max(max(pl.hist(ar1, bins=binax)[0]), max(pl.hist(ar2, bins=binax)[0]))
    pl.clf()
    pl.figure(figsize=(15,10))
    pl.hist(ar1, bins=binax, ec='k', histtype='step',lw=3)
    #pl.hist(ar1, bins=binax, ec='k', label='$ns_B$', histtype='step',lw=3)# hist(pos_ang, bins=40, width=10, align = 'center')
    #pl.hist(ar2, bins=binax, ec='r', label='$ns_A$', histtype='step',lw=3)
    pl.legend(loc='best')
    pl.title('Axis ratio '+f, fontsize=25)
    pl.xlabel('b/a', fontsize=25)
    pl.ylabel('N', fontsize=25)
    pl.ylim(0,ly+10)
    pl.savefig(save_path+'Hist_'+f+'_ax_ratio_2_sersics.png')
    ##pl.show()


def hist_ns(e1, e2, tm1, tm2, save_path,f):
    # histograma exp 
    bine= [i*0.1 for i in range(0, 51, 1)]
    ly=max(max(pl.hist(e1, bins=bine)[0]), max(pl.hist(e2, bins=bine)[0]))
    pl.clf()
    pl.figure(figsize=(15,10))
    pl.hist(e1, bins=bine, ec='k', histtype='step', lw=3)
    #pl.hist(e1, bins=bine, ec='k', label='$ns_B$', histtype='step', lw=3)  #hist(pos_ang, bins=40, width=10, align = 'center')
    #pl.hist(e2, bins=bine, ec='r', label='$ns_A$', histtype='step', lw=3)
    pl.legend(loc='best')
    pl.title('Sersic index '+f, fontsize=25)
    pl.xlabel('ns', fontsize=25)
    pl.ylabel('N', fontsize=25)
    pl.xlim(0,5)
    pl.ylim(0,ly+10)
    pl.savefig(save_path+'Hist_'+f+'_ns_2_sersics.png')
    ##pl.show()


def hist_R(r1, r2, tm1, tm2, save_path,f):
    # plot do R_e com fwhm em arsecs
    R_e_arsec = r1*0.154
    R_e_arsec2 = r2*0.154

    kR=(R_e_arsec>=(0.4))
    kR2=(R_e_arsec2>=(0.4))

    binR = [i*0.2 for i in range(22)]
    lmt_y = max(max(pl.hist(R_e_arsec, bins=binR)[0]),max(pl.hist(R_e_arsec2, bins=binR)[0]))

    majorLocator   = MultipleLocator(1)
    majorFormatter = FormatStrFormatter('%d')
    minorLocator   = MultipleLocator(0.2)
    pl.clf()
    pl.figure(figsize=(15,10))
    ax = pl.subplot(111)
    pl.hist(R_e_arsec, bins=binR, ec='k', histtype='step', lw=3)
    #pl.hist(R_e_arsec, bins=binR, ec='k', label='ns_B', histtype='step', lw=3)
    #pl.hist(R_e_arsec2, bins=binR ,ec='r', label='ns_A', histtype='step', lw=3)
    #p1=pl.axvline(x=(0.44), color='k', linestyle='--', linewidth=1)
    pl.title('Effective radius '+f, fontsize=25)
    pl.legend(loc='best')
    ax.xaxis.set_major_locator(majorLocator)
    ax.xaxis.set_minor_locator(minorLocator)
    ax.xaxis.set_minor_locator(minorLocator)
    pl.xlabel('R_e (arcsec)', fontsize=25)
    pl.ylabel('N', fontsize=25)
    pl.ylim(0, lmt_y+10)
    pl.savefig(save_path+'Hist_'+f+'_R_e_2_sersics.png')
    ##pl.show()



def hist_PA(pa1, pa2, tm1, tm2, save_path,f):
    binP= [i for i in xrange(-100, 101, 5)]
    lmt_y=max(pl.hist(pa1, bins=binP)[0])
    #lmt_y = max(max(pl.hist(pa1, bins=binP)[0]),max(pl.hist(pa2, bins=binP)[0]))
    pl.clf()
    pl.figure(figsize=(15,10))
    p1=pl.hist(pa1, bins=binP, ec='k', histtype='step', lw=3)
    #p1=pl.hist(pa1, bins=binP, ec='k', label='$ns_B$', histtype='step', lw=3) #hist(pos_ang, bins=40, width=10, align = 'center')
    #p2=pl.hist(pa2, bins=binP, ec='r', label= '$ns_A$', histtype='step',lw=3)#, color='orange', ec='k')
    #pl.legend(loc='best')
    pl.title('Position angle '+f, fontsize=25)
    pl.xlabel('PA', fontsize=25)
    pl.ylabel('N', fontsize=25)
    pl.xlim(-100,100)
    pl.ylim(0, lmt_y+10)
    pl.savefig(save_path+'Hist_'+f+'_pos_ang_2_sersics.png')
    ##pl.show()   


