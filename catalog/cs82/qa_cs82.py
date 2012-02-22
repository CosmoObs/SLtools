from __future__ import division
import sys
import os
import errno
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import math
import cmath
from scipy import optimize
import collections
import fit_mag as fmg
import glob
import pyfits
import matplotlib.font_manager
import string

import plot_templates as pltemp

filename = sys.argv[1]

folder_path = os.path.dirname(filename)

# SPATIAL DISTRIBUTION OF VALUES

# reading the data

tiles = np.loadtxt(filename, skiprows=1, dtype='str', usecols=(0,1,11,16,17,18))
tiles = np.char.replace(tiles,'S82','')
tiles = np.char.replace(tiles,'W','')

w4_tiles = tiles[-8:]
cs82_tiles = tiles[0:-8]

print len(cs82_tiles)

# cs82 quadrants

mask_pp = (np.char.count(cs82_tiles[:,0],'p') == 2)
mask_mm = (np.char.count(cs82_tiles[:,0],'p') == 0)
mask_pm = (np.char.count(cs82_tiles[:,0],'p') == 1)*(np.char.find(cs82_tiles[:,0],'p') == 0)
mask_mp = (np.char.count(cs82_tiles[:,0],'p') == 1)*(np.char.find(cs82_tiles[:,0],'p') != 0)

cs82_tiles = np.char.replace(cs82_tiles,'p','')
cs82_tiles = np.char.replace(cs82_tiles,'m','')
cs82_tiles = cs82_tiles.astype(np.float)

cs82_quadrants = []
cs82_quadrants.append(cs82_tiles[mask_mp])
cs82_quadrants.append(cs82_tiles[mask_mm])
cs82_quadrants.append(cs82_tiles[mask_pp])
cs82_quadrants.append(cs82_tiles[mask_pm])

# W4 decs



mask_plus = ((np.char.count(w4_tiles[:,0],'m1') >= 1) & (np.char.rfind(w4_tiles[:,0],'m1') > 2))


w4_tiles = np.char.replace(w4_tiles,'p','')
w4_tiles = np.char.replace(w4_tiles,'m','')
w4_tiles = w4_tiles.astype(np.float)

w4_tiles[:,1] = w4_tiles[:,1] - 360

w4_decs = []
w4_decs.append(w4_tiles[mask_plus])
w4_decs.append(w4_tiles[~mask_plus])


sorted_quadrants = []

for i in range(len(cs82_quadrants)):
    #quadrants[i][:,0] = quadrants[i][:,0]*(-1)*(-1)**(i//2)
    if i < 2:
        cs82_quadrants[i][:,1] = cs82_quadrants[i][:,1] - 360

    sorted_quadrants.append(cs82_quadrants[i][cs82_quadrants[i][:,1].argsort(),:])

cs82_decs=[]

cs82_decs.append(np.vstack((sorted_quadrants[0], sorted_quadrants[2])))
cs82_decs.append(np.vstack((sorted_quadrants[1], sorted_quadrants[3])))


xlabels = ['Comp. Magnitude','Flag = 0 Objects / Total', 'NGal / arcmin^2', 'NStar / arcmin^2']

print range(2,len(cs82_decs[0][0]))

for i in range(2,len(cs82_decs[0][0])):
    fig = plt.figure()
    
    for j in range(2):
        n = int('21'+str(j+1))
        print n
        ax = fig.add_subplot(n)
        ax.bar(cs82_decs[j][:,1], cs82_decs[j][:,i], 0.5)
        ax.bar(w4_decs[j][:,1], w4_decs[j][:,i], 0.5, color = 'r')
        
        if i != 5:
            ax.axhline(np.mean(cs82_decs[j][:,i]), color = 'k', ls = '--', label='CS82 mean = ' + str('%.2f' % np.mean(cs82_decs[j][:,i])))

        ax.set_xlim(-50,50)
        ax.set_ylabel(xlabels[i-2])
        ax.set_xlabel('RA (deg)')
        ax.set_ylim(cs82_decs[1][:,i].min()*0.95,cs82_decs[1][:,i].max()*1.05)
        prop = matplotlib.font_manager.FontProperties(size=12)
        plt.legend(loc='best',prop=prop)

    prop = matplotlib.font_manager.FontProperties(size=12)
    
    plt.savefig(folder_path + '/Plots/' + 'ra_stats'+ str(i-1) + '.png')
    plt.clf()

# HISTOGRAMS NUMBER OF TILES

qa = np.loadtxt(filename, skiprows=1, usecols=(11, 16 ,17, 18))

qa = qa[0:-8]

mn = np.mean(qa,axis=0)
stdev = np.std(qa,axis=0)

bin_sizes = [0.2,np.std(qa[:,1]/4),np.std(qa[:,2]/4),np.std(qa[:,3]/4)]

# the histogram of the data
for i in range(len(qa[0])):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    n, bins, patches = ax.hist(qa[:,i], bins=np.arange(qa[:,i].min(),qa[:,i].max(),bin_sizes[i]),align='left',facecolor='green',label='CS82 mean = ' + str('%.2f' % mn[i]) + '\n CS82 stdev = '+str('%.2f' % stdev[i]))
    print n.max()
    
    prop = matplotlib.font_manager.FontProperties(size=12)
    ax.set_xlabel(xlabels[i])
    ax.set_ylabel('Number of Tiles')
    ax.set_ylim(0,n.max()+2)
    plt.legend(loc='best',prop=prop)
    plt.savefig(folder_path + '/Plots/' + 'histogram_'+ str(i+1) + '.png')
    plt.clf()

# TESTS CORRELATIONS SEEING x NTOT X MAG_COMP

corr = np.loadtxt(filename, skiprows=1, usecols=(3,11,14,8))

corr = corr[0:-8]

# Mag_comp x Seeing

pltemp.plot_data([[corr[:,0],corr[:,1]]], ['.'], ['b'],[''],'Seeing', 'Magnitude', [corr[:,0].max()+0.05,corr[:,0].min()-0.05], [corr[:,1].min()-0.2,corr[:,1].max()+0.2], 'Correlation Completeness Magnitude x Seeing')

plt.savefig(folder_path + '/Plots/' + 'mag_seeing.png')
plt.clf()


# Ngal x Seeing

pltemp.plot_data([[corr[:,0],corr[:,2]]], ['.'], ['b'],[''],'Seeing', 'Ngal', [corr[:,0].max()+0.05,corr[:,0].min()-0.05], [corr[:,2].min()-1000,corr[:,2].max()+1000], 'Correlation Total Number of Galaxies x Seeing')

plt.savefig(folder_path + '/Plots/' + 'ngal_seeing.png')
plt.clf()

# Mag_comp x Ngal

pltemp.plot_data([[corr[:,2],corr[:,1]]], ['.'], ['b'],[''],'Ngal', 'Magnitude', [corr[:,2].min()-1000,corr[:,2].max()+1000], [corr[:,1].min()-0.2,corr[:,1].max()+0.2], 'Correlation Completeness Magnitude x Total Number of Galaxies')

plt.savefig(folder_path + '/Plots/' + 'mag_ngal.png')
plt.clf()

# Mag_comp x b_fit

pltemp.plot_data([[corr[:,3],corr[:,1]]], ['.'], ['b'],[''],'bfit', 'Magnitude', [corr[:,3].min()-0.02,corr[:,3].max()+0.02], [corr[:,1].min()-0.2,corr[:,1].max()+0.2], 'Correlation Completeness Magnitude x b_fit')

plt.savefig(folder_path + '/Plots/' + 'mag_bfit.png')
plt.clf()

# Cortar e imprimir a sublista das tiles que tem magnitude <= 23.0

tiles_23 = np.loadtxt(filename, skiprows=1, dtype='str', usecols=(0,3,11,14))
tiles_23 = tiles_23[(tiles_23[:,2].astype(np.float) <= 23.0)]

tiles_23 = tiles_23[tiles_23[:,1].argsort(),:]
tiles_23 = tiles_23[tiles_23[:,2].argsort(),:]

# print tiles_23[:,4].astype(np.float)

print "| Tile Name | Seeing | Mag_comp | NGal (mag < 99) |"

for i in range(len(tiles_23)):
    print "| " + tiles_23[i][0] + " | " + tiles_23[i][1] + " | " + tiles_23[i][2] + " | " + tiles_23[i][3] + " |"


 

