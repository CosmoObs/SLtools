# ==================================
# Author: Gabriel Bartosch Caminha - gbcaminha@gmail.com
# ==================================
import matplotlib.pyplot as plt
import numpy as np

from math import sqrt

#-------------------------------------------------------------------------------
#especific function to plot the out_put from inversion Function
def plot_hist_inversion():
    ell0, angle0, disp0 = np.loadtxt('conf0.txt', usecols=(0, 1 , 2), unpack=True)
    ell1, angle1, disp1 = np.loadtxt('conf1.txt', usecols=(0, 1 , 2), unpack=True)
    ell2, angle2, disp2 = np.loadtxt('conf2.txt', usecols=(0, 1 , 2), unpack=True)
    ell3, angle3, disp3 = np.loadtxt('conf3.txt', usecols=(0, 1 , 2), unpack=True)

    bins = 15

    plt.figure(1, figsize=(11, 11))

#    plt.subplot(431)
#    plt.hist(ell0, bins, range=[0, 0.3])
#    plt.text(0.1,25, 'v_disp=250', fontsize=15)
#    plt.axis([0,0.31,0,30])

#    plt.subplot(432)
#    plt.hist(angle0, bins, range=[0, 40])
#    plt.axis([0,41, 0, 3.5])

#    plt.subplot(433)
#    plt.hist(disp0, bins, range=[245, 255])
#    plt.axis([245,255,0,12])


########################################
    plt.subplot(434)
    plt.hist(ell1, bins, range=[0.2, 0.3])
    plt.text(0.21, 25, 'v_disp=500', fontsize=15)
    plt.axis([0.2, 0.3, 0, 30])

    plt.subplot(435)
    plt.hist(angle1, bins, range=[17, 23])
    plt.axis([17, 23, 0, 25])

    plt.subplot(436)
    plt.hist(disp1, bins, range=[499, 501])
    plt.axis([499, 501, 0, 20])
########################################
    plt.subplot(437)
    plt.hist(ell2, bins, range=[0.2, 0.3])
    plt.text(0.21,25, 'v_disp=750', fontsize=15)
    plt.axis([0.2, 0.3, 0, 30])

    plt.subplot(438)
    plt.hist(angle2, bins, range=[17, 23])
    plt.axis([17, 23, 0, 25])

    plt.subplot(439)
    plt.hist(disp2, bins, range=[749, 751])
    plt.axis([749, 751, 0, 20])
########################################
    plt.subplot(4, 3, 10)
    plt.hist(ell3, bins, range=[0.2, 0.3])
    plt.text(0.21,25, 'v_disp=1000', fontsize=15)
    plt.axis([0.2,0.3,0,30])
    plt.xlabel('ellipticite')

    plt.subplot(4, 3, 11)
    plt.hist(angle3, bins, range=[17, 23])
    plt.axis([17, 23, 0, 25])
    plt.xlabel('angle_pos')

    plt.subplot(4, 3, 12)
    plt.hist(disp3, bins, range=[999, 1001])
    plt.axis([999, 1001, 0, 20])
    plt.xlabel('v_disp')
########################################

    plt.show()
