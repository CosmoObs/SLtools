#!/usr/bin/env python
# ==================================
# Author: Gabriel Bartosch Caminha - gbcaminha@gmail.com
# ==================================
"""
Package to compute the NFW cross section'
"""

from radial_models import *

from elliptical_models import *

import matplotlib.pyplot as plt

import numpy as np

import constants as cn

def test_ellint():
    ell = 0.45167959731417684
    a = math.sqrt(1.0/(1.0 - ell))
    b = math.sqrt(1.0 - ell)
    tst = {'x1' : 0.3, 'x2' : 0.0, \
           'radial' : 0.223606798, 'theta': 1.1071487177940904, \
           'a' : a, 'b' : b, 'u' : 0.5, \
           'kappa_s' : 0.10157686185488668}

    KG_OBJ1 = KappaGammaEll(NfwLens())
    OBJ1 = LensComputation(KG_OBJ1)
    
    print '\n'
    
    #lens_model2 = EnfwAprox()
    #KG_OBJ2 = KappaGammaEll(EnfwAprox())
    #OBJ2 = LensComputation(KG_OBJ2)

    npts2 = 100
    theta = 0.0
    theta_vec = []
    r_vec = []
    r_vec2 = []

    file_name = 'cc_tang_enfw_ks01_e3.dat'
    f = open(file_name, 'w')
    for i in range(npts2):
        r_vec.append(OBJ1.find_cc_tan(theta, tst['a'], tst['b'], tst['kappa_s'])[0])
        ##r_vec2.append(OBJ2.find_cc_tan(theta, tst['a'], tst['b'], tst['kappa_s'])[0])
        theta_vec.append(theta)
        out_str = str(r_vec[i]*math.cos(theta_vec[i])) + ' ' + str(r_vec[i]*math.sin(theta_vec[i])) + '\n'
        f.write(out_str)
        theta = theta + 2.0*cn.pi/float(npts2-1.0)
    #plt.subplot(1, 1, 1).set_aspect(1)
    #plt.plot(r_vec*np.cos(theta_vec), r_vec*np.sin(theta_vec), '-'  )
    #plt.plot(r_vec2*np.cos(theta_vec), r_vec2*np.sin(theta_vec), '.'  )
    #plt.show()
def func_test(x, a, b, c):
    return a*x + b * c*x**2

if __name__ == '__main__':
    print '---------------------------------------'
    print 'Module to compute the NFW cross section'
    print '---------------------------------------'

    test_ellint()
    #plot_mags()


    
    #test_EnfwAprox()
    #test_PotentialDerivatives()