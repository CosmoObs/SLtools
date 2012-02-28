#!/usr/bin/env python
# ==================================
# Author: Gabriel Bartosch Caminha - gbcaminha@gmail.com
# ==================================
"""
Package to compute the NFW cross section'
"""

import math

import radial_models as rm

import elliptical_models as em

import cross_section_sep as cs

import matplotlib.pyplot as plt

import numpy as np

import constants as cn

from matplotlib.ticker import NullFormatter
from matplotlib import rc

rc('text', usetex=True)
font = {'family' : 'serif', \
        'size' : 17}

rc('font', **font)

def test_PotentialDerivatives():
    tst = {'x1' : 0.3, 'x2' : 0.2, \
            'radial' : 0.36055512754639896, 'theta': 0.5880026035475676, \
            'u' : 0.5, 'a' : 1.414213562, 'b' : 0.707106781, 'u' : 0.5, \
            'kappa_s' : 0.05}
    ENFW = rm.EnfwAprox()
    PD = em.LensPotDrivatives(ENFW)
    #print dir(PD)
    print cs.pot_xx(tst['radial'], tst['theta'], tst['a'], tst['b'], tst['kappa_s'])
    print cs.pot_xy(tst['theta'], tst['a'], tst['b'], tst['kappa_s'])
    print  tst['a']*tst['b']*tst['x1']*\
           cs.j0_function(tst['radial'], tst['theta'], tst['a'], tst['b'], tst['kappa_s'])

    
    print 'pot_x = ', PD.pot_x(tst['x1'], tst['x2'], tst['a'], tst['b'], tst['kappa_s'])
    print 'pot_y = ', PD.pot_y(tst['x1'], tst['x2'], tst['a'], tst['b'], tst['kappa_s'])
    print 'pot_xx = ', PD.pot_xx(tst['x1'], tst['x2'], tst['a'], tst['b'], tst['kappa_s'])
    print 'pot_yy = ', PD.pot_yy(tst['x1'], tst['x2'], tst['a'], tst['b'], tst['kappa_s'])
    print 'pot_xy = ', PD.pot_xy(tst['x1'], tst['x2'], tst['a'], tst['b'], tst['kappa_s'])

    KG = em.KappaGammaEll(ENFW)

    print 'kappa_num = ',  KG.kappa_numeric(tst['x1'], tst['x2'], tst['a'], tst['b'], tst['kappa_s'])
    print 'kappa_ell = ',  KG.kappa_ell(tst['x1'], tst['x2'], tst['a'], tst['b'], tst['kappa_s'])
    print 'gamma1 = ',  KG.gamma1(tst['x1'], tst['x2'], tst['a'], tst['b'], tst['kappa_s'])
    print 'gamma2 = ',  KG.gamma2(tst['x1'], tst['x2'], tst['a'], tst['b'], tst['kappa_s'])
    print 'gamma = ',  KG.gamma(tst['x1'], tst['x2'], tst['a'], tst['b'], tst['kappa_s'])
    print 'mag_tot = ',  KG.mag_total(tst['x1'], tst['x2'], tst['a'], tst['b'], tst['kappa_s'])
    print 'mag_all = ', KG.mag_all(tst['x1'], tst['x2'], tst['a'], tst['b'], tst['kappa_s'])
    print 'mag_inv = ', KG.mag_inv_all(tst['x1'], tst['x2'], tst['a'], tst['b'], tst['kappa_s'])
################################################################################

def plot_critical_lines():
    ell = 0.45167959731417684
    a = math.sqrt(1.0/(1.0 - ell))
    b = math.sqrt(1.0 - ell)
    tst = {'x1' : 0.3, 'x2' : 0.2, \
            'radial' : 0.36055512754639896, 'theta': 0.5880026035475676, \
            'u' : 0.5, 'a' : 1.414213562, 'b' : 0.707106781, 'u' : 0.5, \
            'kappa_s' : 0.1}

    lens_apro_num = em.KappaGammaEll( rm.EnfwAprox() )
    lens_exact =    em.KappaGammaEll( rm.NfwLens()   )
    
    comp_aprox_num  = em.LensComputation(lens_apro_num)
    comp_lens_exact = em.LensComputation(lens_exact)

    npts = 25
    thetha_vec = np.linspace(0.0, cn.pi/2.0, npts)

    aprox_num_vec = []
    exact = []
    aprox_ana_vec = []
    for i in thetha_vec:
        #print i
        aprox_num_vec.append(comp_aprox_num.find_cc_tan(i, tst['a'], tst['b'], \
                             tst['kappa_s'])[0] )

        exact.append(comp_lens_exact.find_cc_tan(i, tst['a'], tst['b'], \
                     tst['kappa_s'])[0])

        aprox_ana_vec.append(cs.critica_curve_tangential(i, tst['a'], tst['b'],\
                             tst['kappa_s']))



    plt.plot(exact*np.cos(thetha_vec), exact*np.sin(thetha_vec), \
             label='ENFW Numeric', linewidth = 3)

    plt.plot(aprox_num_vec*np.cos(thetha_vec), \
             aprox_num_vec*np.sin(thetha_vec), 'm--', label='Aprox Numeric', \
             linewidth = 3)

    plt.plot(aprox_ana_vec*np.cos(thetha_vec), \
             aprox_ana_vec*np.sin(thetha_vec), 'o', label='Aprox Analytic')

    plt.legend( loc='botton left' )
    plt.suptitle(r'\Tex $\kappa_s='+str(tst['kappa_s']) + '$', y = 0.94)
    plt.xlabel('$x_1$')
    plt.ylabel('$x_2$')
    plt.show()
################################################################################
def plot_constant_distortion():

    ell = 0.2
    a = math.sqrt(1.0/(1.0 - ell))
    b = math.sqrt(1.0 - ell)

    tst = {'x1' : 0.3, 'x2' : 0.2, \
            'radial' : 0.36055512754639896, 'theta': 0.5880026035475676, \
            'u' : 0.5, 'a' : a, 'b' : b, 'u' : 0.5, \
            'kappa_s' : 0.15}

    lens_apro_num = em.KappaGammaEll( rm.EnfwAprox() )
    lens_exact =    em.KappaGammaEll( rm.NfwLens()   )
    
    comp_aprox_num  = em.LensComputation(lens_apro_num)
    comp_lens_exact = em.LensComputation(lens_exact)

    npts = 25
    thetha_vec = np.linspace(0.0, cn.pi/2.0, npts)

    raz = 1.0/10.0
    
    cc_tan = []
    cc_rad = []
    dist_p = []
    dist_m = []

    cc_tan2 = []
    cc_rad2 = []
    dist_p2 = []
    dist_m2 = []

    dist_p3 = []
    dist_m3 = []
    cc_tan3 = []
    for i in thetha_vec:
        curves =  comp_aprox_num.find_constant_distortion(i, tst['a'], tst['b'], \
                  tst['kappa_s'], raz)
        cc_tan.append(curves[0])
        cc_rad.append(curves[1])
        dist_m.append(curves[2])
        dist_p.append(curves[3])

        curves2 =  comp_lens_exact.find_constant_distortion(i, tst['a'], tst['b'], \
              tst['kappa_s'], raz)
        cc_tan2.append(curves2[0])
        cc_rad2.append(curves2[1])
        dist_m2.append(curves2[2])
        dist_p2.append(curves2[3])


        dist_p3.append( cs.const_dist_curve(i, tst['a'], tst['b'], \
                        tst['kappa_s'], 1.0/raz)[1] )
        dist_m3.append( cs.const_dist_curve(i, tst['a'], tst['b'], \
                        tst['kappa_s'], -1.0/raz)[1] )
        cc_tan3.append( cs.critica_curve_tangential(i, tst['a'], tst['b'],\
                             tst['kappa_s']) )


    plt.plot( dist_p2*np.cos(thetha_vec), dist_p2*np.sin(thetha_vec), \
              'b-', linewidth = 3)
    plt.plot( dist_m2*np.cos(thetha_vec), dist_m2*np.sin(thetha_vec), \
              'b-', linewidth = 3)
    plt.plot( cc_tan2*np.cos(thetha_vec), cc_tan2*np.sin(thetha_vec), \
              'b-', label='ENFW Numeric', linewidth = 3)


    plt.plot( dist_m*np.cos(thetha_vec), dist_m*np.sin(thetha_vec), \
              'm--', linewidth = 3)
    plt.plot( dist_p*np.cos(thetha_vec), dist_p*np.sin(thetha_vec), \
              'm--', linewidth = 3)
    plt.plot( cc_tan*np.cos(thetha_vec), cc_tan*np.sin(thetha_vec), \
              'm--', label='Aprox Numeric', linewidth = 3)


    plt.plot( dist_m3*np.cos(thetha_vec), dist_m3*np.sin(thetha_vec), \
              'go', linewidth = 3)
    plt.plot( cc_tan3*np.cos(thetha_vec), cc_tan3*np.sin(thetha_vec), \
              'go', label='Aprox Analytic', linewidth = 3)
    plt.plot( dist_p3*np.cos(thetha_vec), dist_p3*np.sin(thetha_vec), \
              'go', linewidth = 3)


    plt.legend(loc='lower left')#'center' 'lower left'
    plt.suptitle(r'\Tex $\kappa_s='+str(tst['kappa_s']) + '$', y = 0.94)
    plt.xlabel('$x_1$')
    plt.ylabel('$x_2$')
    plt.show()
################################################################################
def cross_section_computation():
    ell = 0.0
    a = math.sqrt(1.0/(1.0 - ell))
    b = math.sqrt(1.0 - ell)

    tst = {'x1' : 0.3, 'x2' : 0.2, \
            'radial' : 0.36055512754639896, 'theta': 0.5880026035475676, \
            'u' : 0.5, 'a' : a, 'b' : b, 'u' : 0.5, \
            'kappa_s' : 0.05}

    lens_apro_num = em.KappaGammaEll( rm.EnfwAprox() )
    lens_exact =    em.KappaGammaEll( rm.NfwLens()   )
    
    comp_aprox_num  = em.LensComputation(lens_apro_num)
    comp_lens_exact = em.LensComputation(lens_exact)

    r_lambda = 10.0
    npts = 1001
    kappa_s_vec = thetha_vec = np.linspace(0.05, 0.25, npts)

    
    file_out = open('cross_section.txt', 'w')
    for ks in kappa_s_vec:
        out_str = str('%.10e' % ks) + ' ' + \
                  str('%.10e' % cs.sigma(a, b, ks, r_lambda)) + '\n'
        print out_str
        file_out.write( out_str )

    file_out.close()
################################################################################
def compare_habib():
    file_habib = 'sigma_pnfw_c-e0.dat'
    file_gbclb = 'cross_section.txt'
    habib_values = np.loadtxt(file_habib, unpack=True)
    gbclb_values = np.loadtxt(file_gbclb, unpack=True)

    diff = gbclb_values[1]-habib_values[1]
    diff_rel = diff/gbclb_values[1]




    plt.semilogy(gbclb_values[0], gbclb_values[1], 'g-', linewidth = 3, \
                 label = 'ENFW')
    plt.semilogy(habib_values[0], habib_values[1], 'b--', linewidth = 3, \
                 label = 'PNFW')

    plt.ylabel(r'$\tilde\sigma$')
    plt.suptitle(r'$\varepsilon = 0$', y = 0.94)
    plt.subplots_adjust(hspace=0)
    plt.legend(loc='lower right')#'center' 'lower left'
    plt.show()

    plt.subplot(2, 1, 1)
    plt.plot(habib_values[0], diff)
    plt.ylabel(r'$\Delta\tilde\sigma$')
    plt.xlabel(r'$\kappa_s$')
    plt.subplot(2, 1, 2)
    plt.plot(habib_values[0], diff_rel)
    plt.ylabel(r'$\Delta\tilde\sigma/\tilde\sigma$')
    plt.xlabel(r'$\kappa_s$')

    plt.show()
################################################################################
def plot_magnifications_thetafix():
    ell = 0.45167959731417684
    a = math.sqrt(1.0/(1.0 - ell))
    b = math.sqrt(1.0 - ell)
    tst = {'x1' : 0.3, 'x2' : 0.2, \
            'radial' : 0.36055512754639896, 'theta': 0.5880026035475676, \
            'u' : 0.5, 'a' : 1.414213562, 'b' : 0.707106781, 'u' : 0.5, \
            'kappa_s' : 1.0}

    lens_apro_num = em.KappaGammaEll( rm.EnfwAprox() )
    lens_exact =    em.KappaGammaEll( rm.NfwLens()   )
    
    comp_aprox_num  = em.LensComputation(lens_apro_num)
    comp_lens_exact = em.LensComputation(lens_exact)
    
    npts = 200
    radial_vec = np.linspace(1E-5, 1.5, npts)
    mag_ratio = []
    arg_find = []

    mag_radial = []
    mag_tan = []

    for i in radial_vec:
        mag_radial.append( 1.0/comp_aprox_num.mag_rad_inv_polar(i, 0.0, tst['a'], \
                         tst['b'], tst['kappa_s']) )
        mag_tan.append( 1.0/comp_aprox_num.mag_tan_inv_polar(i, 0.0, tst['a'], \
                         tst['b'], tst['kappa_s']) )
        mag_ratio.append(comp_aprox_num.mag_rad_over_mag_tan(i, 0.0, tst['a'], \
                         tst['b'], tst['kappa_s']))
        arg_find.append(comp_aprox_num.arg_find_constant_distortion(i, 0.0, tst['a'], \
                         tst['b'], tst['kappa_s'], -10) )

    plt.plot(radial_vec, mag_ratio, linewidth = 3, label='mag ratio')
    #plt.plot(radial_vec, arg_find, linewidth = 3)

    plt.plot(radial_vec, mag_radial, linewidth = 3, label='mag radial')
    plt.plot(radial_vec, mag_tan, linewidth = 3, label='mag tan')

    plt.plot([0,1.5], [0,0], color = 'Black')
    plt.plot([0,1.5], [.1,.1], color = 'Black')
    plt.plot([0,1.5], [-.1,-.1], color = 'Black')
    plt.legend( loc='lower left' )
    plt.ylim( (-1, 1) )
    plt.show()

################################################################################
def func_test(x, a, b, c):
    return a*x + b * c*x**2

if __name__ == '__main__':
    print '---------------------------------------'
    print 'Module to compute the NFW cross section'
    print '---------------------------------------'

    
    #plot_critical_lines()
    #plot_constant_distortion()
    #cross_section_computation()
    compare_habib()
    #test_PotentialDerivatives()
    #plot_magnifications_thetafix()