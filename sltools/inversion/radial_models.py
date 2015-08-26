# ==================================
# Author: Gabriel Bartosch Caminha - gbcaminha@gmail.com
# ==================================

import math

from scipy.optimize import brentq

class NfwLens():
    """
    Class to to compute quantities to the elliptical NFW model
    """
    LensModel = 'NfwLens'
    def __init__(self):
        print 'NfwLens: Initializing NfwLens'
################################################################################
    def testing(self, a, b):
        print 'NfwLens: testing function ', a, b
################################################################################
    def f_nfw(self, x):
        """
        Return the f_function defined at Golse & Kneib (2002)

        Input:
         - x_rad   float: radial coordinate
        Output:
         - g_function
        """
        xm1 = x - 1.0
        if xm1 < -1E-6:
            sqrt_one_x2 = math.sqrt(1.0 - x*x)
            return math.atanh( sqrt_one_x2 )/sqrt_one_x2
        elif xm1 > 1E-6:
            sqrt_x2_one = math.sqrt(x*x - 1.0)
            return math.atan( sqrt_x2_one )/sqrt_x2_one
        else:
            return 1.0
################################################################################
    def f_nfw_prima(self, x):
        return ( 1.0 - x*x*self.f_nfw(x) )/( x*(x*x - 1.0) )
################################################################################
    def g_nfw(self, x_rad):
        """
        Return the g_function defined at Golse & Kneib (2002)

        Input:
         - x_rad   float: radial coordinate
        Output:
         - g_function
        """
        xm1 = x_rad - 1.0
        if xm1 < -1E-6:
            c_1 = 1.0 / math.sqrt(1.0 - x_rad**2.0)
            return math.log(x_rad/2.0) + c_1*math.acosh(1.0/x_rad)
        elif xm1 > 1E-6:
            c_1 = 1.0 / math.sqrt(x_rad**2.0 - 1.0)
            return math.log(x_rad/2.0) + c_1*math.acos(1.0/x_rad)
        else:
            return 0.3068528194400547 #1 + math.log(1./2.)
################################################################################
    def kappa(self, x_rad, kappa_s):
        """
        Return the convergence for the aproximation \kappa_s << 1

        Input:
         - x_rad   float: radial coordinate
         - kappa_s float: \kappa_s NFW parameter
        Output:
         - list with file names
        """
        out = 2.0*kappa_s*(1.0 - self.f_nfw(x_rad))/(x_rad*x_rad-1.0)
        return out
################################################################################
    def kappa_prima(self, x_rad, kappa_s):
        """
        kappa derivative with respect to x^2 ( d(kappa(x))/d(x**2) )

        Input:
         - x_rad   float: radial coordinate
         - kappa_s float: \kappa_s NFW parameter
        Output:
         - d(kappa(x))/d(x**2)
        """
        out = kappa_s*( -self.f_nfw_prima(x_rad)*(x_rad**2.0-1.0) - \
                        (1.0-self.f_nfw(x_rad))*2.0*x_rad  )
        out = out/( (x_rad*x_rad-1.0)**2.0 * x_rad ) 
        return out
################################################################################
    def alpha(self, x_rad, kappa_s, r_s):
        """
        angle deflection for the NFW model. See eq (55) of
        http://arxiv.org/abs/astro-ph/0102341

        Input:
         - x_rad   float : radial coordinate
         - kappa_s float : \kappa_s NFW parameter
         - r_s     float : r_s NFW parameter
        Output:
         - angle deflection
        """
        c1 = 4.0*kappa_s*r_s/x_rad
        c2 = math.log(x_rad/2.0) + self.f_nfw(x_rad)
        return c1*c2
################################################################################
    def alpha_minus_x(self, x_rad, kappa_s, r_s):
        """
        angle deflection minus x_rad, useful to compute the eisteins radius

        Input:
         - x_rad   float : radial coordinate
         - kappa_s float : \kappa_s NFW parameter
         - r_s     float : r_s NFW parameter
        Output:
         - angle deflection
        """
        return self.alpha(x_rad, kappa_s, r_s) - x_rad
################################################################################
    def e_radius(self, kappa_s, r_s):
        """
        computs the einstei radius by numericaly solves alpha - x = 0

        Input:
         - kappa_s float : \kappa_s NFW parameter
         - r_s     float : r_s NFW parameter
        Output:
         - angle deflection
        """
        return brentq(self.alpha_minus_x, a = 1E-4, b = 10.0, \
               args = (kappa_s, r_s))
################################################################################
    def gamma(self, x_rad, kappa_s):
        return 2.0*kappa_s*( 2.0*self.g_nfw(x_rad)/(x_rad**2.0) \
               - self.f_nfw(x_rad) )
################################################################################
################################################################################
class EnfwAprox ():
    """
    Class to to compute quantities to the NFW considering \kappa_s << 1'
    """
    LensModel = 'EnfwAprox'
    def __init__(self):
        print 'EnfwAprox: Initializing EnfwAprox'
################################################################################
    def testing(self, a, b):
        print 'EnfwAprox: testing function ', a, b
################################################################################
    def kappa(self, x_rad, kappa_s):
        """
        Return the convergence for the aproximation \kappa_s << 1

        Input:
         - x_rad   float : radial coordinate
         - kappa_s float : \kappa_s NFW parameter
        Output:
         - nfw convergence
        """
        out = -2.0*kappa_s*( 1.0 + math.log(x_rad/2.0) )
        return (out, 1500)[x_rad < 1E-320]
        # if x_rad < 1E-320 return 1500 to avoid the central prifile cusp
        #FIXME
################################################################################
    def kappa_prima(self, x_rad, kappa_s):
        """
        kappa derivative with respect to x^2 ( d(kappa(x))/d(x**2) )

        Input:
         - x_rad   float : radial coordinate
         - kappa_s float : \kappa_s NFW parameter
        Output:
         - d(kappa(x))/d(x**2)
        """
        out = -kappa_s/x_rad/x_rad
        return (out, -1e+200)[x_rad < 1E-100]
        # if x_rad < 1E-100 return -1e+200 to avoid the central prifile cusp
        #FIXME
################################################################################
class SisLens (object):
    """
    Class to to compute quantities to the NFW considering \kappa_s << 1'
    """
    LensModel = 'SisLens'
    def __init__(self):
        1
        #print 'SisLens: Initializing SisLens'
################################################################################
    def kappa(self, x_rad, b_par):
        """
        Return the convergence for SIS model, kappa = b/(2x)

        Input:
         - x_rad  float : radial coordinate
         - b      float : b SIS parameter
        Output:
         - sis convergence
        """
        #print x_rad
        if x_rad < 1E-6:
            return 1500 #to avoid the central prifile cusp
        return b_par/(x_rad*2.0)
        # if x_rad < 1E-320 return 1500 to avoid the central prifile cusp
        #FIXME
################################################################################
    def gamma(self, x_rad, b_par):
        """
        Return the shear for SIS model, kappa = b/(2x)

        Input:
         - x_rad  float : radial coordinate
         - b      float : b SIS parameter
        Output:
         - sis shear
        """
        if x_rad < 1E-6:
            return 1500 #to avoid the central prifile cusp
        return b_par/(x_rad*2.0)
################################################################################