# ==================================
# Author: Gabriel Bartosch Caminha - gbcaminha@gmail.com
# ==================================

import math

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
################################################################################
class EnfwAprox():
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
         - x_rad   float: radial coordinate
         - kappa_s float: \kappa_s NFW parameter
        Output:
         - list with file names
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
         - x_rad   float: radial coordinate
         - kappa_s float: \kappa_s NFW parameter
        Output:
         - d(kappa(x))/d(x**2)
        """
        out = -kappa_s/x_rad/x_rad
        return (out, -1e+200)[x_rad < 1E-100]
        # if x_rad < 1E-100 return -1e+200 to avoid the central prifile cusp
        #FIXME
################################################################################
