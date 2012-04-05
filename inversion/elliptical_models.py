# ==================================
# Author: Gabriel Bartosch Caminha - gbcaminha@gmail.com
# ==================================
"""
Package to compute NFW model quantities
"""

import math

import constants as cn

from scipy import integrate
from scipy.optimize import brentq



class EllipticalLensIntegrals(object):
    """
    Class to compute the 'integrals' related to elliptical lens models
    """
    #FIXME - impruve documentation
    Name = 'EllipticalLensIntegrals'
    def __init__(self, LensModel):
        print 'EllipticalLensIntegrals: Initializing EllipticalLensIntegrals'
        if not ('kappa' in dir(LensModel)):
            print 'EllipticalLensIntegrals: ERROR kappa does not existis in', \
                   LensModel.__name__
        if not ('kappa_prima' in dir(LensModel)):
            print 'EllipticalLensIntegrals: ERROR kappa_prima does not ' \
                  'existis in', LensModel.__name__
        print 'EllipticalLensIntegrals: considering lens model:',\
               LensModel.LensModel
        self.kappa_prima = LensModel.kappa_prima
        self.kappa = LensModel.kappa
################################################################################
    def testing(self, a_in, b_in):
        print 'EnfwLens: testing function ', a_in, b_in
################################################################################
    def m_subst(self, u_subs, x_1, x_2, a_subs, b_subs):
        """
        Function to make the coordinate substitution to compute the J_n and K_n
        integrals

        Input:
         - u_subs  float: new parameter
         - x_1 float: cartesian x1 coordinate (x)
         - x_2 float: cartesian x2 coordinate (y)
         - a_subs  float: parameter to rescale the x1 coordinate
         - b_subs  float: parameter to rescale the x2 coordinate
        Output:
         - sqrt(u*(x1*x1/(1.0 - (1.0 - a*a )*u) + x2*x2/(1.0 - (1.0 - b*b )*u)))
        """
        c1 = x_1*x_1/( 1.0 - (1.0 - a_subs*a_subs )*u_subs )
        c2 = x_2*x_2/( 1.0 - (1.0 - b_subs*b_subs )*u_subs )
        return math.sqrt( u_subs*(c1 + c2) )
################################################################################
    def j_argument(self, u_subs, x_1, x_2, a_subs, b_subs, kappa_s, n):
        """
        Argument for the integral J_n

        Input:
         - u_subs       float: new parameter
         - x_1      float: cartesian x1 coordinate (x)
         - x_2      float: cartesian x2 coordinate (y)
         - a_subs       float: parameter to rescale the x1 coordinate
         - b_subs       float: parameter to rescale the x2 coordinate
         - kappa_s float: \kappa_s NFW parameter
         - n       float: integral index
        Output:
         - argument
        """
        #FIXME - description for the funciton outpt
        kappa_subst = \
               self.kappa(self.m_subst(u_subs, x_1, x_2, a_subs, b_subs), kappa_s)
        c1 = ( 1.0 - (1.0 - b_subs*b_subs )*u_subs)**(0.5 + n)
        c2 = ( 1.0 - (1.0 - a_subs*a_subs )*u_subs)**(1.5 - n)
        return kappa_subst/(c1*c2)
################################################################################
    def k_argument(self, u_subs, x_1, x_2, a_subs, b_subs, kappa_s, n):
        """
        Argument for the integral K_n

        Input:
         - u_subs       float: new parameter
         - x_1      float: cartesian x1 coordinate (x)
         - x_2      float: cartesian x2 coordinate (y)
         - a_subs       float: parameter to rescale the x1 coordinate
         - b_subs       float: parameter to rescale the x2 coordinate
         - kappa_s float: \kappa_s NFW parameter
         - n       float: integral index
        Output:
         - argument
        """
        #FIXME - description for the funciton outpt
        kappa_prima_subst = \
         self.kappa_prima(self.m_subst(u_subs, x_1, x_2, a_subs, b_subs), kappa_s)
        c1 = 1.0/( 1.0 - (1.0 - b_subs*b_subs )*u_subs)**(0.5 + n)
        c2 = 1.0/( 1.0 - (1.0 - a_subs*a_subs )*u_subs)**(2.5 - n)
        return u_subs*kappa_prima_subst*c1*c2
################################################################################
    def j_func(self, x_1, x_2, a_subs, b_subs, kappa_s, n):
        return integrate.quad(self.j_argument, 0.0, 1.0, \
                              args=(x_1, x_2, a_subs, b_subs, kappa_s, n))
################################################################################
    def k_func(self, x_1, x_2, a_subs, b_subs, kappa_s, n):
        return integrate.quad(self.k_argument, 0.0, 1.0, \
                              args=(x_1, x_2, a_subs, b_subs, kappa_s, n))
################################################################################
################################################################################
class LensPotDrivatives(object):
    """
    Class to compute the lens potential derivatives from the 'integrals' related
    to elliptical lens models
    """
    #FIXME - impruve documentation
    def __init__(self, LensModel):
        print 'LensPotDrivatives: Initializing LensPotDrivatives'
        self.ell_int_obj = EllipticalLensIntegrals(LensModel)
        self.j_func = self.ell_int_obj.j_func
        self.k_func = self.ell_int_obj.k_func
################################################################################
    def pot_x(self, x_1, x_2, a_subs, b_subs, kappa_s):
        return \
           a_subs*b_subs*x_1*self.j_func(x_1, x_2, a_subs, b_subs, kappa_s, 0.0)[0]
################################################################################
    def pot_y(self, x_1, x_2, a_subs, b_subs, kappa_s):
        return \
           a_subs*b_subs*x_2*self.j_func(x_1, x_2, a_subs, b_subs, kappa_s, 1.0)[0]
################################################################################
    def pot_xx(self, x_1, x_2, a_subs, b_subs, kappa_s):
        c1 = 2.0*a_subs*b_subs*x_1*x_1* \
             self.k_func(x_1, x_2, a_subs, b_subs, kappa_s, 0.0)[0]
        c2 = a_subs*b_subs* \
             self.j_func(x_1, x_2, a_subs, b_subs, kappa_s, 0.0)[0]
        return c1 + c2
################################################################################
    def pot_yy(self, x_1, x_2, a_subs, b_subs, kappa_s):
        c1 = 2.0*a_subs*b_subs*x_2*x_2* \
             self.k_func(x_1, x_2, a_subs, b_subs, kappa_s, 2.0)[0]
        c2 = a_subs*b_subs*\
             self.j_func(x_1, x_2, a_subs, b_subs, kappa_s, 1.0)[0]
        return c1 + c2
################################################################################
    def pot_xy(self, x_1, x_2, a_subs, b_subs, kappa_s):
        return 2.0*a_subs*b_subs*x_1*x_2* \
               self.k_func(x_1, x_2, a_subs, b_subs, kappa_s, 1.0)[0]
################################################################################
################################################################################
class KappaGammaEll(object):
    def __init__(self, LensModel):
        self.pot_derivatives_obj = LensPotDrivatives(LensModel)
        self.pot_xx = self.pot_derivatives_obj.pot_xx
        self.pot_yy = self.pot_derivatives_obj.pot_yy
        self.pot_xy = self.pot_derivatives_obj.pot_xy
        self.kappa_analitic = LensModel.kappa
################################################################################
    def kappa_ell(self, x_1, x_2, a_subs, b_subs, kappa_s):
        x_ell = math.sqrt((x_1/a_subs)**2 + (x_2/b_subs)**2)
        return self.kappa_analitic(x_ell, kappa_s)
################################################################################
    def kappa_numeric(self, x_1, x_2, a_subs, b_subs, kappa_s):
        c1 = self.pot_xx(x_1, x_2, a_subs, b_subs, kappa_s) + \
             self.pot_yy(x_1, x_2, a_subs, b_subs, kappa_s)
        return c1/2.0
################################################################################
    def gamma1(self, x_1, x_2, a_subs, b_subs, kappa_s):
        return 0.5*( self.pot_xx(x_1, x_2, a_subs, b_subs, kappa_s) - \
                      self.pot_yy(x_1, x_2, a_subs, b_subs, kappa_s) )
################################################################################
    def gamma2(self, x_1, x_2, a_subs, b_subs, kappa_s):
        return self.pot_xy(x_1, x_2, a_subs, b_subs, kappa_s)
################################################################################
    def gamma(self, x_1, x_2, a_subs, b_subs, kappa_s):
        gamma1 = self.gamma1(x_1, x_2, a_subs, b_subs, kappa_s)
        gamma2 = self.gamma2(x_1, x_2, a_subs, b_subs, kappa_s)
        return math.sqrt(gamma1**2.0 + gamma2**2.0)
################################################################################
    def mag_total(self, x_1, x_2, a_subs, b_subs, kappa_s):
        kappa = self.kappa_ell(x_1, x_2, a_subs, b_subs, kappa_s)
        gamma = self.gamma(x_1, x_2, a_subs, b_subs, kappa_s)
        out =1.0/( (1.0 - kappa)**2 - gamma**2  )
        return out
################################################################################
    def mag_tan(self, x_1, x_2, a_subs, b_subs, kappa_s):
        kappa = self.kappa_ell(x_1, x_2, a_subs, b_subs, kappa_s)
        gamma = self.gamma(x_1, x_2, a_subs, b_subs, kappa_s)
        out = 1.0/(1.0 - kappa - gamma)
        return out
################################################################################
    def mag_tan_inv_polar(self, radial, theta, a_subs, b_subs, kappa_s):
        x_1 = radial*math.cos(theta)
        x_2 = radial*math.sin(theta)
        kappa = self.kappa_ell(x_1, x_2, a_subs, b_subs, kappa_s)
        gamma = self.gamma(x_1, x_2, a_subs, b_subs, kappa_s)
        out = (1.0 - kappa - gamma)
        return out
################################################################################
    def mag_rad(self, x_1, x_2, a_subs, b_subs, kappa_s):
        kappa = self.kappa_in(x_1, x_2, a_subs, b_subs, kappa_s)
        gamma = self.gamma(x_1, x_2, a_subs, b_subs, kappa_s)
        out = 1.0/(1.0 - kappa + gamma)
        return out
################################################################################
    def mag_rad_inv_polar(self, radial, theta, a_subs, b_subs, kappa_s):
        x_1 = radial*math.cos(theta)
        x_2 = radial*math.sin(theta)
        kappa = self.kappa_ell(x_1, x_2, a_subs, b_subs, kappa_s)
        gamma = self.gamma(x_1, x_2, a_subs, b_subs, kappa_s)
        out = (1.0 - kappa + gamma)
        return out
################################################################################
    def mag_all(self, x_1, x_2, a_subs, b_subs, kappa_s):
        kappa = self.kappa_ell(x_1, x_2, a_subs, b_subs, kappa_s)
        gamma = self.gamma(x_1, x_2, a_subs, b_subs, kappa_s)

        magtan = 1.0/(1.0 - kappa - gamma)
        magrad = 1.0/(1.0 - kappa + gamma)
        magtot = magtan*magrad        

        return magtot, magtan, magrad
################################################################################
    def mag_inv_all(self, x_1, x_2, a_subs, b_subs, kappa_s):
        kappa = self.kappa_ell(x_1, x_2, a_subs, b_subs, kappa_s)
        gamma = self.gamma(x_1, x_2, a_subs, b_subs, kappa_s)

        magtan_inv = (1.0 - kappa - gamma)
        magrad_inv = (1.0 - kappa + gamma)
        magtot_inv = magtan_inv*magrad_inv

        return magtot_inv, magtan_inv, magrad_inv
################################################################################
class LensComputation(object):
    def __init__(self, KappaGamma):
        self.mag_tan_inv_polar = KappaGamma.mag_tan_inv_polar
        self.mag_rad_inv_polar = KappaGamma.mag_rad_inv_polar
        self.mag_inv_all = KappaGamma.mag_inv_all
    def find_cc_tan(self, theta, a_subs, b_subs, kappa_s):
        out = brentq(self.mag_tan_inv_polar, a = 1E-4, b = 10.0, \
                     args = (theta, a_subs, b_subs, kappa_s), full_output=True,\
                     xtol = 1E-12)
        return out
################################################################################
    def find_cc_rad(self, theta, a_subs, b_subs, kappa_s):
        out = brentq(self.mag_rad_inv_polar, a = 1E-4, b = 10.0, \
                     args = (theta, a_subs, b_subs, kappa_s), full_output=True,\
                     xtol = 1E-12)
        return out
################################################################################
    def mag_rad_over_mag_tan(self, radial, theta, a_subs, b_subs, kappa_s):
        c1 = self.mag_tan_inv_polar(radial, theta, a_subs, b_subs, kappa_s)
        c2 = self.mag_rad_inv_polar(radial, theta, a_subs, b_subs, kappa_s)
        return c1/c2
################################################################################
    def arg_find_constant_distortion(self, radial, theta, a_subs, b_subs, \
                                     kappa_s, raz):
        c1 = self.mag_rad_over_mag_tan(radial, theta, a_subs, b_subs, kappa_s)
        return  c1 - raz


################################################################################
    def find_constant_distortion(self, theta, a_subs, b_subs, kappa_s, raz):
        cc_rad = self.find_cc_rad(theta, a_subs, b_subs, kappa_s)[0]
        cc_tan = self.find_cc_tan(theta, a_subs, b_subs, kappa_s)[0]
        c1 = brentq(self.arg_find_constant_distortion, a = cc_rad+1E-5, b = cc_tan, \
                    args = (theta, a_subs, b_subs, kappa_s, -raz), \
                    full_output=True, xtol = 1E-12)[0]
        c2 = brentq(self.arg_find_constant_distortion, a = cc_tan, b = 100, \
                     args = (theta, a_subs, b_subs, kappa_s, raz), \
                     full_output=True, xtol = 1E-12)[0]
        #print  cc_rad, c1, cc_tan, c2
        return cc_tan, cc_rad, c1, c2
################################################################################
    def arg_sigma_rad(self, radial, theta, a_subs, b_subs, kappa_s):
        x_1 = radial*math.cos(theta)
        x_2 = radial*math.sin(theta)
        out = self.mag_inv_all(x_1, x_2, a_subs, b_subs, kappa_s)[0]

        #out = out if out > 0.0 else -out
        return math.fabs(out)*radial
################################################################################
    def sigma_radial(self, theta, a_subs, b_subs, kappa_s, raz):
        cc_tan, cc_rad, c1, c2 = \
              self.find_constant_distortion(theta, a_subs, b_subs, kappa_s, raz)
        return integrate.quad( self.arg_sigma_rad, a = c1, b = c2, \
                               args = (theta, a_subs, b_subs, kappa_s), \
                               points= [cc_tan])[0]
################################################################################
    def sigma(self, a_subs, b_subs, kappa_s, raz):
        return integrate.quad( self.sigma_radial, a = 0.0, b = cn.pi/2.0, \
                               args = (a_subs, b_subs, kappa_s, raz) )[0]*4.0
################################################################################
