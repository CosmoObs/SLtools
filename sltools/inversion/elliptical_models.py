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

    def __init__(self, lensmodel_object):
        print 'EllipticalLensIntegrals: Initializing EllipticalLensIntegrals'
        if not ('kappa' in dir(lensmodel_object)):
            print 'EllipticalLensIntegrals: ERROR kappa does not existis in', \
                   lensmodel_object.__name__
        if not ('kappa_prima' in dir(lensmodel_object)):
            print 'EllipticalLensIntegrals: ERROR kappa_prima does not ' \
                  'existis in', lensmodel_object.__name__
        print 'EllipticalLensIntegrals: considering lens model:', \
               lensmodel_object.LensModel
        self.kappa_prima = lensmodel_object.kappa_prima
        self.kappa = lensmodel_object.kappa
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
        c_1 = x_1*x_1/( 1.0 - (1.0 - a_subs*a_subs )*u_subs )
        c_2 = x_2*x_2/( 1.0 - (1.0 - b_subs*b_subs )*u_subs )
        return math.sqrt( u_subs*(c_1 + c_2) )
################################################################################
    def j_argument(self, u_subs, x_1, x_2, a_subs, b_subs, kappa_s, n_index):
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
        c_1 = ( 1.0 - (1.0 - b_subs*b_subs )*u_subs)**(0.5 + n_index)
        c_2 = ( 1.0 - (1.0 - a_subs*a_subs )*u_subs)**(1.5 - n_index)
        return kappa_subst/(c_1*c_2)
################################################################################
    def k_argument(self, u_subs, x_1, x_2, a_subs, b_subs, kappa_s, n_index):
        """
        Argument for the integral K_n

        Input:
         - u_subs  float : new parameter
         - x_1     float : cartesian x1 coordinate (x)
         - x_2     float : cartesian x2 coordinate (y)
         - a_subs  float : parameter to rescale the x1 coordinate
         - b_subs  float : parameter to rescale the x2 coordinate
         - kappa_s float : \kappa_s NFW parameter
         - n_index float : integral index
        Output:
         - argument
        """
        #FIXME - description for the funciton outpt

        kappa_prima_subst = \
         self.kappa_prima(self.m_subst(u_subs, x_1, x_2, a_subs, b_subs), kappa_s)
        c_1 = 1.0/( 1.0 - (1.0 - b_subs*b_subs )*u_subs)**(0.5 + n_index)
        c_2 = 1.0/( 1.0 - (1.0 - a_subs*a_subs )*u_subs)**(2.5 - n_index)
        return u_subs*kappa_prima_subst*c_1*c_2
################################################################################
    def j_func(self, x_1, x_2, a_subs, b_subs, kappa_s, n_index):
        """
        Function J_n defined at 2.29 of GBC MSc Thesis.
        http://cbpfindex.cbpf.br/publication_pdfs/\
        Caminha-Mestrado.2009_06_10_13_00_02.pdf

        Input:
         - x_1     float : cartesian x1 coordinate (x)
         - x_2     float : cartesian x2 coordinate (y)
         - a_subs  float : parameter to rescale the x1 coordinate
         - b_subs  float : parameter to rescale the x2 coordinate
         - kappa_s float : \kappa_s NFW parameter
         - n_index       float : integral index
        Output:
         - integral J_n see eq 2.29 of GBC MSc Thesis
        """
        return integrate.quad(self.j_argument, 0.0, 1.0, \
                              args=(x_1, x_2, a_subs, b_subs, kappa_s, n_index))
################################################################################
    def k_func(self, x_1, x_2, a_subs, b_subs, kappa_s, n_index):
        """
        Function K_n defined at 2.29 of GBC MSc Thesis.
        http://cbpfindex.cbpf.br/publication_pdfs/\
        Caminha-Mestrado.2009_06_10_13_00_02.pdf

        Input:
         - x_1     float : cartesian x1 coordinate (x)
         - x_2     float : cartesian x2 coordinate (y)
         - a_subs  float : parameter to rescale the x1 coordinate
         - b_subs  float : parameter to rescale the x2 coordinate
         - kappa_s float : \kappa_s NFW parameter
         - n_index float : integral index
        Output:
         - integral K_n see eq 2.30 of GBC MSc Thesis
        """
        return integrate.quad(self.k_argument, 0.0, 1.0, \
                              args=(x_1, x_2, a_subs, b_subs, kappa_s, n_index))
################################################################################
################################################################################
class LensPotDrivatives(object):
    """
    Class to compute the lens potential derivatives from the 'integrals' related
    to elliptical lens models
    """
    #FIXME - impruve documentation
    def __init__(self, lensmodel_object):
        print 'LensPotDrivatives: Initializing LensPotDrivatives'
        self.ell_int_obj = EllipticalLensIntegrals(lensmodel_object)
        self.j_func = self.ell_int_obj.j_func
        self.k_func = self.ell_int_obj.k_func
################################################################################
    def pot_x(self, x_1, x_2, a_subs, b_subs, kappa_s):
        """
        Potential derivative in repect to 'x'

        Input:
         - x_1     float : cartesian x1 coordinate (x)
         - x_2     float : cartesian x2 coordinate (y)
         - a_subs  float : parameter to rescale the x1 coordinate
         - b_subs  float : parameter to rescale the x2 coordinate
         - kappa_s float : \kappa_s NFW parameter
        Output:
         - d(pot)/dx
        """
        return \
           a_subs*b_subs*x_1*self.j_func(x_1, x_2, a_subs, b_subs, kappa_s, 0.0)[0]
################################################################################
    def pot_y(self, x_1, x_2, a_subs, b_subs, kappa_s):
        """
        Potential derivative in repect to 'y'

        Input:
         - x_1     float : cartesian x1 coordinate (x)
         - x_2     float : cartesian x2 coordinate (y)
         - a_subs  float : parameter to rescale the x1 coordinate
         - b_subs  float : parameter to rescale the x2 coordinate
         - kappa_s float : \kappa_s NFW parameter
        Output:
         - d(pot)/dy
        """
        return \
           a_subs*b_subs*x_2*self.j_func(x_1, x_2, a_subs, b_subs, kappa_s, 1.0)[0]
################################################################################
    def pot_xx(self, x_1, x_2, a_subs, b_subs, kappa_s):
        """
        Potential second derivative in repect to 'x'

        Input:
         - x_1     float : cartesian x1 coordinate (x)
         - x_2     float : cartesian x2 coordinate (y)
         - a_subs  float : parameter to rescale the x1 coordinate
         - b_subs  float : parameter to rescale the x2 coordinate
         - kappa_s float : \kappa_s NFW parameter
        Output:
         - d^2(pot)/dx^2
        """
        c_1 = 2.0*a_subs*b_subs*x_1*x_1* \
             self.k_func(x_1, x_2, a_subs, b_subs, kappa_s, 0.0)[0]
        c_2 = a_subs*b_subs* \
             self.j_func(x_1, x_2, a_subs, b_subs, kappa_s, 0.0)[0]
        return c_1 + c_2
################################################################################
    def pot_yy(self, x_1, x_2, a_subs, b_subs, kappa_s):
        """
        Potential second derivative in repect to 'y'

        Input:
         - x_1     float : cartesian x1 coordinate (x)
         - x_2     float : cartesian x2 coordinate (y)
         - a_subs  float : parameter to rescale the x1 coordinate
         - b_subs  float : parameter to rescale the x2 coordinate
         - kappa_s float : \kappa_s NFW parameter
        Output:
         - d^2(pot)/dy^2
        """
        c_1 = 2.0*a_subs*b_subs*x_2*x_2* \
             self.k_func(x_1, x_2, a_subs, b_subs, kappa_s, 2.0)[0]
        c_2 = a_subs*b_subs*\
             self.j_func(x_1, x_2, a_subs, b_subs, kappa_s, 1.0)[0]
        return c_1 + c_2
################################################################################
    def pot_xy(self, x_1, x_2, a_subs, b_subs, kappa_s):
        """
        Potential second derivative in repect to 'x' and 'y'

        Input:
         - x_1     float : cartesian x1 coordinate (x)
         - x_2     float : cartesian x2 coordinate (y)
         - a_subs  float : parameter to rescale the x1 coordinate
         - b_subs  float : parameter to rescale the x2 coordinate
         - kappa_s float : \kappa_s NFW parameter
        Output:
         - d^2(pot)/dxdy
        """
        return 2.0*a_subs*b_subs*x_1*x_2* \
               self.k_func(x_1, x_2, a_subs, b_subs, kappa_s, 1.0)[0]
################################################################################
class KappaGammaEll(object):
    """
    Class to compute the convergence and the shear using a LensModel object
    """
    def __init__(self, lensmodel_object):
        self.pot_derivatives_obj = LensPotDrivatives(lensmodel_object)
        self.pot_xx = self.pot_derivatives_obj.pot_xx
        self.pot_yy = self.pot_derivatives_obj.pot_yy
        self.pot_xy = self.pot_derivatives_obj.pot_xy
        self.kappa_analitic = lensmodel_object.kappa
################################################################################
    def kappa_ell(self, x_1, x_2, a_subs, b_subs, kappa_s):
        """
        Compute the elliptical convergance from the radial model in LensModel
        object

        Input:
         - x_1     float : cartesian x1 coordinate (x)
         - x_2     float : cartesian x2 coordinate (y)
         - a_subs  float : parameter to rescale the x1 coordinate
         - b_subs  float : parameter to rescale the x2 coordinate
         - kappa_s float : \kappa_s NFW parameter
        Output:
         - kappa_radial(math.sqrt((x_1/a_subs)**2 + (x_2/b_subs)**2))
        """
        x_ell = math.sqrt((x_1/a_subs)**2 + (x_2/b_subs)**2)
        return self.kappa_analitic(x_ell, kappa_s)
################################################################################
    def kappa_numeric(self, x_1, x_2, a_subs, b_subs, kappa_s):
        """
        Compute the elliptical convergance from the potential derivatives.
        Usefull only for numerical checks

        Input:
         - x_1     float : cartesian x1 coordinate (x)
         - x_2     float : cartesian x2 coordinate (y)
         - a_subs  float : parameter to rescale the x1 coordinate
         - b_subs  float : parameter to rescale the x2 coordinate
         - kappa_s float : \kappa_s NFW parameter
        Output:
         - 0.5*(pot_xx + pot_yy)
        """
        c_1 = self.pot_xx(x_1, x_2, a_subs, b_subs, kappa_s) + \
             self.pot_yy(x_1, x_2, a_subs, b_subs, kappa_s)
        return c_1/2.0
################################################################################
    def gamma1(self, x_1, x_2, a_subs, b_subs, kappa_s):
        """
        First shear component

        Input:
         - x_1     float : cartesian x1 coordinate (x)
         - x_2     float : cartesian x2 coordinate (y)
         - a_subs  float : parameter to rescale the x1 coordinate
         - b_subs  float : parameter to rescale the x2 coordinate
         - kappa_s float : \kappa_s NFW parameter
        Output:
         - 0.5*(pot_xx - pot_yy)
        """
        return 0.5*( self.pot_xx(x_1, x_2, a_subs, b_subs, kappa_s) - \
                      self.pot_yy(x_1, x_2, a_subs, b_subs, kappa_s) )
################################################################################
    def gamma2(self, x_1, x_2, a_subs, b_subs, kappa_s):
        """
        Second shear component

        Input:
         - x_1     float : cartesian x1 coordinate (x)
         - x_2     float : cartesian x2 coordinate (y)
         - a_subs  float : parameter to rescale the x1 coordinate
         - b_subs  float : parameter to rescale the x2 coordinate
         - kappa_s float : \kappa_s NFW parameter
        Output:
         - pot_xy
        """
        return self.pot_xy(x_1, x_2, a_subs, b_subs, kappa_s)
################################################################################
    def gamma(self, x_1, x_2, a_subs, b_subs, kappa_s):
        """
        Total shear

        Input:
         - x_1     float : cartesian x1 coordinate (x)
         - x_2     float : cartesian x2 coordinate (y)
         - a_subs  float : parameter to rescale the x1 coordinate
         - b_subs  float : parameter to rescale the x2 coordinate
         - kappa_s float : \kappa_s NFW parameter
        Output:
         - sqrt(gamma1**2.0 + gamma2**2.0)
        """
        gamma1 = self.gamma1(x_1, x_2, a_subs, b_subs, kappa_s)
        gamma2 = self.gamma2(x_1, x_2, a_subs, b_subs, kappa_s)
        return math.sqrt(gamma1**2.0 + gamma2**2.0)
################################################################################
    def mag_total(self, x_1, x_2, a_subs, b_subs, kappa_s):
        """
        Total magnification

        Input:
         - x_1     float : cartesian x1 coordinate (x)
         - x_2     float : cartesian x2 coordinate (y)
         - a_subs  float : parameter to rescale the x1 coordinate
         - b_subs  float : parameter to rescale the x2 coordinate
         - kappa_s float : \kappa_s NFW parameter
        Output:
         - 1/( (1 - kappa)**2 - gamma**2 )
        """
        kappa = self.kappa_ell(x_1, x_2, a_subs, b_subs, kappa_s)
        gamma = self.gamma(x_1, x_2, a_subs, b_subs, kappa_s)
        out = 1.0/( (1.0 - kappa)**2 - gamma**2  )
        return out
################################################################################
    def mag_tan(self, x_1, x_2, a_subs, b_subs, kappa_s):
        """
        Tangential magnification

        Input:
         - x_1     float : cartesian x1 coordinate (x)
         - x_2     float : cartesian x2 coordinate (y)
         - a_subs  float : parameter to rescale the x1 coordinate
         - b_subs  float : parameter to rescale the x2 coordinate
         - kappa_s float : \kappa_s NFW parameter
        Output:
         - 1/( 1 - kappa - gamma )
        """
        kappa = self.kappa_ell(x_1, x_2, a_subs, b_subs, kappa_s)
        gamma = self.gamma(x_1, x_2, a_subs, b_subs, kappa_s)
        out = 1.0/(1.0 - kappa - gamma)
        return out
################################################################################
    def mag_tan_inv_polar(self, radial, theta, a_subs, b_subs, kappa_s):
        """
        Inverse of tangential magnification in polar coordinates. It is useful
        for some applications.

        Input:
         - radial  float : radial coordinate
         - theta   float : polar coordinate
         - a_subs  float : parameter to rescale the x1 coordinate
         - b_subs  float : parameter to rescale the x2 coordinate
         - kappa_s float : \kappa_s NFW parameter
        Output:
         - ( 1 - kappa - gamma )
        """
        x_1 = radial*math.cos(theta)
        x_2 = radial*math.sin(theta)
        kappa = self.kappa_ell(x_1, x_2, a_subs, b_subs, kappa_s)
        gamma = self.gamma(x_1, x_2, a_subs, b_subs, kappa_s)
        out = (1.0 - kappa - gamma)
        return out
################################################################################
    def mag_rad(self, x_1, x_2, a_subs, b_subs, kappa_s):
        """
        Radial magnification

        Input:
         - x_1     float : cartesian x1 coordinate (x)
         - x_2     float : cartesian x2 coordinate (y)
         - a_subs  float : parameter to rescale the x1 coordinate
         - b_subs  float : parameter to rescale the x2 coordinate
         - kappa_s float : \kappa_s NFW parameter
        Output:
         - 1/( 1 - kappa + gamma )
        """
        kappa = self.kappa_ell(x_1, x_2, a_subs, b_subs, kappa_s)
        gamma = self.gamma(x_1, x_2, a_subs, b_subs, kappa_s)
        out = 1.0/(1.0 - kappa + gamma)
        return out
################################################################################
    def mag_rad_inv_polar(self, radial, theta, a_subs, b_subs, kappa_s):
        """
        Inverse of radial magnification in polar coordinates. It is useful 
        for some applications.

        Input:
         - radial  float : radial coordinate
         - theta   float : polar coordinate
         - a_subs  float : parameter to rescale the x1 coordinate
         - b_subs  float : parameter to rescale the x2 coordinate
         - kappa_s float : \kappa_s NFW parameter
        Output:
         - ( 1 - kappa + gamma )
        """
        x_1 = radial*math.cos(theta)
        x_2 = radial*math.sin(theta)
        kappa = self.kappa_ell(x_1, x_2, a_subs, b_subs, kappa_s)
        gamma = self.gamma(x_1, x_2, a_subs, b_subs, kappa_s)
        out = (1.0 - kappa + gamma)
        return out
################################################################################
    def mag_all(self, x_1, x_2, a_subs, b_subs, kappa_s):
        """
        Computes the total, tangential and radial magnifications. It is faster
        use this functions to get the three magnifications than the three
        functions for each one. 

        Input:
         - x_1     float : cartesian x1 coordinate (x)
         - x_2     float : cartesian x2 coordinate (y)
         - a_subs  float : parameter to rescale the x1 coordinate
         - b_subs  float : parameter to rescale the x2 coordinate
         - kappa_s float : \kappa_s NFW parameter
        Output:
         - mag_tot, mag_tan, mag_rad
        """
        kappa = self.kappa_ell(x_1, x_2, a_subs, b_subs, kappa_s)
        gamma = self.gamma(x_1, x_2, a_subs, b_subs, kappa_s)

        magtan = 1.0/(1.0 - kappa - gamma)
        magrad = 1.0/(1.0 - kappa + gamma)
        magtot = magtan*magrad        

        return magtot, magtan, magrad
################################################################################
    def mag_inv_all(self, x_1, x_2, a_subs, b_subs, kappa_s):
        """
        Computes the inverse of the total, tangential and radial magnifications.
        It is faster use this functions to get the three magnifications than the
        three functions for each one. 

        Input:
         - x_1     float : cartesian x1 coordinate (x)
         - x_2     float : cartesian x2 coordinate (y)
         - a_subs  float : parameter to rescale the x1 coordinate
         - b_subs  float : parameter to rescale the x2 coordinate
         - kappa_s float : \kappa_s NFW parameter
        Output:
         - 1/mag_tot, 1/mag_tan, 1/mag_rad
        """
        kappa = self.kappa_ell(x_1, x_2, a_subs, b_subs, kappa_s)
        gamma = self.gamma(x_1, x_2, a_subs, b_subs, kappa_s)

        magtan_inv = (1.0 - kappa - gamma)
        magrad_inv = (1.0 - kappa + gamma)
        magtot_inv = magtan_inv*magrad_inv

        return magtot_inv, magtan_inv, magrad_inv
################################################################################
class LensComputation(object):
    """
    Class to compute complex quantities related to the lensing
    """
    def __init__(self, kappagamma_object):
        self.mag_tan_inv_polar = kappagamma_object.mag_tan_inv_polar
        self.mag_rad_inv_polar = kappagamma_object.mag_rad_inv_polar
        self.mag_inv_all = kappagamma_object.mag_inv_all
    def find_cc_tan(self, theta, a_subs, b_subs, kappa_s):
        """
        Find the radial coordinate for the tangential critical curve, for a 
        fixed angular coordinate

        Input:
         - theta   float : polar coordinate
         - a_subs  float : parameter to rescale the x1 coordinate
         - b_subs  float : parameter to rescale the x2 coordinate
         - kappa_s float : \kappa_s NFW parameter
        Output:
         - r_tangential_critical_curve
        """
        out = brentq(self.mag_tan_inv_polar, a = 1E-4, b = 10.0, \
                     args = (theta, a_subs, b_subs, kappa_s), full_output=True,\
                     xtol = 1E-12)
        return out
################################################################################
    def find_cc_rad(self, theta, a_subs, b_subs, kappa_s):
        """
        Find the radial coordinate for the radial critical curve, for a 
        fixed angular coordinate

        Input:
         - theta   float : polar coordinate
         - a_subs  float : parameter to rescale the x1 coordinate
         - b_subs  float : parameter to rescale the x2 coordinate
         - kappa_s float : \kappa_s NFW parameter
        Output:
         - r_radial_critical_curve
        """
        out = brentq(self.mag_rad_inv_polar, a = 1E-4, b = 10.0, \
                     args = (theta, a_subs, b_subs, kappa_s), full_output=True,\
                     xtol = 1E-12)
        return out
################################################################################
    def mag_rad_over_mag_tan(self, radial, theta, a_subs, b_subs, kappa_s):
        """
        Computes the ratio between the radial and tangential magnifications.
        Usefull to compute the curves of constant distortion

        Input:
         - radial  float : radial coordinate
         - theta   float : polar coordinate
         - a_subs  float : parameter to rescale the x1 coordinate
         - b_subs  float : parameter to rescale the x2 coordinate
         - kappa_s float : \kappa_s NFW parameter
        Output:
         - magnification_radial / magnification_tangential
        """
        c_1 = self.mag_tan_inv_polar(radial, theta, a_subs, b_subs, kappa_s)
        c_2 = self.mag_rad_inv_polar(radial, theta, a_subs, b_subs, kappa_s)
        return c_1/c_2
################################################################################
    def arg_find_constant_distortion(self, radial, theta, a_subs, b_subs, \
                                     kappa_s, raz):
        """
        Function to be used by the root finder to find the curves of constant
        distortion

        Input:
         - radial  float : radial coordinate
         - theta   float : polar coordinate
         - a_subs  float : parameter to rescale the x1 coordinate
         - b_subs  float : parameter to rescale the x2 coordinate
         - kappa_s float : \kappa_s NFW parameter
         - raz     float : tangential over radial magnification
        Output:
         - magnification_radial / magnification_tangential
        """
        mag_ratio = self.mag_rad_over_mag_tan(radial, theta, a_subs, b_subs, kappa_s)
        return  mag_ratio - 1.0/raz
################################################################################
    def find_constant_distortion(self, theta, a_subs, b_subs, kappa_s, raz):
        """
        Find the radial coordinate for the curves of constant distortion. Since
        this function computes the radial and tangential critical curves, this
        quantities are outpts too.

        Input:
         - theta   float : polar coordinate
         - a_subs  float : parameter to rescale the x1 coordinate
         - b_subs  float : parameter to rescale the x2 coordinate
         - kappa_s float : \kappa_s NFW parameter
         - raz     float : tangential over radial magnification
        Output:
         - r_tangential_cc, r_radial_cc, r_constdist_minus, r_constdist_plus 
        """
        cc_rad = self.find_cc_rad(theta, a_subs, b_subs, kappa_s)[0]
        cc_tan = self.find_cc_tan(theta, a_subs, b_subs, kappa_s)[0]
        c_1 = brentq(self.arg_find_constant_distortion, a = cc_rad+1E-5, b = cc_tan, \
                    args = (theta, a_subs, b_subs, kappa_s, -raz), \
                    full_output=True, xtol = 1E-12)[0]
        c_2 = brentq(self.arg_find_constant_distortion, a = cc_tan, b = 100, \
                     args = (theta, a_subs, b_subs, kappa_s, raz), \
                     full_output=True, xtol = 1E-12)[0]
        return cc_tan, cc_rad, c_1, c_2
################################################################################
    def arg_sigma_rad(self, radial, theta, a_subs, b_subs, kappa_s):
        """
        Function to be integrated along the radial coordinate to compute the 
        cross-section

        Input:
         - radial  float : radial coordinate
         - theta   float : polar coordinate
         - a_subs  float : parameter to rescale the x1 coordinate
         - b_subs  float : parameter to rescale the x2 coordinate
         - kappa_s float : \kappa_s NFW parameter
        Output:
         - 1/total_magnifications * radial
        """
        x_1 = radial*math.cos(theta)
        x_2 = radial*math.sin(theta)
        out = self.mag_inv_all(x_1, x_2, a_subs, b_subs, kappa_s)[0]
        return math.fabs(out)*radial
################################################################################
    def sigma_radial(self, theta, a_subs, b_subs, kappa_s, raz):
        """
        Function to be integrated along the angular coordinate to compute the 
        cross-section

        Input:
         - theta   float : polar coordinate
         - a_subs  float : parameter to rescale the x1 coordinate
         - b_subs  float : parameter to rescale the x2 coordinate
         - kappa_s float : \kappa_s NFW parameter
         - raz     float : tangential over radial magnification
        Output:
         - integrate(arg_sigma_rad)
        """
        cc_tan, cc_rad, r_min, r_max = \
              self.find_constant_distortion(theta, a_subs, b_subs, kappa_s, raz)
        del(cc_rad)
        return integrate.quad( self.arg_sigma_rad, a = r_min, b = r_max, \
                               args = (theta, a_subs, b_subs, kappa_s), \
                               points= [cc_tan])[0]
################################################################################
    def sigma(self, a_subs, b_subs, kappa_s, raz):
        """
        Computes the deformations cross-section.

        Input:
         - a_subs  float : parameter to rescale the x1 coordinate
         - b_subs  float : parameter to rescale the x2 coordinate
         - kappa_s float : \kappa_s NFW parameter
         - raz     float : tangential over radial magnification
        Output:
         - deformation cross-sections, see eq 2.31 of GBC MSc Thesis
        """
        return integrate.quad( self.sigma_radial, a = 0.0, b = cn.pi/2.0, \
                               args = (a_subs, b_subs, kappa_s, raz) )[0]*4.0
################################################################################
