# ==================================
# Author: Gabriel Bartosch Caminha - gbcaminha@gmail.com
# ==================================

import math

from radial_models import NfwLens

class PellModels(object):
################################################################################
    def __init__(self, radial_lens_object):
        #print 'PellModels: Initializing', radial_lens_object.LensModel
        self.kappa_radial = radial_lens_object.kappa
        self.gamma_radial = radial_lens_object.gamma
################################################################################
    def phi_ell(self, x1, x2, a1, a2):
        if math.fabs(x1) < 1.0E-10:
            return math.pi/2.0
        return math.atan(x2/x1*math.sqrt(a2/a1))
################################################################################
    def A_subs(self, a1, a2):
        return 0.5*(a1 + a2)
################################################################################
    def B_subs(self, a1, a2):
        return 0.5*(a1 - a2)
################################################################################
    def ell_coordenate(self, x_1, x_2, a_1, a_2):
        return math.sqrt(a_1*x_1**2 + a_2*x_2**2)
################################################################################
    def kappa_pell(self, x_1, x_2, a_1, a_2, kappa_s):
        x_e = self.ell_coordenate(x_1, x_2, a_1, a_2)
        p_e = self.phi_ell(x_1, x_2, a_1, a_2)
        c_1 = self.A_subs(a_1, a_2)*self.kappa_radial(x_e, kappa_s)
        c_2 = self.B_subs(a_1, a_2)*self.gamma_radial(x_e, kappa_s)* \
              math.cos(2.0*p_e)
        return c_1 - c_2
################################################################################
    def mag_tan_inv_pol(self, radial, theta, a, b, kappa_s):
        x1 = radial*math.cos(theta)
        x2 = radial*math.sin(theta)
        phiell = self.phi_ell(x1, x2, a1, a2)
        A = self.A_subs(a1, a2)
        B = self.B_subs(a1, a2)
        print "mag_tan_inv_pol not finished"
################################################################################

def kappa_sigma_test(x_1, x_2, ell, b_sis, t0 = 0.0):

    t = math.atan(x_2/x_1)
    r = math.sqrt( x_1**2.0 + x_2**2.0 )
    #t0 = math.pi/4;
    
    cons_1 = (-1.0 + ell*math.cos(2.0*(t - t0)))/math.cos(t)
    cons_2 = ell*(ell - math.cos(2.0*(t - t0))) /math.fabs(math.cos(t))
    cons_3 = 2.0*r*(1.0 - ell*math.cos(2.0*(t - t0)))**(1.5) \
             /math.fabs(math.cos(t))

    return -b_sis*(cons_1 + cons_2)/cons_3




