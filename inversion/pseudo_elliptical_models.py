# ==================================
# Author: Gabriel Bartosch Caminha - gbcaminha@gmail.com
# ==================================

import math

from radial_models import NfwLens

class PnfwLens(object):
    nfw_obj = NfwLens()
    def __init__(self):
        print 'PnfwLens: Initializing PnfwLens'
    def phi_ell(x1, x2, a1, a2):
        arg = x2/x1*math.sqrt(a2/a1)
        return math.atan(arg)
    def A_subs(a1, a2):
        return 0.5*(a1+a2)
    def B_subs(a1, a2):
        return 0.5*(a1 - a2)
    def kappa_pell(self, x1, x2, a1, a2, kappa_s)
        phiell = self.phi_ell(x1, x2, a1, a2)
        A = self.A_subs(a1, a2)
        B = self.B_subs(a1, a2)
        print 'kappa_pell: functions not finished'
    def mag_tan_inv_pol(self, radial, theta, a, b, kappa_s)
        x1 = radial*math.cos(theta)
        x2 = radial*math.sin(theta)
        phiell = self.phi_ell(x1, x2, a1, a2)
        A = self.A_subs(a1, a2)
        B = self.B_subs(a1, a2)



