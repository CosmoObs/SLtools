# ==================================
# Author: Gabriel Bartosch Caminha - gbcaminha@gmail.com
# ==================================
"""
Package to deal with lens sub-halos
"""


import math

import numpy as np

from scipy.integrate import quad

import functions as fu

import constants as cs

class SubHalo(object):

    def __init__(self):
        print 'SubHalo: Initializing SubHalo object'

    def giocoli_mass_function(self, m_sub, m_z0):
        n_m0 = 10.0**(-3.03) # normalozation factos
        #m_z0 = 1E14 # halo final mass
        #m_sub = 1.0 # substructure mass

        xi = m_sub/m_z0

        alpha = -0.9
        beta = 12.2715

        const_1 = n_m0 * m_sub**alpha * math.exp( -beta*xi**3.0 )
        return const_1

    def un_mass_function(self, m_sub):
        m_z0 = 1.0
        xi = m_sub/m_z0

        n_m0 = 0.18 # normalozation factos

class Gao2004(object):

    def __init__(self):
        print "Gao2004: Initializing Gao2004 object"

    def sub_abundance(self, m_sub):
        """
        Sub_hale abundance givem by equation (1)

        Input:
         - m_sub  float : sub halo mass in units of (m_sol/h)
        Output:
         - 10^(-3.2) msub^(-1.9)
        """
        return self.const_1 * (m_sub)**(-1.9)


    const_1 = 10.0**(-3.2)

class Giocoli2011(object):

    def __init__(self):
        print "Giocoli2011: Initializing Giocoli2011"

    def spatial_distribution(self, x, c_vir):
        cons_1 = (1.0 + self.alpha_coma*c_vir)*x**self.beta_coma
        cons_2 = (1.0 + self.alpha_coma*c_vir*x*x)
        return cons_1/cons_2

    def diff_spatial_distribution(self, x, c_vir):

        cons_1 = (1.0 + self.alpha_coma*c_vir)
        cons_2 = (1.0 + self.alpha_coma*c_vir*x*x)
        cons_3 = cons_2*cons_2

        cons_4 = self.beta_coma*x**(self.beta_coma - 1.0)*cons_2 - \
                 2.0*self.alpha_coma*c_vir*x**(self.beta_coma + 1.0)
        return cons_1/cons_3 * cons_4

    alpha_coma = 0.244
    beta_coma = 2.75

    def int_diff_spatial_distribution(self, x, c_vir):
        return quad( self.diff_spatial_distribution, a = 0.0, b =x, \
                     args = (c_vir) )[0]

    def randon_3d_position(self, c_vir, x_min, x_max):
        n_pts = 1000
        radial_vec = fu.randon_function(self.diff_spatial_distribution, n_pts, \
                     x_min, x_max, c_vir)

        theta_vec = np.random.rand(n_pts)*2.0*cs.pi
        phi_vec = np.arccos( 2.0*np.random.rand(n_pts) - 1.0)
        return radial_vec, theta_vec, phi_vec

    def randon_projected_position(self, c_vir, x_min, x_max):
        n_pts = 15000
        radial_vec = fu.randon_function(self.diff_spatial_distribution, n_pts, \
                     x_min, x_max, c_vir)

        theta_vec = np.random.rand(n_pts)*2.0*cs.pi
        phi_vec = np.arccos( 2.0*np.random.rand(n_pts) - 1.0)
        #phi_vec = cs.pi*np.random.rand(n_pts)
        radial_vec = radial_vec * np.sin(phi_vec)
        return radial_vec, theta_vec

    def sub_mass_func(self, m_sub, m_vir, c_vir, z_l):
        alpha = -0.9
        betta = 12.2715
        A = 9.33E-4
        c_men = c_vir #c_vir #FIXME

        cons_1 = A*math.sqrt(1.0 + z_l) * c_men/c_vir * m_sub**alpha
        cons_2 = math.exp( -betta*(m_sub/m_vir)**3.0 )
        return cons_1 * cons_2 



