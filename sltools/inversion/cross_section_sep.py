#!/usr/bin/env python
# ==================================
# Author: Gabriel Bartosch Caminha - gbcaminha@gmail.com
# ==================================

import math

from scipy import integrate

import constants as cn

import numpy as np

def kappa_radial(radial, kappa_s):
    return -2.0*kappa_s*( 1.0 + math.log(radial/2.0) )
################################################################################
def g_function(u_subs, phi, a_subs, b_subs):
    const1 = (np.cos(phi)**2.0)/( 1.0 - (1.0 - a_subs**2.0 )*u_subs )
    const2 = (np.sin(phi)**2.0)/( 1.0 - (1.0 - b_subs**2.0 )*u_subs )
    #print math.sqrt(u_subs*const1*const2), 'AA'
    return np.sqrt(u_subs*( const1 + const2 ) )
################################################################################
def a0_function(a_subs, b_subs):
    return 2.0/a_subs/(a_subs+b_subs)
################################################################################
def a1_function(a_subs, b_subs):
    return 2.0/b_subs/(a_subs+b_subs)
################################################################################
def b_argument(u_subs, phi, a_subs, b_subs, n_ind):
    const1 = 1.0/( 1.0 - (1.0 - b_subs**2.0 )*u_subs )**(0.5 + float(n_ind))
    const2 = 1.0/( 1.0 - (1.0 - a_subs**2.0 )*u_subs )**(1.5 - float(n_ind))
    const3 = np.log( g_function(u_subs, phi, a_subs, b_subs) )
    return const1*const2*const3
################################################################################
def b_function(phi, a_subs, b_subs, n_ind):
    return integrate.fixed_quad(b_argument, 0.0, 1.0, \
                          args=(phi, a_subs, b_subs, n_ind), n = 50)[0]
################################################################################
def d_argument(u_subs, phi, a_subs, b_subs, n_ind):
    const1 = 1.0/( 1.0 - (1.0 - b_subs**2.0 )*u_subs )**(0.5 + float(n_ind))
    const2 = 1.0/( 1.0 - (1.0 - a_subs**2.0 )*u_subs )**(2.5 - float(n_ind))
    const3 = u_subs * g_function(u_subs, phi, a_subs, b_subs)**(-2.0)
    return const1*const2*const3
################################################################################
def d_function(phi, a_subs, b_subs, n_ind):
    return integrate.fixed_quad(d_argument, 0.0, 1.0, \
                          args=(phi, a_subs, b_subs, n_ind), n = 50)[0]
################################################################################
def j0_function(radial, phi, a_subs, b_subs, kappa_s):
    const1 = kappa_radial(radial, kappa_s)*a0_function(a_subs, b_subs)
    const2 = 2.0*kappa_s*b_function(phi, a_subs, b_subs, 0.0)
    return const1 - const2
################################################################################
def p_function(phi, a_subs, b_subs, kappa_s):
    const1 = 2.0*a_subs*b_subs*kappa_s
    const2 = b_function(phi, a_subs, b_subs, 0.0) + \
             math.cos(phi)**2.0 * d_function(phi, a_subs, b_subs, 0.0)
    return 1.0 + const1*const2
################################################################################
def q_function(phi, a_subs, b_subs, kappa_s):
    const1 = 2.0*a_subs*b_subs*kappa_s
    const2 = b_function(phi, a_subs, b_subs, 1.0) + \
             math.sin(phi)**2.0 * d_function(phi, a_subs, b_subs, 2.0)
    return 1.0 + const1*const2
################################################################################
def kappa_mag_zero_minus(phi, a_subs, b_subs, kappa_s):
    a1_eval = a1_function(a_subs, b_subs)
    a0_eval = a0_function(a_subs, b_subs)
    p_eval = p_function(phi, a_subs, b_subs, kappa_s)
    q_eval = q_function(phi, a_subs, b_subs, kappa_s)
    pot_xy_eval = pot_xy(phi, a_subs, b_subs, kappa_s)

    const1 = a1_eval*p_eval + a0_eval*q_eval
    const2 = (a1_eval*p_eval - a0_eval*q_eval)**2.0
    const3 = 4.0*a0_eval*a1_eval*pot_xy_eval*pot_xy_eval
    const4 = 2.0*a_subs*b_subs*a0_eval*a1_eval
    return ( const1 - math.sqrt( const2 + const3 ) )/const4
################################################################################
def kappa_mag_zero_plus(phi, a_subs, b_subs, kappa_s):
    a1_eval = a1_function(a_subs, b_subs)
    a0_eval = a0_function(a_subs, b_subs)
    p_eval = p_function(phi, a_subs, b_subs, kappa_s)
    q_eval = q_function(phi, a_subs, b_subs, kappa_s)
    pot_xy_eval = pot_xy(phi, a_subs, b_subs, kappa_s)

    const1 = a1_eval*p_eval + a0_eval*q_eval
    const2 = (a1_eval*p_eval - a0_eval*q_eval)**2.0
    const3 = 4.0*a0_eval*a1_eval*pot_xy_eval*pot_xy_eval
    const4 = 2.0*a_subs*b_subs*a0_eval*a1_eval
    return ( const1 + math.sqrt( const2 + const3 ) )/const4
################################################################################
def critica_curve_tangential(phi, a_subs, b_subs, kappa_s):
    k_m = kappa_mag_zero_minus(phi, a_subs, b_subs, kappa_s)
    return 2.0*math.exp( -k_m/2.0/kappa_s - 1.0 )
################################################################################
def critica_curve_radial(phi, a_subs, b_subs, kappa_s):
    k_m = kappa_mag_zero_plus(phi, a_subs, b_subs, kappa_s)
    return 2.0*math.exp( -k_m/2.0/kappa_s - 1.0 )
################################################################################
def pot_xx(radial, phi, a_subs, b_subs, kappa_s):
    p_eval = p_function(phi, a_subs, b_subs, kappa_s)
    a0_eval = a0_function(a_subs, b_subs)
    return 1.0 - p_eval + a_subs*b_subs*kappa_radial(radial, kappa_s)*a0_eval
################################################################################
def pot_yy(radial, phi, a_subs, b_subs, kappa_s):
    q_eval = q_function(phi, a_subs, b_subs, kappa_s)
    a1_eval = a1_function(a_subs, b_subs)
    return 1.0 - q_eval + a_subs*b_subs*kappa_radial(radial, kappa_s)*a1_eval
################################################################################
def pot_xy(phi, a_subs, b_subs, kappa_s):
    d_eval = d_function( phi, a_subs, b_subs, 1.0)
    return -a_subs*b_subs*kappa_s*math.sin(2.0*phi)*d_eval
################################################################################
def zeta(a_subs, b_subs):
    a0_eval = a0_function(a_subs, b_subs)
    a1_eval = a1_function(a_subs, b_subs)
    c1 = a_subs*b_subs
    return c1*(a0_eval + a1_eval), c1*(a0_eval - a1_eval)
################################################################################
def xi(phi, a_subs, b_subs, kappa_s):
    b0_eval = b_function(phi, a_subs, b_subs, 0.0)
    b1_eval = b_function(phi, a_subs, b_subs, 1.0)
    d0_eval = d_function(phi, a_subs, b_subs, 0.0)
    d2_eval = d_function(phi, a_subs, b_subs, 2.0)

    c1 = 2.0*a_subs*b_subs*kappa_s
    c2 = b0_eval + math.cos(phi)**2.0*d0_eval
    c3 = b1_eval + math.sin(phi)**2.0*d2_eval
    return c1*(c2 + c3), c1*(c2 - c3)
################################################################################
def q_lambda(r_lambda):
    return (r_lambda + 1.0)/(r_lambda - 1.0)
################################################################################
def abc_tilda(phi, a_subs, b_subs, kappa_s, r_lambda):
    zeta_plus, zeta_minus = zeta(a_subs, b_subs)
    xi_plus, xi_minus = xi(phi, a_subs, b_subs, kappa_s)
    q_eval = q_lambda(r_lambda)

    a_tilda = zeta_plus**2.0 - q_eval**2.0 * zeta_minus**2.0
    c_1 = -4.0*zeta_plus
    c_2 = -2.0*zeta_plus*xi_plus
    c_3 =  2.0*(q_eval**2.0)*zeta_minus*xi_minus
    b_tilda = c_1 + c_2 + c_3

    c_4 = 4.0*(1.0 + xi_plus)
    c_5 = xi_plus**2.0
    c_6 = -(q_eval**2.0)*(xi_minus**2.0)
    c_7 = -4.0*(q_eval**2.0)*pot_xy(phi, a_subs, b_subs, kappa_s)**2.0
    c_tilda = c_4 + c_5 + c_6 + c_7

    return a_tilda, b_tilda, c_tilda
################################################################################
def kappa_dist_const(phi, a_subs, b_subs, kappa_s, r_lambda):
    a_t, b_t, c_t = abc_tilda(phi, a_subs, b_subs, kappa_s, r_lambda)

    c_1 = math.sqrt(b_t**2.0 - 4.0*a_t*c_t)
    return (-b_t + c_1)/(2.0*a_t), (-b_t - c_1)/(2.0*a_t)
################################################################################
def const_dist_curve(phi, a_subs, b_subs, kappa_s, r_lambda):
    kd_plus, kd_minus = np.fabs(kappa_dist_const(phi, a_subs, b_subs, kappa_s, r_lambda))
   # print " "
   # print kd_plus, kd_minus, phi, kappa_s, r_lambda
   # print " "
    c_1 = 2.0*math.exp( -kd_plus/2.0/kappa_s - 1.0 )
    c_2 = 2.0*math.exp( -kd_minus/2.0/kappa_s - 1.0 )
   # print c_1, c_2
    return c_1, c_2
################################################################################
def arg_sigma_radial_num(radial, phi, a_subs, b_subs, kappa_s):
    a0_eval = a0_function(a_subs, b_subs)
    a1_eval = a1_function(a_subs, b_subs)

    p_eval = p_function(phi, a_subs, b_subs, kappa_s)
    q_eval = q_function(phi, a_subs, b_subs, kappa_s)
    potxy_eval = pot_xy(phi, a_subs, b_subs, kappa_s)

    alpha = a_subs**2.0 * b_subs**2.0 * a0_eval * a1_eval
    betta = a_subs * b_subs * (p_eval*a1_eval + q_eval*a0_eval)
    gamma = p_eval*q_eval - potxy_eval**2.0
    kappa_eval = kappa_radial(radial, kappa_s)

    return math.fabs(alpha*kappa_eval**2.0 - betta*kappa_eval + gamma)*radial
################################################################################
def sigma_radial_num(phi, a_subs, b_subs, kappa_s, r_lambda):
    r_min = const_dist_curve(phi, a_subs, b_subs, kappa_s, -r_lambda)[1]
    r_max = const_dist_curve(phi, a_subs, b_subs, kappa_s, r_lambda)[1]
    r_cc = critica_curve_tangential(phi, a_subs, b_subs, kappa_s)
    out = integrate.quad(arg_sigma_radial_num, a = r_min, b = r_max, \
                          args = (phi, a_subs, b_subs, kappa_s), \
                          points = [r_cc])
    return  out[0]
################################################################################
def sigma_num(a_subs, b_subs, kappa_s, r_lambda):
    out = integrate.quad(sigma_radial_num, a = 0.0, b = cn.pi/2.0, \
                          args = (a_subs, b_subs, kappa_s, r_lambda))
    return out[0]*4.0
################################################################################
def anl_integrals(radial, kappa_s):
    #out1 = int r*kappa; out2 = int r*kappa**2
    out1 = -0.5*kappa_s*(radial**2.0)*(1.0 + 2.0*math.log(radial/2.0))
    out2 = (kappa_s**2.0)*(radial**2.0)*\
           (1.0 + 2.0*math.log(radial/2.0)*(1.0 + math.log(radial/2.0)))
    return out1, out2
################################################################################
def det_int(radial, phi, a_subs, b_subs, kappa_s):
    int_1, int_2 = anl_integrals(radial, kappa_s)

    a0_eval = a0_function(a_subs, b_subs)
    a1_eval = a1_function(a_subs, b_subs)

    p_eval = p_function(phi, a_subs, b_subs, kappa_s)
    q_eval = q_function(phi, a_subs, b_subs, kappa_s)
    potxy_eval = pot_xy(phi, a_subs, b_subs, kappa_s)

    alpha = a_subs**2.0 * b_subs**2.0 * a0_eval * a1_eval
    betta = a_subs * b_subs * (p_eval*a1_eval + q_eval*a0_eval)
    gamma = p_eval*q_eval - potxy_eval**2.0
    return alpha*int_2 - betta*int_1 + gamma*radial**2.0/2.0
################################################################################
def sigma_radial(phi, a_subs, b_subs, kappa_s, r_lambda):
    r_min = const_dist_curve(phi, a_subs, b_subs, kappa_s, -r_lambda)[1]
    r_max = const_dist_curve(phi, a_subs, b_subs, kappa_s, r_lambda)[1]
    r_cc = critica_curve_tangential(phi, a_subs, b_subs, kappa_s)

    c_min = det_int(r_min, phi, a_subs, b_subs, kappa_s)
    c_max = det_int(r_max, phi, a_subs, b_subs, kappa_s)
    c_cc  = det_int(r_cc , phi, a_subs, b_subs, kappa_s)
    return c_max + c_min - 2.0*c_cc
################################################################################
def sigma(a_subs, b_subs, kappa_s, r_lambda):
    return integrate.quad(sigma_radial, a = 0.0, b = cn.pi/2.0, \
                          args = (a_subs, b_subs, kappa_s, r_lambda))[0]*4.0
################################################################################
if __name__ == '__main__':

    #print a1_function(1.2, 1.2)
    #print b_argument(0.5, 1.0, 0.5, 1.5, 1.0)
    #print b_function(2.0, 0.5, 1.5, 1.0)
    #print d_function(2.0, 0.5, 1.5, 1.0)
    #print kappa_mag_zero_minus(2.0, 0.5, 1.5, 1.0)
    #print kappa_mag_zero_plus(2.0, 0.5, 1.5, 1.0)
    print critica_curve_tangential(0.00000001, 1.1, 1.0, 1.0)
