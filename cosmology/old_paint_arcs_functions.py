from scipy.integrate import dblquad, quad
from scipy import constants
import math
import re
import string
import numpy


def proper_dist_int(z,omega_m,omega_lambda):
	'''Computes the integrand for the proper distance integral for a LCDM model with omega_k = 0.
	
	Input:
	- z <float>: redshift
	- h_0 <float>: Hubble constant / 100.
	- omega_m <float>: matter density (relative to the critical density)
	- omega_lambda <float>: density associated to the cosmological constant (relative to the critical density)
	
	Output:
	- <float>: integrand for the proper distance integral
	'''
	
	return 1./math.sqrt(omega_m*((1.+z)**3)+omega_lambda)


def proper_dist(z1,z2,h_0,omega_m,omega_lambda):

	'''Computes the angular distance for a LCDM model with omega_k = 0
	
	Input:
	- z1 <float>: observer redshift
	- z2 <float>: objtect redshift
	- h_0 <float>: Hubble constant / 100.
	- omega_m <float>: matter density (relative to the critical density)
	- omega_lambda <float>: density associated to the cosmological constant (relative to the critical density)
	
	Output:
	- <float>: proper distance in Mpc
	'''	

	c_h_0 = 3000./h_0 #(c_h_0 in Mpc)


	return c_h_0 * quad(lambda z: proper_dist_int(z,omega_m,omega_lambda),z1,z2)[0]

def angular_diameter_distance(z1,z2,h_0,omega_m,omega_lambda):
	'''Computes the angular diameter distance in Mpc using a LCDM model with omega_k = 0.
	
	Input:
	- z1 <float>: observer redshift
	- z2 <float>: objtect redshift
	- h_0 <float>: Hubble constant / 100.
	- omega_m <float>: matter density (relative to the critical density)
	- omega_lambda <float>: density associated to the cosmological constant (relative to the critical density)
	
	Output:
	- <float>: angular diameter distance in Mpc
	'''	

	return proper_dist(z1,z2,h_0,omega_m,omega_lambda)/ (1.+z2)
