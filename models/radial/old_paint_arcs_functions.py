from scipy.integrate import dblquad, quad
from scipy import constants
import math
import re
import string
import numpy


def get_sersic_b_n(m):

	''' Finds the b_n coefficient of the Sersic profile.
	This expression is more accurate thant the usual b_n = 1.9992*n - 0.3271
	Reference: Ciotti, L., Bertin, G. 1999, A&A, 352, 447
	MacArthur, L.A., Courteau, S., & Holtzman, J.A. 2003, ApJ, 582, 689
	b_n is chosen so that half the total luminosity predicted by the Sersic profile comes from r <= r_e .

	@param m Sersic index 

	@return b_n coefficient of the Sersic profile 	
	'''

	if m == 0.0: return -0.3271
	b_n = 2.0*m - 1.0/3.0 + 4.0/m/405.0 + 46.0/m/m/25515.0 \
		+ 131.0/m/m/m/1148175.0  \
                - 2194697.0/m/m/m/m/30690717750.0

	return b_n


def einstein_radius_sis(m200,z1,z2,h_0,omega_m,omega_lambda):
	'''Computes the Einstein radius in arcseconds for a singular isothermal sphere lens at z=z1
	and source at z=z2, considering a LCDM model with omega_k = 0.
	
	Input:
	- m200 <float>: mass within r200 (radius that encloses 200 times the critical density)
	- z1 <float>: lens redshift
	- z2 <float>: source redshift
	- h_0 <float>: Hubble constant / 100.
	- omega_m <float>: matter density (relative to the critical density)
	- omega_lambda <float>: density associated to the cosmological constant (relative to the critical density)
	
	Output:
	- <float>: Einstein radius in arcseconds
	'''	

	sigma = velocity_dispersion(m200,z1,h_0,omega_m,omega_lambda)
	#sigma must be in km/s

	dfrac = angular_diameter_distance(z1,z2,h_0,omega_m,omega_lambda)/angular_diameter_distance(0,z2,h_0,omega_m,omega_lambda)

	return (4.*math.pi*dfrac*(sigma/300000)**2)*206264.81
	

def velocity_dispersion(m200,z,h_0,omega_m,omega_lambda):

	'''Computes the velocity dispersion in km/s for a singular isothermal sphere at redshift z
	using a LCDM model with omega_k = 0.
	
	Input:
	- m200 <float>: mass within r200 (radius that encloses 200 times the critical density)
	- z <float>: redshift
	- h_0 <float>: Hubble constant / 100.
	- omega_m <float>: matter density (relative to the critical density)
	- omega_lambda <float>: density associated to the cosmological constant (relative to the critical density)
	
	Output:
	- <float>: velocity dispersion in km/s
	'''	

	e_z = 1./(proper_dist_int(z,omega_m,omega_lambda))**2

	return 475.*(m200*h_0*math.sqrt(200.)*e_z/(10e15))**(1./3.)
