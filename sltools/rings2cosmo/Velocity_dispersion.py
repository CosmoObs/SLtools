import numpy as np
import scipy as sp
import scipy.integrate as it

from astropy import constants as const
from astropy import units as u
from astropy.cosmology import FlatLambdaCDM

# Cosmology used:
cosmo = FlatLambdaCDM(H0 = 67.3, Om0 = 0.315)

# Physical constants:
c = (const.c).to(u.km/u.second) #light velocity in Km/s This is bc velocity dispersion is measured in (km/s)^2
clight = c.value

def Lambda_function(x):
    """Eq. (15) from arXiv:0907.4992v2

    Args:
        x (float): parameter

    Returns:
        float: ratio of gamma funcions.
    """
    return sp.special.gamma((x - 1) / 2) / sp.special.gamma(x / 2)

### brightness profile
def brightness_3D (delta, r):
    return r**(-delta)




##Radial Sigma for constant Beta Eq. (16) from arXiv:0907.4992v2
def radial_sigma_cte_beta (r, z_S, z_L, theta_E, alpha, delta, beta, gamma):

    #distances
    DS = cosmo.angular_diameter_distance(z_S).value
    DL = cosmo.angular_diameter_distance(z_L).value
    DLS = cosmo.angular_diameter_distance_z1z2(z_L, z_S).value
    #phisical einstein radius
    RE = theta_E*DL
    #relation alpha-delta-xi
    xi = alpha + delta - 2
    #####constants####
    ## relativistical term
    term_RE  = (1 / (1 + gamma)) * ((clight ** 2) / 2) * (DS / DLS) * theta_E
    #dynamics
    term_dynamic = 2/(((np.sqrt(np.pi)))*(Lambda_function(alpha))*(xi - (2*beta)))
    ####term integrable####
    integral_term = r/RE 
    
    return term_RE*term_dynamic*(integral_term**(2 - alpha))
###########################################################################################################################################################################


###### mean sigma LOS #####
#weignt function
def SDSS_seeing_weight (R, z_L, seeing_atm, theta_ap):

    #distances
    DL = cosmo.angular_diameter_distance(z_L).value
    #seeing information
    chi = theta_ap/seeing_atm
    tilde_sigma = seeing_atm * \
        np.sqrt(1 + (chi ** 2) / 4 + (chi ** 4) / 40)  # Eq. (20)
    
    r_angular = R/DL
    x = r_angular/tilde_sigma
    return np.exp(-(x**2)/2)

## the normalization term or in another words the denominator part of sigma star
def Weight_normalization(seeing_atm, theta_ap,z_L, theta_E, delta):
    ## \tilde{\sigma}
    chi = theta_ap/seeing_atm

    tilde_sigma = seeing_atm * \
        np.sqrt(1 + (chi ** 2) / 4 + (chi ** 4) / 40)  # Eq. (20)
    
    #R_E
    DL = cosmo.angular_diameter_distance(z_L).value
    RE = DL*theta_E

    #term B\left ( \frac{\delta - 1}{2}, \frac{1}{2} \right )
    #x = (delta - 1)/2
    euler_beta = Lambda_function(delta)*np.sqrt(np.pi)

    #term  \Gamma\left ( \frac{3 - \delta}{2} \right )
    y = (3 - delta)/2
    euler_gamma = sp.special.gamma(y)

    #term \left(\frac{2\tilde{\sigma}^2_{atm}}{\theta_E^2}\right)
    tilde_sigma_by_theta_E = 2*((tilde_sigma/theta_E)**2)

    return euler_beta*(RE**(2*y))*(tilde_sigma_by_theta_E**y)*euler_gamma/2



#this is the function that I have to integrate, (the up part of sigma star before the integration) I do intend to writhe this in a better way, but this is what I have for today
def integrando (r, R, z_L,z_S,theta_E, seeing_atm, theta_ap, alpha, beta, delta,gamma):
    
    #integral terms
    anisotropy = 1 - (((R/r)**2)*beta)
    

    omega = SDSS_seeing_weight (R, z_L, seeing_atm, theta_ap)

    nu = r**(-delta)
    sigma = radial_sigma_cte_beta (r, z_S, z_L, theta_E, alpha, delta, beta, gamma)

    projection = (r**2) - (R**2)

    spherical_coordinates = r*R

    return 2*nu*sigma*omega*anisotropy*spherical_coordinates/(np.sqrt(projection))

#here I ain going to defnetely perform the sigma star calculation
def sigma_star(z_L,z_S,theta_E, seeing_atm, theta_ap, alpha, beta, delta,gamma):
    integral = it.dblquad(integrando, 0, np.inf, lambda r: r, lambda r: np.inf, \
                          args=(z_L,z_S, theta_E,seeing_atm, theta_ap, alpha, beta, delta,gamma), epsabs=1.49e-03)
    sigma_star = integral[0]/Weight_normalization(seeing_atm, theta_ap,z_L, theta_E, delta)
    return np.sqrt(sigma_star)

    