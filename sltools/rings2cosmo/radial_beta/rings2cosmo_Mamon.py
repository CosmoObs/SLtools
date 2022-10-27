import emcee
import numpy as np
import scipy as sp

from scipy.optimize import minimize
from scipy import special
from astropy import constants as const
from astropy import units as u
from astropy.cosmology import FlatLambdaCDM
from multiprocessing import Pool

# Cosmology used:
# From Wikipedia (https://en.wikipedia.org/wiki/Lambda-CDM_model)
Hubble = 67.74
cosmo = FlatLambdaCDM(H0 = Hubble, Om0=0.3089)

# Physical constants:
c = (const.c).to(u.km/u.second)
clight = c.value


###########################################################################################################################
################### Functions used in the model of Sigma #######################
###########################################################################################################################
def Lambda_function(x):
    """Eq. (15) from arXiv:0907.4992v2

    Args:
        x (float): parameter

    Returns:
        float: ratio of gamma funcions.
    """
    return special.gamma((x - 1) / 2) / special.gamma(x / 2)

#reduced mass within the virial radius
def virial_mass (DS, DL, DLS, theta_E, gamma, alpha):
    gamma_term = 1.0/(1.0 + gamma)
    dinamics = Lambda_function(alpha)*(np.sqrt(np.pi))
    distance_ratio = DS/DLS
    Ering_ratio = (DL**(alpha - 2))/(theta_E**(1 - alpha))
    
    return (clight**2)*gamma_term*distance_ratio*Ering_ratio/dinamics

#virial radius usin r200
def virial_radii (DS, DL, DLS, theta_E, gamma, alpha):
    MuE = virial_mass (DS, DL, DLS, theta_E, gamma, alpha)
    D = 200
    return (2*MuE/(D*Hubble))**(1/alpha)

#b=DL*tilde_sigma dispersion used in the wight function
def b_func (DL, seeing_atm, theta_ap):
    chi = theta_ap/seeing_atm
    tilde_sigma = seeing_atm *np.sqrt(1.0 + (chi ** 2.0) / 4.0 + (chi ** 4.0) / 40.0)  # Eq. (20)
    
    return tilde_sigma*DL
#normalization term
def Weight_normalization(seeing_atm, theta_ap,DL, theta_E, delta):
    ## \tilde{\sigma}
    chi = theta_ap/seeing_atm

    tilde_sigma = seeing_atm * \
        np.sqrt(1 + (chi ** 2) / 4 + (chi ** 4) / 40)  # Eq. (20)
    
    RE = DL*theta_E

    #term B\left ( \frac{\delta - 1}{2}, \frac{1}{2} \right )
    #x = (delta - 1)/2
    euler_beta = Lambda_function(delta)*np.sqrt(np.pi)

    #term  \Gamma\left ( \frac{3 - \delta}{2} \right )
    y = (3 - delta)/2
    euler_gamma = special.gamma(y)

    #term \left(\frac{2\tilde{\sigma}^2_{atm}}{\theta_E^2}\right)
    tilde_sigma_by_theta_E = 2*((tilde_sigma/theta_E)**2)

    return euler_beta*(RE**(2*y))*(tilde_sigma_by_theta_E**y)*euler_gamma/2

#Einstein ring dependence function
def Relativistic_term (DS, DLS, theta_E, gamma):
    return (clight**2.0)*(DS/DLS)*theta_E/(2.0*(1.0+gamma))

#parameter from Mamon's model
def mamon_parameters (DL, DS, DLS, theta_E, alpha, delta,gamma):
    xi = delta + alpha - 2
    R_E = theta_E*DL
    
    rel_term = Relativistic_term (DS, DLS, theta_E, gamma)
    dim = 2.0/(np.sqrt(np.pi)*Lambda_function(alpha)*xi*(xi - 1)*(R_E**(2 - alpha)))
    
    return rel_term*dim
    

###########################################################################################################################
################### Velocity Dispersion #######################
###########################################################################################################################
import VDmod

def VD_scalar(z_S, z_L, theta_E, seeing_atm, theta_ap, virial_frac, alpha, delta, gamma):
    DS = cosmo.angular_diameter_distance(z_S).value
    DL = cosmo.angular_diameter_distance(z_L).value
    DLS = cosmo.angular_diameter_distance_z1z2(z_L, z_S).value
    xi = delta + alpha - 2
    
    b = b_func (DL, seeing_atm, theta_ap)
    ra = virial_frac*virial_radii (DS, DL, DLS, theta_E, gamma, alpha)
    integral = VDmod.integralVD(b,alpha,delta,ra)

    
    C_mamon_term = mamon_parameters (DL, DS, DLS, theta_E, alpha, delta,gamma)
    weigth_term = Weight_normalization(seeing_atm, theta_ap,DL, theta_E, delta)
    vel = C_mamon_term*integral/weigth_term
    
    return np.sqrt(np.abs(vel))


vel = np.vectorize(VD_scalar)

###########################################################################################################################

# Goodness of fit of a statistical model
def log_likelihood(theta, z_S, z_L, velDisp, velDispErr, theta_E, seeing_atm, theta_ap, virial_frac):
    """log(Eq. (25)) from arXiv:0907.4992v2 

    Args:
        theta (list): list of parameters [alpha, delta, gamma]
        z_S (float): source redshift
        z_L (float): lens redshift
        velDisp (float): velocity dispersion
        velDispErr (float): velocity dispersion error (std. dev.)
        theta_E (float): Einstein radius (in radians)
        seeing_atm (float): Atmospheric seeing (in radians)
        theta_ap (float): Aperture size (in radians)

    Returns:
        float: Log-likehood function given model for velocity dispersion and the one measured.
    """
    alpha, delta, gamma = theta
    model = vel(z_S, z_L, theta_E, seeing_atm,
                theta_ap, virial_frac, alpha, delta, gamma)
    return - 0.5*np.sum((velDisp - model) ** 2 / (velDispErr ** 2) + np.log(2 * np.pi * velDispErr ** 2))



def log_prior(theta, alpha_0, eps_alpha_0, delta_0, eps_delta_0):
    """Gaussian priors

    Args:
        theta (list): list of parameters [alpha, delta, gamma]
        alpha_0 (float): expected value for alpha
        eps_alpha_0 (float): variance of alpha
         
        delta_0 (float): expected value for delta
        eps_delta_0 (float): variance of delta

    Returns:
        float: Sum of log of priors for alpha, and delta.
    """
    alpha, delta, gamma = theta
    n_sigma = 5
    if (alpha_0[0] - n_sigma * eps_alpha_0[0] < alpha < alpha_0[0] + n_sigma * eps_alpha_0[0]) and \
        (delta_0[0] - n_sigma * eps_delta_0[0] < delta < delta_0[0] + n_sigma * eps_delta_0[0]):
        log_prior_alpha = - 0.5 * \
            np.sum((alpha - alpha_0)**2 / eps_alpha_0 **
                   2 + np.log(2 * np.pi * eps_alpha_0**2))
        log_prior_delta = - 0.5 * \
            np.sum((delta - delta_0) ** 2 / eps_delta_0 **
                   2 + np.log(2 * np.pi * eps_delta_0**2))
        return log_prior_alpha + log_prior_delta
    else:
        return - np.inf

def log_probability(theta, z_S, z_L, velDisp, velDispErr, theta_E, seeing_atm, theta_ap, virial_frac,
                    alpha_0, eps_alpha_0, delta_0, eps_delta_0):
    """Log of probability of interest

    Args:
        theta (list): list of parameters [alpha, delta, gamma]
        z_S (float): source redshift
        z_L (float): lens redshift
        velDisp (float): velocity dispersion
        velDispErr (float): velocity dispersion error (std. dev.)
        theta_E (float): Einstein radius (in radians)
        seeing_atm (float): Atmospheric seeing (in radians)
        theta_ap (float): Aperture size (in radians)
        alpha_0 (float): expected value for alpha
        eps_alpha_0 (float): variance of alpha
        delta_0 (float): expected value for delta
        eps_delta_0 (float): variance of delta

    Returns:
        float: Eq. (27) from arXiv:0907.4992v2
    """
    lp = log_prior(theta, alpha_0, eps_alpha_0, delta_0, eps_delta_0)
    if not np.isfinite(lp):
        return - np.inf
    else:
        return lp + log_likelihood(theta, z_S, z_L, velDisp, velDispErr, theta_E, seeing_atm, theta_ap, virial_frac)




# Minimizations and sampling methods


def minimization_loglikelihood(z_S, z_L, velDisp, velDispErr, theta_E, seeing_atm, theta_ap, virial_frac,
                               seed=42, alpha_ini=2.0, delta_ini=2.4, gamma_ini=1.0):
    """Maximization of Likehood function

    Args:
        z_S (float): source redshift
        z_L (float): lens redshift
        velDisp (float): velocity dispersion
        velDispErr (float): velocity dispersion error (std. dev.)
        theta_E (float): Einstein radius (in radians)
        seeing_atm (float): Atmospheric seeing (in radians)
        theta_ap (float): Aperture size (in radians)
        seed (float): random seed for reproducibility purposes. Default: 42.
        alpha_ini (float, optional): Initial guess for alpha. Default: 2.0.
        delta_ini (float, optional): Initial guess for delta. Default: 2.4.
        gamma_ini (float, optional): Initial guess for gamma. Default: 1.0.

    Returns:
        list: list of alpha, delta, and gamma obtained from minimization of likelihood function.
    """
    np.random.seed(seed)
    nll = lambda *args: - log_likelihood(*args)

    initial = np.array([alpha_ini,delta_ini, gamma_ini]) + \
        1e-5 * np.random.randn(3)

    soln = minimize(nll, initial, args=(z_S, z_L, velDisp,
                    velDispErr, theta_E, seeing_atm, theta_ap, virial_frac))  # , method='Nelder-Mead', tol=1e-10)

    alpha_ml, delta_ml, gamma_ml = soln.x
    return alpha_ml, delta_ml, gamma_ml

def minimization_logprobability(z_S, z_L, velDisp, velDispErr, theta_E, seeing_atm, theta_ap, virial_frac,
                                seed=42, alpha_ini=2.0, delta_ini=2.4, gamma_ini=1.0,
                                alpha_0_value=2.0, eps_alpha_0_value=0.08,
                                delta_0_value=2.4, eps_delta_0_value=0.11):
    """Maximization of Log Probability function 

    Args:
        z_S (float): source redshift
        z_L (float): lens redshift
        velDisp (float): velocity dispersion
        velDispErr (float): velocity dispersion error (std. dev.)
        theta_E (float): Einstein radius (in radians)
        seeing_atm (float): Atmospheric seeing (in radians)
        theta_ap (float): Aperture size (in radians)
        seed (float): random seed for reproducibility purposes.
        alpha_ini (float, optional): Initial guess for alpha. Defaults to 2.0.
        delta_ini (float, optional): Initial guess for delta. Defaults to 2.4.
        gamma_ini (float, optional): Initial guess for gamma. Defaults to 1.0.
        alpha_0_value (float, optional): expected value for alpha. Defaults to 2.0.
        eps_alpha_0_value (float, optional): variance for alpha. Defaults to 0.08.
        delta_0_value (float, optional): expected value for delta. Defaults to 2.4.
        eps_delta_0_value (float, optional): variance for delta. Defaults to 0.11.

    Returns:
        list: list of alpha, delta, and gamma obtained from minimization of log-probability function.
    """
    alpha_0 = np.repeat(alpha_0_value, len(z_S))
    eps_alpha_0 = np.repeat(eps_alpha_0_value, len(z_S))

    delta_0 = np.repeat(delta_0_value, len(z_S))
    eps_delta_0 = np.repeat(eps_delta_0_value, len(z_S))

    np.random.seed(seed)
    nll_2 = lambda *args: - log_probability(*args)
    initial = np.array([alpha_ini, delta_ini, gamma_ini]) + \
        1e-5 * np.random.randn(3)

    #soln_2 = minimize(nll_2, initial, args=(z_S, z_L, velDisp, velDispErr, theta_E,
    #                  seeing_atm, theta_ap, virial_frac, alpha_0, eps_alpha_0, delta_0, eps_delta_0, ))  # , method='Nelder-Mead', tol=1e-10)
    ####################################################
    #this contrain is to enshure that none parameter is <0, cus the log function does not handle negative values
    abstol = 1e-30
    cons = ({'type': 'ineq', 'fun': lambda x: x - abstol})
    soln_2 = minimize(nll_2, initial, args=(z_S, z_L, velDisp, velDispErr, theta_E,
                    seeing_atm, theta_ap, virial_frac, alpha_0, eps_alpha_0, delta_0, eps_delta_0, ), constraints = cons)  # , method='Nelder-Mead', tol=1e-10)
    
    alpha_ml2, delta_ml2, gamma_ml2 = soln_2.x
    

    return float(alpha_ml2), float(delta_ml2), float(gamma_ml2)



def logprobability_sampling(z_S, z_L, velDisp, velDispErr, theta_E, seeing_atm, theta_ap, virial_frac,
                            seed=42, alpha_ini=2.0, delta_ini=2.4, gamma_ini=1.0,
                            alpha_0_value=2.0, eps_alpha_0_value=0.08,
                            delta_0_value=2.4, eps_delta_0_value=0.11,
                            n_dim=3, n_walkers=64, n_burn=500, n_steps=10000, progress=True, processes=1):
    """Sampling logprobability function with emcee

    Args:
        z_S (float): source redshift
        z_L (float): lens redshift
        velDisp (float): velocity dispersion
        velDispErr (float): velocity dispersion error (std. dev.)
        theta_E (float): Einstein radius (in radians)
        seeing_atm (float): Atmospheric seeing (in radians)
        theta_ap (float): Aperture size (in radians)
        seed (float): random seed for reproducibility purposes.
        alpha_ini (float, optional): Initial guess for alpha. Defaults to 2.0.
        delta_ini (float, optional): Initial guess for delta. Defaults to 2.4.
        gamma_ini (float, optional): Initial guess for gamma. Defaults to 1.0.
        alpha_0_value (float, optional): expected value for alpha. Defaults to 2.0.
        eps_alpha_0_value (float, optional): variance for alpha. Defaults to 0.08.
        delta_0_value (float, optional): expected value for delta. Defaults to 2.4.
        eps_delta_0_value (float, optional): variance for delta. Defaults to 0.11.
        n_dim (int, optional): number of parameters in the model (r and p). Defaults to 4.
        n_walkers (int, optional): number of MCMC walkers. Defaults to 64.
        n_burn (int, optional): "burn-in" period to let chains stabilize. Defaults to 500.
        n_steps (int, optional): number of MCMC steps to take after burn-in. Defaults to 10000.
        progress (bool, optional): Show progress bar. Defaults to True.
        processes (int, optional): Number of processes in parallel. Defaults to 1.
    """
    alpha_0 = np.repeat(alpha_0_value, len(z_S))
    eps_alpha_0 = np.repeat(eps_alpha_0_value, len(z_S))

    delta_0 = np.repeat(delta_0_value, len(z_S))
    eps_delta_0 = np.repeat(eps_delta_0_value, len(z_S))

    with Pool(processes=processes) as pool:

        sampler = emcee.EnsembleSampler(n_walkers, n_dim, log_probability, args=(z_S, z_L, velDisp, velDispErr, theta_E,
                                                                                 seeing_atm, theta_ap, virial_frac, alpha_0,
                                                                                 eps_alpha_0, delta_0, eps_delta_0, ), pool=pool)
        np.random.seed(seed)
        solu = minimization_loglikelihood(z_S, z_L, velDisp, velDispErr, theta_E, seeing_atm, theta_ap, virial_frac,
                                          seed, alpha_ini, delta_ini, gamma_ini)
        p0 = solu + 1e-3 * np.random.randn(n_walkers, n_dim)

        # Run n_burn steps as a burn-in:
        print('Running burn-in ...')
        pos, prob, state = sampler.run_mcmc(p0, n_burn, progress=progress)

        # Reset the chain to remove the burn-in samples:
        sampler.reset()

        # Starting from the final position in the burn-in chain, sample for n_steps steps:
        print('Sampling ...')
        sampler.run_mcmc(pos, n_steps, rstate0=state, progress=progress)

    return sampler
##########################################################################################################################################


