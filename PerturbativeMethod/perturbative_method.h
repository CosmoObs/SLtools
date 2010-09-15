/** @file
* Example of doxygen documentation for C functions. 
*/

/** @package perturbative_method
*  Package to compute quantities related to the Perturbative Method discussed in  http://adsabs.harvard.edu/abs/2007MNRAS.382L..58A and further works 
*
*  Detailed descrition FIXME
*
*/



#ifndef PERTURBATIVEMETHOD_H
#define PERTURBATIVEMETHOD_H


#include <cmath> 

/*********************************************************************************************************************/
/** \brief \f$ \phi(r) \f$ type function 
*
* \f$ r \f$ is the radial coordinate and pot_params[] is a vector that contains all the potential parameters (ex: \f$ \kappa_s \f$ and \f$ \rho_s \f$ from NFW profile)
*
*/

typedef double (*phi_type) (double r, double pot_params[]);
/** \brief \f$ \psi(r,\theta) \f$ type function 
*
* \f$ r \f$ and \f$ \theta \f$ are the radial and angular coordinate, pert_params[] is a vector that contains all the perturbation parameters
*
*/
typedef double (*psi_type) (double r, double theta, double pert_params[]);

/** \brief \f$ f_{n}(\theta)\equiv \frac{1}{n!} \left[\frac{\partial^n\psi(r,\theta)}{\partial r^{n}}\right]_{r=r_e} \f$ type function 
*
* \f$ \theta \f$ is the angular coordinate, pert_params[] is a vector that contains all the perturbation parameters
*
*/
typedef double (*f_type)   (double theta, double pert_params[]);
/*********************************************************************************************************************/


/** \brief A struct to define elliptical sources
* 
*
*  x0 (\f$ x_0 \f$) -> position at x axis
*
*  y0 (\f$ y_0\f$) ->     position at y axis
*
*  R0 (\f$ R_0\f$)->     characteristic size
*
*  eta0 (\f$ \eta_0\f$)->    parameter relater to the ellipticity ( e = sqrt(2*eta0) )
*
*  theta0 (\f$ \theta_0\f$)->  inclination to the main axis
*/
typedef struct {
  double x0 ;     //position at x axis
  double y0 ;     //position at y axis
  double R0 ;     //characteristic size
  double eta0 ;   //parameter relater to the ellipticity
  double theta0 ; //inclination to the main axis
} elliptical_source;



/**  \f$ \overline{f_1}(\theta) = f_1(\theta) + x0\cos(\theta) + y0\sin(\theta) \f$
*
*  \param f1_in function related to the perturbed potential 
*  \param source_in a struct that define elliptical sources
*  \param theta angular coordinate in the lens plane 
*  \param pert_params[] a vector that contains all the perturbation parameters
*  \return \f$ f_1(\theta) + x0\cos(\theta) + y0\sin(\theta) \f$
*
*  \sa f_type, elliptical_source
*/
double f1_bar(f_type f1_in, elliptical_source source_in, double theta, double pert_params[]){
  return f1_in(theta, pot_params) + source_in.x0*cos(theta) + source_in.y0*sin(theta);
}

/**  \f$ \overline{\frac{d f_1(\theta)}{d \theta}} \equiv \frac{d f_0(\theta)}{d \theta} - x0\sin(\theta) + y0\cos(\theta) \f$
*
*  \param Df0Dtheta_in function related to the perturbed potential
*  \param source_in a struct that define elliptical sources
*  \param theta angular coordinate in the lens plane 
*  \param pert_params[] a vector that contains all the perturbation parameters
*  \return \f$ \frac{d f_0(\theta)}{d \theta} - x0\sin(\theta) + y0\cos(\theta) \f$
*
*  \sa f_type, elliptical_source
*/
double Df0Dtheta_bar(f_type Df0Dtheta_in, elliptical_source source_in, double theta, double pot_params[]){
  return Df0Dtheta_in(theta, pot_params) - source_in.x0*sin(theta) +  source_in.y0*cos(theta);
}


/**  Argument of square root: \f$ R_{0}^2(1-\eta_0\cos(\theta - \theta_0))-(1-\eta_0^2) \overline{\frac{d f_1(\theta)}{d \theta}} \f$ 
*
*  \param Df0Dtheta_in function related to the perturbed potential
*  \param source_in a struct that define elliptical sources
*  \param theta angular coordinate in the lens plane 
*  \param pert_params[] a vector that contains all the perturbation parameters
*  \return \f$ R_{0}^2(1-\eta_0\cos(2\theta - 2\theta_0))-(1-\eta_0^2) \overline{\frac{d f_1(\theta)}{d \theta}}^2 \f$
*
*  \sa f_type, elliptical_source
*/
double arg_sqrt(f_type Df0Dtheta_in, elliptical_source source_in, double theta, double pot_params[]){
  double Df0Dtheta_bar_tmp = Df0Dtheta_bar(Df0Dtheta_in, source_in, theta, pot_params);
  double Df0Dtheta_bar_tmp2 = Df0Dtheta_bar_tmp*Df0Dtheta_bar_tmp;
  double theta_bar = theta - source_in.theta0;
  double S = 1.0 - source_in.eta0*cos(2.0*theta_bar);

  double R0_tmp2 = source_in.R0*source_in.R0;
  double eta0_tmp2 = source_in.eta0*source_in.eta0;

  return R0_tmp2*S - (1.0-eta0_tmp2)*Df0Dtheta_bar_tmp2;
}

/**  
*
*  \param Df0Dtheta_in function related to the perturbed potential
*  \param f1_in function related to the perturbed potential
*  \param kappa2 \f$  \kappa_2 = 1 - \left[\frac{d^2 \phi_0(r)}{dr^2}\right]_{r=r_e} \f$
*  \param source_in a struct that define elliptical sources
*  \param theta angular coordinate in the lens plane 
*  \param pert_params[] a vector that contains all the perturbation parameters
*  \return \f$ \frac{1}{\kappa_2}\left\{ \overline{f_1}(\theta) + \sin(2\theta - 2\theta_0)\frac{\eta_0}{1-\eta_0\cos(2\theta - 2\theta_0)}\frac{d\overline{f_0}(\theta)}{d\theta} + \frac{\sqrt{R_{0}^2(1-\eta_0\cos(2\theta - 2\theta_0))-(1-\eta_0^2) \overline{\frac{d f_1(\theta)}{d \theta}}^2}}{1-\eta_0\cos(2\theta - 2\theta_0)} \right\}\f$
*
*  \sa f_type, elliptical_source
*/
double dr_plus(f_type f1_in, f_type Df0Dtheta_in, double kappa2, elliptical_source source_in, double theta, double pot_params[]){
  double Df0Dtheta_bar_tmp = Df0Dtheta_bar(Df0Dtheta_in, source_in, theta, pot_params);
  double Df0Dtheta_bar_tmp2 = Df0Dtheta_bar_tmp*Df0Dtheta_bar_tmp;

  double f1_bar_tmp = f1_bar(f1_in, source_in, theta, pot_params);

  double theta_bar = theta - source_in.theta0;
  double S = 1.0 - source_in.eta0*cos(2.0*theta_bar);

  double R0_tmp2 = source_in.R0*source_in.R0;
  double eta0_tmp2 = source_in.eta0*source_in.eta0;

  double sqrt_arg = R0_tmp2*S - (1.0-eta0_tmp2)*Df0Dtheta_bar_tmp2;

  double quant_tmp = sin(2.0*theta_bar)* (source_in.eta0/S) * Df0Dtheta_bar_tmp;


  return 1.0/kappa2 * (f1_bar_tmp + quant_tmp + sqrt(sqrt_arg)/S);
}

/**  
*
*  \param Df0Dtheta_in function related to the perturbed potential
*  \param f1_in function related to the perturbed potential
*  \param kappa2 \f$  \kappa_2 = 1 - \left[\frac{d^2 \phi_0(r)}{dr^2}\right]_{r=r_e} \f$
*  \param source_in a struct that define elliptical sources
*  \param theta angular coordinate in the lens plane 
*  \param pert_params[] a vector that contains all the perturbation parameters
*  \return \f$ \frac{1}{\kappa_2}\left\{ \overline{f_1}(\theta) + \sin(2\theta - 2\theta_0)\frac{\eta_0}{1-\eta_0\cos(2\theta - 2\theta_0)}\frac{d\overline{f_0}(\theta)}{d\theta} - \frac{\sqrt{R_{0}^2(1-\eta_0\cos(2\theta - 2\theta_0))-(1-\eta_0^2) \overline{\frac{d f_1(\theta)}{d \theta}}^2}}{1-\eta_0\cos(2\theta - 2\theta_0)} \right\}\f$
*
*  \sa f_type, elliptical_source
*/
double dr_minus(f_type f1_in, f_type Df0Dtheta_in, double kappa2, elliptical_source source_in, double theta, double pot_params[]){
  double Df0Dtheta_bar_tmp = Df0Dtheta_bar(Df0Dtheta_in, source_in, theta, pot_params);
  double Df0Dtheta_bar_tmp2 = Df0Dtheta_bar_tmp*Df0Dtheta_bar_tmp;

  double f1_bar_tmp = f1_bar(f1_in, source_in, theta, pot_params);

  double theta_bar = theta - source_in.theta0;
  double S = 1.0 - source_in.eta0*cos(2.0*theta_bar);

  double R0_tmp2 = source_in.R0*source_in.R0;
  double eta0_tmp2 = source_in.eta0*source_in.eta0;

  double sqrt_arg = R0_tmp2*S - (1.0-eta0_tmp2)*Df0Dtheta_bar_tmp2;

  double quant_tmp = sin(2.0*theta_bar)* (source_in.eta0/S) * Df0Dtheta_bar_tmp;


  return 1.0/kappa2 * (f1_bar_tmp + quant_tmp - sqrt(sqrt_arg)/S);
}


#endif
