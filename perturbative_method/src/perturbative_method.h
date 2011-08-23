/** @file
* This module compute or define several quantities related to the Perturbative Approach as
* - structure for the elliptical source
* - \f$ \bar{f}_n(\theta)\f$  and their derivatives
* - \f$ \Delta(\eta_0,\theta) \f$ (argument of the square root, it determines the dark regions) 
* - \f$ dr_{\pm} \f$ (perturbed position)
* - \f$ r_{\rm crit}(\theta)\f$  Radial coordinate of the Tangential Critical Curve
* - \f$ r_{\rm caust}(\theta)\f$ Radial coordinate of the Tangential Caustic
* - Parametric equation of the tangential caustic
* - Parametric Equation of an elliptical source not aligned to main axis
* - Mean value of the radial distance of the arcs
* - Mean width of the image at some angular position
* - Unitary area of the image, i.e, is the argument of the integral \f$ W(\theta)\bar{r}(\theta)\f$ 
*
* All of these expressions are described in .../sltools/PerturbativeMethod/writeups/Report_on_Perturbative_Method.pdf
*/

/** @package perturbative_method
*  To compute quantities related to the Perturbative Method discussed in  http://adsabs.harvard.edu/abs/2007MNRAS.382L..58A and further works 
*
*/


#ifndef PERTURBATIVEMETHOD_H
#define PERTURBATIVEMETHOD_H


#include <math.h> 

/*********************************************************************************************************************/
/** \brief \f$ \phi(r) \f$ type function 
*
* \f$ r \f$ is the radial coordinate and pert_params[] is a vector that contains all the potential parameters (ex: \f$ \kappa_s \f$ and \f$ \rho_s \f$ from NFW profile)
*
*/
typedef double (*phi_type) (double r, double pert_params[]);

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
// typedef double (*f_type)   (double theta, double pert_params[]);
typedef double (*f_type)   (double theta, double pert_params[],double _r_e); 
//(Gabriel this is correct?)
/*********************************************************************************************************************/


/** \brief A struct to define elliptical sources
* 
*
*  x0 (\f$ x_0 \f$) -> position at x axis
*
*  y0 (\f$ y_0\f$) ->     position at y axis
*
*  R0 (\f$ R_0\f$) ->     characteristic size
*
*  eta0 (\f$ \eta_0\f$) ->    parameter relater to the ellipticity
*
*  theta0 (\f$ \theta_0\f$) ->  inclination to the main axis
*/
typedef struct {
  double x0 ;     //position at x axis
  double y0 ;     //position at y axis
  double R0 ;     //characteristic size
  double eta0 ;   //parameter relater to the ellipticity
  double theta0 ; //inclination to the main axis
} elliptical_source;


double f1_bar(f_type f1_in, elliptical_source source_in, double theta, double pert_params[], double _r_e);

double Df0Dtheta_bar(f_type Df0Dtheta_in, elliptical_source source_in, double theta, double pert_params[], double _r_e);

double arg_sqrt(f_type Df0Dtheta_in, elliptical_source source_in, double theta, double pert_params[], double _r_e);

double dr_plus(f_type f1_in, f_type Df0Dtheta_in, double kappa2, elliptical_source source_in, double theta, double pert_params[],double _r_e);

double dr_minus(f_type f1_in, f_type Df0Dtheta_in, double kappa2, elliptical_source source_in, double theta, double pert_params[],double _r_e);

double r_crit(f_type f1_in, f_type D2f0Dtheta2_in, double kappa2, double theta, double pert_params[],double _r_e);

double caustic_y1(f_type Df0Dtheta_in, f_type D2f0Dtheta2_in, double theta, double pert_params[],double _r_e);

double caustic_y2(f_type Df0Dtheta_in, f_type D2f0Dtheta2_in, double theta, double pert_params[],double _r_e);

double r_caust(f_type Df0Dtheta_in, f_type D2f0Dtheta2_in, double theta, double pert_params[],double _r_e);

double y1_leq(f_type f1_in, f_type Df0Dtheta_in, double x, double theta, double kappa2, double pert_params[],double _r_e);

double y2_leq(f_type f1_in, f_type Df0Dtheta_in, double x, double theta, double kappa2, double pert_params[],double _r_e);

double x_th(f_type f1_in, f_type D2f0Dtheta2_in, double kappa2, double theta, double pert_params[], double r_th, double _r_e);

double y1_th(f_type f1_in, f_type Df0Dtheta_in, f_type D2f0Dtheta2_in, double theta, double kappa2, double pert_params[],double _r_e, double r_th);

double y2_th(f_type f1_in, f_type Df0Dtheta_in, f_type D2f0Dtheta2_in, double theta, double kappa2, double pert_params[],double _r_e, double r_th);

double r_mean(f_type f1_in, f_type Df0Dtheta_in, double kappa2, elliptical_source source_in, double theta, double pert_params[],double _r_e);

double w_mean(f_type f1_in, f_type Df0Dtheta_in, double kappa2, elliptical_source source_in, double theta, double pert_params[],double _r_e);

double wr_theta(f_type f1_in, f_type Df0Dtheta_in, double kappa2, elliptical_source source_in, double theta, double pert_params[],double _r_e);

double y1_src(elliptical_source source_in, double theta);

double slc_y2_src(elliptical_source source_in, double theta);

#endif
