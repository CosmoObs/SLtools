/** @file
* Example of doxygen documentation for C functions FIXME. 
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

//double _r_e = 1.0;

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
  return f1_in(theta, pert_params) + source_in.x0*cos(theta) + source_in.y0*sin(theta);
}


// Note for Gabriel: In the original version you wrote \overline{\frac{d f_1(\theta)}{d \theta}}, but the correct is d f_0 (\theta)
/**  \f$ \overline{\frac{d f_0(\theta)}{d \theta}} \equiv \frac{d f_0(\theta)}{d \theta} - x0\sin(\theta) + y0\cos(\theta) \f$
*
*  \param Df0Dtheta_in function related to the perturbed potential
*  \param source_in a struct that define elliptical sources
*  \param theta angular coordinate in the lens plane 
*  \param pert_params[] a vector that contains all the perturbation parameters
*  \return \f$ \frac{d f_0(\theta)}{d \theta} - x0\sin(\theta) + y0\cos(\theta) \f$
*
*  \sa f_type, elliptical_source
*/
double Df0Dtheta_bar(f_type Df0Dtheta_in, elliptical_source source_in, double theta, double pert_params[]){
  return Df0Dtheta_in(theta, pert_params) + (source_in.y0*cos(theta) - source_in.x0*sin(theta))*_r_e ;
}


/**  Argument of square root: \f$ R_{0}^2(1-\eta_0\cos(\theta - \theta_0))-(1-\eta_0^2) \overline{\frac{d f_1(\theta)}{d \theta}}^2 \f$ 
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
  double Df0Dtheta_bar_tmp = Df0Dtheta_bar(Df0Dtheta_in, source_in, theta, pot_params)/_r_e;
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
  double Df0Dtheta_bar_tmp = Df0Dtheta_bar(Df0Dtheta_in, source_in, theta, pot_params)/_r_e;
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
  double Df0Dtheta_bar_tmp = Df0Dtheta_bar(Df0Dtheta_in, source_in, theta, pot_params)/_r_e;
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


/**  Radius of the Tangential Critical Curve (r_crit)
*
*  To plot the tangential critical curve, use  polar coordinates i.e.\f$ x_{c_1}=r_{\mathrm{crit}}*\cos{\theta},x_{c_2}=r_{\mathrm{crit}}*\sin{\theta}\f$  
*  \param D2f0Dtheta2_in function related to the perturbed potential (second derivative)
*  \param f1_in function related to the perturbed potential
*  \param kappa2 \f$  \kappa_2 = 1 - \left[\frac{d^2 \phi_0(r)}{dr^2}\right]_{r=r_e} \f$
*  \param theta angular coordinate in the lens plane.
*  \return \f$ r_{\mathrm{crit}}= r_{\mathrm{E}} +\frac{1}{\kappa_2}\left[f_1+\frac{1}{r_{\mathrm{E}}}\frac{d^2 f_0}{d \theta^2}\right]\f$
**/


double r_crit(f_type f1_in, f_type D2f0Dtheta2_in, double kappa2, double theta, double pot_params[]){

  double f1_temp=f1_in(theta,pot_params);
  double D2f0Dtheta2_temp= D2f0Dtheta2_in(theta,pot_params)/_r_e;
  double r_t=_r_e+(1.0/kappa2)*(f1_temp+D2f0Dtheta2_temp);
  return r_t;
}


/** Tangential Caustic Line (First Component)
* 
*  \param D2f0Dtheta2_in function related to the perturbed potential (second derivative)
*  \param f1_in function related to the perturbed potential
*  \param theta angular coordinate in the lens plane 
*   This caustic line, in polar coordinates, is defined by :
*   \return \f$ y_1=\frac{1}{r_{\mathrm{E}}}\frac{d^2f_0}{d \theta^2}\cos{\theta}+\frac{1}{r_{\mathrm{E}}}\sin{\theta}\f$
**/
double caustic_y1(f_type Df0Dtheta_in, f_type D2f0Dtheta2_in, double theta, double pot_params[]){
    double Df0Dtheta_tmp = Df0Dtheta_in(theta,pot_params)/_r_e;
    double D2f0Dtheta2_temp= D2f0Dtheta2_in(theta,pot_params)/_r_e;
    double y1_c=(1.0/_r_e)*D2f0Dtheta2_temp*cos(theta)+(1.0/_r_e)*Df0Dtheta_tmp*sin(theta);
    return y1_c;
}


/** Tangential Caustic Line (Second Component)
*  \param D2f0Dtheta2_in function related to the perturbed potential (second derivative)
*  \param f1_in function related to the perturbed potential
*  \param theta angular coordinate in the lens plane 
*   This caustic line, in polar coordinates, is defined by :
*   \return \f$ y_2=\frac{1}{r_{\mathrm{E}}}\frac{d^2f_0}{d \theta^2}\sin{\theta}+\frac{1}{r_{\mathrm{E}}}\cos{\theta}\f$
**/
double caustic_y2(f_type Df0Dtheta_in, f_type D2f0Dtheta2_in, double theta, double pot_params[]){
    double Df0Dtheta_tmp = Df0Dtheta_in(theta,pot_params)/_r_e;
    double D2f0Dtheta2_temp= D2f0Dtheta2_in(theta,pot_params)/_r_e;
    double y2_c=(1.0/_r_e)*D2f0Dtheta2_temp*sin(theta)-(1.0/_r_e)*Df0Dtheta_tmp*cos(theta);
    return y2_c;
}
// Plotting a circular source

double y1_src(elliptical_source source_in, double theta){
      double theta_bar = theta - source_in.theta0;  
      return source_in.R0*cos(theta_bar)/sqrt(1.-source_in.eta0) + source_in.x0;
}


double y2_src(elliptical_source source_in, double theta){
      double theta_bar = theta - source_in.theta0;  
      return source_in.R0*sin(theta_bar)/sqrt(1.+source_in.eta0) + source_in.y0;
}



/** A functions to print the arcs
*
*  \param source_in a struct that define elliptical sources
*  \param Df0Dtheta_in function related to the perturbed potential
*  \param f1_in function related to the perturbed potential
*  \param Df0Dtheta_in function related to the perturbed potential
*  \param pert_params[] a vector that contains all the perturbation parameters
*  \param kappa2 \f$  \kappa_2 = 1 - \left[\frac{d^2 \phi_0(r)}{dr^2}\right]_{r=r_e} \f$
*  \param npts number of theta divisions 
*  \param file_in a file to write the arc's contours
*  \return nothing
*
*  \sa f_type, dr_plus, dr_minus
*/
void plot_arcs(elliptical_source source_in, f_type f1_in, f_type Df0Dtheta_in, double pert_params[], double kappa2, int npts=500, FILE *file_in=NULL){
  double theta = 0.0;
  double r_p = 0.0;
  double r_m = 0.0;

  for(int i=0;i<=npts;i++){

    if(arg_sqrt(Df0Dtheta_in, source_in, theta, pert_params)>0.0){
      r_p = _r_e+dr_plus(f1_in, Df0Dtheta_in, kappa2, source_in, theta, pert_params);
      r_m = _r_e+dr_minus(f1_in, Df0Dtheta_in, kappa2, source_in, theta, pert_params);
      if(file_in==NULL){
        printf("%E %E\n",r_p*cos(theta),r_p*sin(theta));
        printf("%E %E\n",r_m*cos(theta),r_m*sin(theta));
      } else {
        fprintf(file_in,"%E %E\n",r_p*cos(theta),r_p*sin(theta));
        fprintf(file_in,"%E %E\n",r_m*cos(theta),r_m*sin(theta));
      }
    }
    theta+= 6.283185308/npts;
  }

}






#endif
