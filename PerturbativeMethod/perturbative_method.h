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



/**  \f$ \overline{f}_1(\theta) = f_1(\theta) + x0\cos(\theta) + y0\sin(\theta) \f$
*
*  \param f1_in function related to the perturbed potential 
*  \param source_in a struct that define elliptical sources
*  \param theta angular coordinate in the lens plane 
*  \param pert_params[] a vector that contains all the perturbation parameters
*  \return \f$ f_1(\theta) + x0\cos(\theta) + y0\sin(\theta) \f$
*
*  \sa f_type, elliptical_source
*/
double f1_bar(f_type f1_in, elliptical_source source_in, double theta, double pert_params[], double _r_e){
  return f1_in(theta, pert_params,_r_e) + source_in.x0*cos(theta) + source_in.y0*sin(theta);
}

/**  \f$ \frac{d \bar{f}_0(\theta)}{d \theta} \equiv \frac{d f_0(\theta)}{d \theta} - (x0\sin(\theta) - y0\cos(\theta))R_{_{\mathrm{E}}} \f$
*
*  \param Df0Dtheta_in function related to the perturbed potential
*  \param source_in a struct that define elliptical sources
*  \param theta angular coordinate in the lens plane 
*  \param pert_params[] a vector that contains all the perturbation parameters
*   \param _r_e is the einstein radius (If is not write, the code assumes \f$ R_{_{\mathrm{E}}}=1 \f$ )
*  \return \f$ \frac{d f_0(\theta)}{d \theta} - (x0\sin(\theta) - y0\cos(\theta))R_{_{\mathrm{E}}} \f$
*
*  \sa f_type, elliptical_source
*/
double Df0Dtheta_bar(f_type Df0Dtheta_in, elliptical_source source_in, double theta, double pert_params[], double _r_e){
  return Df0Dtheta_in(theta, pert_params,_r_e) + (source_in.y0*cos(theta) - source_in.x0*sin(theta))*_r_e ;
}

/**  Argument of square root: \f$ \Delta(\eta_0,\theta) = R_{0}^2(1-\eta_0\cos{2(\theta - \theta_0)})-(1-\eta_0^2)\left[\frac{1}{R_{_{\mathrm{E}}}}\dfrac{d\bar{f}_0}{d \theta} \right]^2 \f$ 
*
*  \param Df0Dtheta_in function related to the perturbed potential
*  \param source_in a struct that define elliptical sources
*  \param theta angular coordinate in the lens plane 
*  \param pert_params[] a vector that contains all the perturbation parameters
*  \param _r_e a Einstein Radius (If is not write, the code assumes \f$ R_{_{\mathrm{E}}}=1 \f$)
*  \return \f$  S\,R_{0}^2-(1-\eta_0^2)\left[\frac{1}{R_{_{\mathrm{E}}}}\dfrac{d\bar{f}_0}{d \theta} \right]^2,\qquad S=1-\eta_0\cos{2\tilde{\theta}},\quad \tilde{\theta}=\theta-\theta_0  \f$
*
*  \sa f_type, elliptical_source
*/
double arg_sqrt(f_type Df0Dtheta_in, elliptical_source source_in, double theta, double pert_params[], double _r_e){
  double Df0Dtheta_bar_tmp = Df0Dtheta_bar(Df0Dtheta_in, source_in, theta, pert_params,_r_e)/_r_e;
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
\return \f$ \frac{1}{\kappa_2}\left\{ \overline{f}_1(\theta) + \dfrac{\eta_0\sin{2\tilde{\theta}}}{S}\left(\dfrac{1}{R_{_\mathrm{E}}}\dfrac{\bar{f}_0(\theta)}{d\theta}\right)+\dfrac{1}{S}\sqrt{\Delta(\eta_0,\theta)}\right\}\f$*
*  \sa f_type, elliptical_source
*/
double dr_plus(f_type f1_in, f_type Df0Dtheta_in, double kappa2, elliptical_source source_in, double theta, double pert_params[],double _r_e){
  double Df0Dtheta_bar_tmp = Df0Dtheta_bar(Df0Dtheta_in, source_in, theta, pert_params,_r_e)/_r_e;
  double Df0Dtheta_bar_tmp2 = Df0Dtheta_bar_tmp*Df0Dtheta_bar_tmp;

  double f1_bar_tmp = f1_bar(f1_in, source_in, theta, pert_params,_r_e);

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
*
*  \return \f$ \frac{1}{\kappa_2}\left\{ \overline{f}_1(\theta) + \dfrac{\eta_0\sin{2\tilde{\theta}}}{S}\left(\dfrac{1}{R_{_\mathrm{E}}}\dfrac{\bar{f}_0(\theta)}{d\theta}\right)-\dfrac{1}{S}\sqrt{\Delta(\eta_0,\theta)}\right\}\f$
*
*  \sa f_type, elliptical_source
*/
double dr_minus(f_type f1_in, f_type Df0Dtheta_in, double kappa2, elliptical_source source_in, double theta, double pert_params[],double _r_e){
  double Df0Dtheta_bar_tmp = Df0Dtheta_bar(Df0Dtheta_in, source_in, theta, pert_params,_r_e)/_r_e;
  double Df0Dtheta_bar_tmp2 = Df0Dtheta_bar_tmp*Df0Dtheta_bar_tmp;

  double f1_bar_tmp = f1_bar(f1_in, source_in, theta, pert_params, _r_e);

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
*  \return \f$ R_{\mathrm{crit}}= R_{_\mathrm{E}} +\dfrac{1}{\kappa_2}\left[f_1+\dfrac{1}{R_{_\mathrm{E}}}\dfrac{d^2 f_0}{d \theta^2}\right]\f$
**/


double r_crit(f_type f1_in, f_type D2f0Dtheta2_in, double kappa2, double theta, double pert_params[],double _r_e){

  double f1_temp=f1_in(theta,pert_params,_r_e);
  double D2f0Dtheta2_temp= D2f0Dtheta2_in(theta,pert_params,_r_e)/_r_e;
  double r_t=_r_e+(1.0/kappa2)*(f1_temp+D2f0Dtheta2_temp);
  return r_t;
}


/** Tangential Caustic Line (First Component)
* 
*  \param D2f0Dtheta2_in function related to the perturbed potential (second derivative)
*  \param f1_in function related to the perturbed potential
*  \param theta angular coordinate in the lens plane 
*   This caustic line, in polar coordinates, is defined by :
*   \return \f$ y_1=\dfrac{1}{R_{_{\mathrm{E}}}}\frac{d^2f_0}{d \theta^2}\cos{\theta}+\dfrac{1}{R_{_{\mathrm{E}}}}\frac{df_0}{d \theta}\sin{\theta}\f$
**/
double caustic_y1(f_type Df0Dtheta_in, f_type D2f0Dtheta2_in, double theta, double pert_params[],double _r_e){
    double Df0Dtheta_tmp = Df0Dtheta_in(theta,pert_params,_r_e)/_r_e;
    double D2f0Dtheta2_temp= D2f0Dtheta2_in(theta,pert_params,_r_e)/_r_e;
    double y1_c=D2f0Dtheta2_temp*cos(theta)+Df0Dtheta_tmp*sin(theta);
    return y1_c;
}


/** Tangential Caustic Line (Second Component)
*  \param D2f0Dtheta2_in function related to the perturbed potential (second derivative)
*  \param f1_in function related to the perturbed potential
*  \param theta angular coordinate in the lens plane 
*   This caustic line, in polar coordinates, is defined by :
*   \return \f$ y_2=\dfrac{1}{R_{_\mathrm{E}}}\frac{d^2f_0}{d \theta^2}\sin{\theta}-\frac{1}{R_{_\mathrm{E}}}\frac{df_0}{d \theta}\cos{\theta}\f$
**/
double caustic_y2(f_type Df0Dtheta_in, f_type D2f0Dtheta2_in, double theta, double pert_params[],double _r_e){
    double Df0Dtheta_tmp = Df0Dtheta_in(theta,pert_params,_r_e)/_r_e;
    double D2f0Dtheta2_temp= D2f0Dtheta2_in(theta,pert_params,_r_e)/_r_e;
    double y2_c=D2f0Dtheta2_temp*sin(theta)-Df0Dtheta_tmp*cos(theta);
    return y2_c;
}


/** First component of the parametric equation for the elliptical souce not allignet to main axis
*   This caustic line, in polar coordinates, is defined by :
*
*  \f$ y_{1s}= \dfrac{R_s}{\sqrt{1-\eta_s}}\cos{\theta},\qquad y_{2s}= \dfrac{R_s}{\sqrt{1+\eta_s}}\sin{\theta} \f$
*
*  \param elliptical_source source parameters 
*  \param theta angular coordinate in the source plane 
*
*   \return \f$ y^{\prime}_{1,src}= x_0+y_{1s}\cos{\theta_s}-y_{2s}\sin{\theta_s} \f$ 
**/


double y1_src(elliptical_source source_in, double theta){
      double theta_0 = source_in.theta0;  
      double y1s=source_in.R0*cos(theta)/sqrt(1.-source_in.eta0);
      double y2s=source_in.R0*sin(theta)/sqrt(1.+source_in.eta0);
      double y1_linha_src=y1s*cos(theta_0)-y2s*sin(theta_0)+source_in.x0;
      return y1_linha_src;
}


/** Second component of the parametric equation for the elliptical souce not allignet to main axis
*   This caustic line, in polar coordinates, is defined by :
*
*  \f$ y_{1s}= \dfrac{R_s}{\sqrt{1-\eta_s}}\cos{\theta},\qquad y_{2s}= \dfrac{R_s}{\sqrt{1+\eta_s}}\sin{\theta} \f$
*  \param elliptical_source source parameters 
*  \param theta  angular coordinate in the source plane 
*   \return \f$ y^{\prime}_{2,src}= y_0+y_{1s}\sin{\theta_s}+y_{2s}\cos{\theta_s} \f$ 
**/

double y2_src(elliptical_source source_in, double theta){
      double theta_0 = source_in.theta0;  
      double y1s=source_in.R0*cos(theta)/sqrt(1.-source_in.eta0);
      double y2s=source_in.R0*sin(theta)/sqrt(1.+source_in.eta0);
      double y2_linha_src=y1s*sin(theta_0)+y2s*cos(theta_0)+source_in.y0 ;
      return y2_linha_src;
}

/** Mean value of the radial distance of the arcs (Eq. 2.68 of the Report) 
*   This function must be documented well
*  \param elliptical_source (source parameters)
*  \param theta  angular coordinate in the source plane 
*   \return \f$ \bar{r}(\theta) \f$ 
**/

double bar_r(f_type f1_in, f_type Df0Dtheta_in, double kappa2, elliptical_source source_in, double theta, double pert_params[],double _r_e){
  return _r_e+(dr_plus(f1_in,Df0Dtheta_in, kappa2, source_in,theta, pert_params,_r_e)+dr_minus(f1_in, Df0Dtheta_in, kappa2, source_in,theta, pert_params, _r_e))/2.0;
}

/** Width of the image (Eq. 2.69 of the Repor)
*
*  \param elliptical_source (source parameters)
*  \param theta  angular coordinate in the source plane 
*   \return \f$ W(\theta) \f$ 
**/

double w_mean(f_type f1_in, f_type Df0Dtheta_in, double kappa2, elliptical_source source_in, double theta, double pert_params[],double _r_e){
  return (dr_plus(f1_in,Df0Dtheta_in, kappa2, source_in, theta, pert_params, _r_e)-dr_minus(f1_in, Df0Dtheta_in, kappa2, source_in,theta, pert_params, _r_e));
}

/** Unitary area of the image, i.e, \f$ W(\theta)\bar{r}(\theta)\f$ (see Eq. 2.87 of the Report)
*
*  \param elliptical_source (source parameters)
*  \param theta  angular coordinate in the source plane 
*   \return \f$ W(\theta) \f$ 
**/
double wr_theta(f_type f1_in, f_type Df0Dtheta_in, double kappa2, elliptical_source source_in, double theta, double pert_params[],double _r_e){
  double wr_temp=bar_r(f1_in,Df0Dtheta_in, kappa2, source_in,theta, pert_params,_r_e)*w_mean(f1_in, Df0Dtheta_in, kappa2, source_in, theta, pert_params, _r_e);
  return wr_temp;  
}
#endif
