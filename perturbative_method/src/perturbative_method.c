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

#include "perturbative_method.h"
#include "config.h"

/**  \f$ \overline{f}_1(\theta) = f_1(\theta) + x0\cos(\theta) + y0\sin(\theta) \f$
*
*  \param f1_in function related to the perturbed potential 
*  \param source_in a struct that define elliptical sources
*  \param theta angular coordinate in the lens plane 
*  \param pert_params[] a vector containing all the perturbation parameters
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
*  \param pert_params[] a vector containing all the perturbation parameters
*   \param _r_e is the einstein radius (\f$ R_{_{\mathrm{E}}}\f$ )
*  \return \f$ \frac{d f_0(\theta)}{d \theta} - (x0\sin(\theta) - y0\cos(\theta))R_{_{\mathrm{E}}} \f$
*
*  \sa f_type, elliptical_source
*/
double Df0Dtheta_bar(f_type Df0Dtheta_in, elliptical_source source_in, double theta, double pert_params[], double _r_e){
  return Df0Dtheta_in(theta, pert_params,_r_e) + (source_in.y0*cos(theta) - source_in.x0*sin(theta))*_r_e ;
}

/**  Argument of square root: \f$ \Delta(\eta_0,\theta) = R_{0}^2(1-\eta_0\cos{2(\theta - \theta_0)})-(1-\eta_0^2)\left[\frac{1}{R_{_{\mathrm{E}}}}\frac{d\bar{f}_0}{d \theta} \right]^2 \f$ 
*
*  \param Df0Dtheta_in function related to the perturbed potential
*  \param source_in a struct that define elliptical sources
*  \param theta angular coordinate in the lens plane 
*  \param pert_params[] a vector containing all the perturbation parameters
*   \param _r_e is the einstein radius (\f$ R_{_{\mathrm{E}}}\f$ )
*  \return \f$  S\,R_{0}^2-(1-\eta_0^2)\left[\frac{1}{R_{_{\mathrm{E}}}}\frac{d\bar{f}_0}{d \theta} \right]^2,\qquad S=1-\eta_0\cos{2\tilde{\theta}},\quad \tilde{\theta}=\theta-\theta_0  \f$
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
*  \param pert_params[] a vector containing all the perturbation parameters
*   \param _r_e is the einstein radius (\f$ R_{_{\mathrm{E}}}\f$ )
\return \f$ \frac{1}{\kappa_2}\left\{ \overline{f}_1(\theta) + \frac{\eta_0\sin{2\tilde{\theta}}}{S}\left(\frac{1}{R_{_\mathrm{E}}}\frac{\bar{f}_0(\theta)}{d\theta}\right)+\frac{1}{S}\sqrt{\Delta(\eta_0,\theta)}\right\}\f$*
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
*  \param pert_params[] a vector containing all the perturbation parameters
*   \param _r_e is the einstein radius (\f$ R_{_{\mathrm{E}}}\f$ )
*
*  \return \f$ \frac{1}{\kappa_2}\left\{ \overline{f}_1(\theta) + \frac{\eta_0\sin{2\tilde{\theta}}}{S}\left(\frac{1}{R_{_\mathrm{E}}}\frac{\bar{f}_0(\theta)}{d\theta}\right)-\frac{1}{S}\sqrt{\Delta(\eta_0,\theta)}\right\}\f$
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
*  \param pert_params[] a vector containing all the perturbation parameters
*   \param _r_e is the einstein radius (\f$ R_{_{\mathrm{E}}}\f$ ) of the circular model
*  \return \f$ R_{\mathrm{crit}}= R_{_\mathrm{E}} +\frac{1}{\kappa_2}\left[f_1+\frac{1}{R_{_\mathrm{E}}}\frac{d^2 f_0}{d \theta^2}\right]\f$
* \sa f_type
*/
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
*  \param pert_params[] a vector containing all the perturbation parameters
*   \param _r_e is the einstein radius (\f$ R_{_{\mathrm{E}}}\f$ ) of the circular model
*   This caustic line, in polar coordinates, is defined by :
*   \return \f$ y_1=\frac{1}{R_{_{\mathrm{E}}}}\frac{d^2f_0}{d \theta^2}\cos{\theta}+\frac{1}{R_{_{\mathrm{E}}}}\frac{df_0}{d \theta}\sin{\theta}\f$
* \sa f_type
*/
double caustic_y1(f_type Df0Dtheta_in, f_type D2f0Dtheta2_in, double theta, double pert_params[],double _r_e){
    double Df0Dtheta_tmp = Df0Dtheta_in(theta,pert_params,_r_e)/_r_e;
    double D2f0Dtheta2_temp= D2f0Dtheta2_in(theta,pert_params,_r_e)/_r_e;
    double y1_c=D2f0Dtheta2_temp*cos(theta)+Df0Dtheta_tmp*sin(theta);
    return y1_c;
}


/** Tangential Caustic Line (Second Component)
*  \param D2f0Dtheta2_in function related to the perturbed potential (second derivative in relation to \f$ theta \f$)
*  \param f1_in function related to the perturbed potential
*  \param theta angular coordinate in the lens plane 
*  \param pert_params[] a vector containing all the perturbation parameters
*   \param _r_e is the einstein radius (\f$ R_{_{\mathrm{E}}}\f$ ) of the circular model
*   This caustic line, in polar coordinates, is defined by :
*   \return \f$ y_2=\frac{1}{R_{_\mathrm{E}}}\frac{d^2f_0}{d \theta^2}\sin{\theta}-\frac{1}{R_{_\mathrm{E}}}\frac{df_0}{d \theta}\cos{\theta}\f$
**/
double caustic_y2(f_type Df0Dtheta_in, f_type D2f0Dtheta2_in, double theta, double pert_params[],double _r_e){
    double Df0Dtheta_tmp = Df0Dtheta_in(theta,pert_params,_r_e)/_r_e;
    double D2f0Dtheta2_temp= D2f0Dtheta2_in(theta,pert_params,_r_e)/_r_e;
    double y2_c=D2f0Dtheta2_temp*sin(theta)-Df0Dtheta_tmp*cos(theta);
    return y2_c;
}

double r_caust(f_type Df0Dtheta_in, f_type D2f0Dtheta2_in, double theta, double pert_params[],double _r_e){
    double y1_temp=caustic_y1( Df0Dtheta_in, D2f0Dtheta2_in, theta, pert_params, _r_e);  
    double y2_temp=caustic_y2( Df0Dtheta_in, D2f0Dtheta2_in, theta, pert_params, _r_e);  
    double y2_c=sqrt(pow(y1_temp,2)+pow(y2_temp,2));
    return y2_c;
}


/** First components of the lens equation for the perturbative method
*
*  \f$ y_{1_s}= (\kappa_2 x-f_1)\cos{\theta}+\frac{1}{R_{\rm E}}\frac{\partial f_0}{\partial \theta}\sin{\theta}\f$
*
*  \param x : radial coordinate in the lens plane 
*  \param theta : angular coordinate in the lens plane 
*  \param Df0Dtheta_in :function related to the perturbed potential
*  \param elliptical_source :source parameters 
*  \param kappa_2 :parameter related with the convergence of the circular model
*  \param _r_e : Einstein radius of the lens model
*  \param sflag : Integer to decide wheter the source is needed (sflag=1) or not (sflag=0)	
*   \return \f$ y^{1_s}(\theta) \f$ 
*  \sa  f_type, elliptical_source
*/

double y1_leq(f_type f1_in, f_type Df0Dtheta_in, double x, double theta, double kappa2, double pert_params[],double _r_e){
  double f1_temp, df0dtheta_temp;
 f1_temp= f1_in(theta,pert_params,_r_e);
  df0dtheta_temp= Df0Dtheta_in(theta,pert_params,_r_e);
   return (kappa2*x-f1_temp)*cos(theta)+(1./_r_e)*df0dtheta_temp*sin(theta);
  
}


/** Second components of the lens equation for the perturbative method
*
*  \f$ y_{2_s}= (\kappa_2 x-f_1)\sin{\theta}-\frac{1}{R_{\rm E}}\frac{\partial f_0}{\partial \theta}\cos{\theta}\f$
*
*  \param x : radial coordinate in the lens plane 
*  \param theta : angular coordinate in the lens plane 
*  \param Df0Dtheta_in :function related to the perturbed potential
*  \param elliptical_source :source parameters 
*  \param kappa_2 :parameter related with the convergence of the circular model
*  \param _r_e : Einstein radius of the lens model
*  \param sflag : Integer to decide wheter the source is needed (sflag=1) or not (sflag=0)	
*   \return \f$ y_{2_s}(\theta) \f$ 
*  \sa  f_type, elliptical_source
*/

double y2_leq(f_type f1_in, f_type Df0Dtheta_in, double x, double theta, double kappa2, double pert_params[],double _r_e){
  double f1_temp, df0dtheta_temp;
  f1_temp=f1_in(theta,pert_params,_r_e);
  df0dtheta_temp=Df0Dtheta_in(theta,pert_params,_r_e);
  
  return (kappa2*x-f1_temp)*sin(theta)-(1./_r_e)*df0dtheta_temp*cos(theta);
  
}



/**  Radial compontent (\f$ x_{\rm th} \f$ ) in the lens plane of the Constant Distortion Curves (take care: \f$ r_{\rm th}=x_{\rm th}+R_{\rm E}\f$ ) 
*
*  To plot the tangential critical curve, use  polar coordinates i.e.\f$ x_{c_1}=r_{\mathrm{crit}}*\cos{\theta},x_{c_2}=r_{\mathrm{crit}}*\sin{\theta}\f$  
*  \param D2f0Dtheta2_in function related to the perturbed potential (second derivative)
*  \param f1_in function related to the perturbed potential
*  \param kappa2 \f$  \kappa_2 = 1 - \left[\frac{d^2 \phi_0(r)}{dr^2}\right]_{r=r_e} \f$
*  \param theta angular coordinate in the lens plane.
*  \param pert_params[] a vector containing all the perturbation parameters
*  \param r_th : length-to-width threshold (could be positive or negative)
*  \param _r_e is the einstein radius (\f$ R_{_{\mathrm{E}}}\f$ ) of the circular model
*  \return \f$ R_{\mathrm{crit}}= R_{_\mathrm{E}} +\frac{1}{\kappa_2}\left[f_1+\frac{1}{R_{_\mathrm{E}}}\frac{d^2 f_0}{d \theta^2}\right]\f$
* \sa f_type
*/
double x_th(f_type f1_in, f_type D2f0Dtheta2_in, double kappa2, double theta, double pert_params[], double r_th, double _r_e){
  double q=1.-1/r_th;
  double x_cdc=r_crit(f1_in,D2f0Dtheta2_in,kappa2,theta, pert_params, _r_e)/q-_r_e;
  return x_cdc;
}

/** First components of the constant distortion curves in the source 
*
*  \f$ y_{1{\rm th}}= (\kappa_2 x-f_1)\cos{\theta}+\frac{1}{R_{\rm E}}\frac{\partial f_0}{\partial \theta}\sin{\theta}\f$
*
*  \param theta : angular coordinate in the lens plane 
*  \param f1_in : function related to the perturbed potential
*  \param Df0Dtheta_in :function related to the perturbed potential
*  \param D2f0Dtheta2_in :function related to the perturbed potential
*  \param elliptical_source :source parameters 
*  \param kappa_2 :parameter related with the convergence of the circular model
*  \param _r_e : Einstein radius of the lens model
*  \param sflag : Integer to decide wheter the source is needed (sflag=1) or not (sflag=0)	
*   \return \f$ y^{1_s}(\theta) \f$ 
*  \sa  f_type, elliptical_source
*/

double y1_th(f_type f1_in, f_type Df0Dtheta_in, f_type D2f0Dtheta2_in, double theta, double kappa2, double pert_params[],double _r_e, double r_th){
  double x_cor=x_th(f1_in,D2f0Dtheta2_in, kappa2, theta, pert_params, r_th, _r_e);
  return y1_leq(f1_in,Df0Dtheta_in, x_cor, theta, kappa2, pert_params,_r_e);
  
}

/** Second components of the constant distortion curves in the source 
*
*  \f$ y_{2 {\rm th}}= (\kappa_2 x_{\rm th}-f_1)\sin{\theta}-\frac{1}{R_{\rm E}}\frac{\partial f_0}{\partial \theta}\cos{\theta}\f$
*
*  \param theta : angular coordinate in the lens plane 
*  \param Df0Dtheta_in :function related to the perturbed potential
*  \param D2f0Dtheta2_in :function related to the perturbed potential  
*  \param kappa_2 :parameter related with the convergence of the circular model
*  \param _r_e : Einstein radius of the lens model
*   \return \f$ y^{1_s}(\theta) \f$ 
*  \sa  f_type, x_th
*/

double y2_th(f_type f1_in, f_type Df0Dtheta_in, f_type D2f0Dtheta2_in, double theta, double kappa2, double pert_params[],double _r_e, double r_th){
  double x_cor=x_th(f1_in,D2f0Dtheta2_in, kappa2, theta, pert_params, r_th, _r_e);
  return y2_leq(f1_in,Df0Dtheta_in, x_cor, theta, kappa2, pert_params,_r_e);
  
}


/** Mean value of the radial distance of the arcs (Eq. 2.68 of the Report) 
*
*  \param f_1in: function related to the first derivative of the perturbed potential
*  \param Df0Dtheta_in: function related to the perturbative potential (first derivative in relation to \f$ theta \f$)
*  \param elliptical_source (source parameters)
*  \param theta  angular coordinate in the source plane 
*  \param pert_params[] a vector containing all the perturbation parameters
*   \param _r_e is the einstein radius (\f$ R_{_{\mathrm{E}}}\f$ ) of the circular model
*   \return \f$ \bar{r}(\theta) = R_{_\mathrm{E}}+\left(\frac{x_{+}(\theta)+x_{-}(\theta)}{2}\right)\f$ 
*   \sa f_type, elliptical_source,  
*/

double r_mean(f_type f1_in, f_type Df0Dtheta_in, double kappa2, elliptical_source source_in, double theta, double pert_params[],double _r_e){
  return _r_e+(dr_plus(f1_in,Df0Dtheta_in, kappa2, source_in,theta, pert_params,_r_e)+dr_minus(f1_in, Df0Dtheta_in, kappa2, source_in,theta, pert_params, _r_e))/2.0;
}

/** Mean width of the image at some angular position  (Eq. 2.69 of the Repor)
*
*  \param f_1in: function related to the first derivative of the perturbed potential
*  \param Df0Dtheta_in: function related to the perturbative potential (first derivative in relation to \f$ theta \f$)
*   \param kappa2: corresponds to \f$ \kappa_2 \f$
*  \param elliptical_source (source parameters)
*  \param theta  angular coordinate in the source plane 
*  \param pert_params[] a vector containing all the perturbation parameters
*   \param _r_e is the einstein radius (\f$ R_{\rm E}\f$ ) of the circular model
*   \return \f$ W(\theta)=x_{+}(\theta)-x_{-}(\theta) \f$ 
*   \sa f_type, elliptical_source, r_mean, w_mean 
*
*/

double w_mean(f_type f1_in, f_type Df0Dtheta_in, double kappa2, elliptical_source source_in, double theta, double pert_params[],double _r_e){
  return (dr_plus(f1_in,Df0Dtheta_in, kappa2, source_in, theta, pert_params, _r_e)-dr_minus(f1_in, Df0Dtheta_in, kappa2, source_in,theta, pert_params, _r_e));
}

/** Unitary area of the image, i.e, is the argument of the integral \f$ W(\theta)\bar{r}(\theta)\f$ (see Eq. 2.87 of the Report)
*  \param f_1in: function related to the first derivative of the perturbed potential
*  \param Df0Dtheta_in: function related to the perturbative potential (first derivative in relation to \f$ theta \f$)
*   \param kappa2: corresponds to \f$ \kappa_2 \f$
*  \param elliptical_source (source parameters)
*  \param theta  angular coordinate in the source plane 
*  \param pert_params[] a vector containing all the perturbation parameters
*   \param _r_e is the einstein radius (\f$ R_{\rm E}\f$ ) of the circular model
*   \return \f$ W(\theta)\bar{r}(\theta) \f$ 
*   \sa f_type, elliptical_source, r_mean, w_mean
**/
double wr_theta(f_type f1_in, f_type Df0Dtheta_in, double kappa2, elliptical_source source_in, double theta, double pert_params[],double _r_e){
  double wr_temp=r_mean(f1_in,Df0Dtheta_in, kappa2, source_in,theta, pert_params,_r_e)*w_mean(f1_in, Df0Dtheta_in, kappa2, source_in, theta, pert_params, _r_e);
  return wr_temp;  
}

/** First component of the parametric equation for the elliptical souce not allignet to main axis
*
*  \f$ y_{1s}= \frac{R_s}{\sqrt{1-\eta_s}}\cos{\theta},\qquad y_{2s}= \frac{R_s}{\sqrt{1+\eta_s}}\sin{\theta} \f$
*
*  \param f1_in function related to the perturbed potential
*  \param elliptical_source source parameters 
*  \param theta angular coordinate in the source plane 
*   \return \f$ y^{\prime}_{1,src}= x_0+y_{1s}\cos{\theta_s}-y_{2s}\sin{\theta_s} \f$ 
*  \sa elliptical_source
*/

double y1_src(elliptical_source source_in, double theta){
      double theta_0 = source_in.theta0;  
      double y1s=source_in.R0*cos(theta)/sqrt(1.-source_in.eta0);
      double y2s=source_in.R0*sin(theta)/sqrt(1.+source_in.eta0);
      double y1_linha_src=y1s*cos(theta_0)-y2s*sin(theta_0)+source_in.x0;
      return y1_linha_src;
}




/**
 * pm_y2_src:
 * @source_in: source_in.
 * @k: theta theta. 
 * 
 *  
 * Returns: a double.
 */
double y2_src(elliptical_source source_in, double theta){
      double theta_0 = source_in.theta0;  
      double y1s=source_in.R0*cos(theta)/sqrt(1.-source_in.eta0);
      double y2s=source_in.R0*sin(theta)/sqrt(1.+source_in.eta0);
      double y2_linha_src=y1s*sin(theta_0)+y2s*cos(theta_0)+source_in.y0 ;
      return y2_linha_src;
}

double slc_test_func(){
    return 0.0;
}

