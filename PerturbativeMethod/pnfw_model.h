/** @file
* Example of doxygen documentation for C functions FIXME. 
*/

/** @package nfw_model
*  Package to compute quantities related to NFW models.
*
*  Detailed descrition FIXME
*
*/
#ifndef NFW_MODEL_H
#define NFW_MODEL_H

//to compile: g++ -Wall nfw_pert.cpp
/*******************************************************************************
 Lensing Functions related to the NFW model
*******************************************************************************/
/***************************************************************************************************************/
//! Function useful to define the convergence and angle deflection of the NFW model.
//!
//! \f$ \mathcal{F}(x)=\left\{\begin{array}{lc} \dfrac{1}{\sqrt{1-x^2}}\tanh^{-1}{(\sqrt{1-x^2})} & (x<1)\\  \\ 1 & (x=1)\\ \\ \dfrac{1}{\sqrt{x^2-1}}\tan^{-1}{(\sqrt{x^2-1})} & (x>1) \end{array}\right.  \f$, 
/*!   \param x: dimensionless coordinate   \return \f$ \mathcal{F}(x)\f$ */

double F(double X)
{
    double f;
    double x2;
    double arg;
    
    x2 = X*X;
    
    if (X < 1.0 - 1e-8)
    {
        arg = sqrt(1.0-x2);
        
        f = 1.0/arg * atanh(arg);
    }
    else if (X > 1.0 + 1e-8)
    {
        arg = sqrt(x2 - 1.0);
        
        f = 1.0/arg * atan(arg);
    }
    else
    {
        f = 1.0;
    }
    
    return f;
}
/***************************************************************************************************************/
//! Convergence of the NFW model.
//!
//! \f$ \kappa(r)=2\,\kappa_s\dfrac{(1-\mathcal{F}(x))}{x^2-1}\f$, 
//!
//! where \f$ \kappa_s \f$ is the characteristic convergence and \f$ \mathcal{F}(x)\f$ is defined in F(double x)
/*!   \param r : radial distance,   \param pot_params[] : NFW lens parameters,   \return \f$ \kappa(r)\f$ */
double conv_nfw(double r, double pot_params[])
{
    double K, ks=pot_params[1];
    double rs=pot_params[0],X=r/rs;
    
    if(X == 1.0)
    {
        return 2.0/3.0*ks;
    }
    else
    {
        K = 2.0*ks*(1.0 - F(X))/(X*X - 1.0);
        return K;
    }    
}
/***************************************************************************************************************/
//! Angle deflection of the NFW model.
//!
//! \f$ \alpha(r)=   2\,\kappa_s r_s \dfrac{\log{x/2}+\mathcal{F}(x)}{x}\f$, 
//!  where \f$ \kappa_s \f$ is the characteristic convergence, \f$ r_s \f$ is the scale radius, \f$ x=r/r_s \f$  and \f$ \mathcal{F}(x)\f$ is defined in F(double x)
/*!
  \param r : radial distance,\param pot_params[] : NFW lens parameters 
  \return \f$ \alpha(r)\f$
*/
double alpha_nfw(double r, double pot_params[])
{
    double alpha,rs=pot_params[0],ks=pot_params[1];
    double X=r/rs;
      alpha = 4.0*ks*rs*(log(X/2.0)+F(X))/X;
      return alpha;
}
/***************************************************************************************************************/
//! Shear of the NFW model.
//!
//! \f$ \gamma(r)= \dfrac{\alpha(r)}{r}-\kappa(r) \f$, 
//!
//!  where \f$ \alpha(r) \f$ is angle deflection and \f$ \kappa(r) \f$ is the convergence.
//!
/*!   \param r : radial distance,   \param pot_params[] : NFW lens parameters, \return \f$ \gamma(r)\f$ */
double shear_nfw(double r, double pot_params[]){
      double gamma,rs=pot_params[0];
      double X=r/rs;
      gamma=alpha_nfw(r, pot_params)/r-conv_nfw(r, pot_params);
      return gamma;
}
/********************************************************************************************************************/
//! Solution for the Circular NFW Lens 
//!
//! \f$ \kappa_2=2-2\kappa(R_{\mathrm{E}}) \f$,
//!
//! where \f$ \kappa(R_{\mathrm{E}})\f$ is the convergence for evaluated at the Einstein Radius
//!
/*!   \param r_e : Einstein Radius, \param pot_params[] : NFW lens parameters,\return \f$ \kappa_2 \f$ */

double kappa2_nfw(double pot_params[]){
  double k2,r=_r_e;

  k2=2-2*conv_nfw(r,pot_params);

  return k2;
}

/*********************************************************************************************************************/
//!  First perturbative field for the PNFW model as Perturbation.
//!
//!  \f$ \psi_E(r)=\phi_E(\xi)-\phi_0(r) \f$ 
//!
//!  \f$ \xi=r\sqrt{a_{1\eta}\cos^2{\theta}+a_{2\eta}\sin^2{\theta}}\f$ 
//!        (\f$ \xi_E\f$ is the elliptical variable \f$ \xi(r,\theta)\f$ for \f$ r=R_{E}\f$) 
//!
//!  \f$ f_{1}(\theta) \equiv \left[\dfrac{\partial\psi_E(r,\theta)}{\partial r}\right]_{r=R_{E}} = \dfrac{\xi_{\mathrm{E}}}{R_{\mathrm{E}}}\alpha(\xi_{\mathrm{E}})-\alpha(R_{\mathrm{E}})\f$
/*!
  \param theta : angular coordinate in the lens plane,
  \param pot_params[] : PNFW lens parameters
  \return \f$ f_1(\theta)\f$
*/
double f1_nfw(double theta, double pot_params[]){
    double eta=pot_params[2];
    double a_eta1=1.0-eta,a_eta2=1.0/(1.0-eta); // If you want work with the standar parameterization
//     double a_eta1=1.0-eta,a_eta2=1.0+eta;   // If you want work with the Angle Deflection Method
    double xie=_r_e*sqrt(a_eta1*pow(cos(theta),2)+a_eta2*pow(sin(theta),2));
    double f1=(xie/_r_e)*alpha_nfw(xie,pot_params)-alpha_nfw(_r_e,pot_params);
    return f1;
}
/*********************************************************************************************************************/
//!  Second perturbative field for the PNFW model as Perturbation (see PerturbedPNFW.pdf) 
//!
//!  \f$ \psi_E(r)=\phi_E(\xi)-\phi_0(r) \f$
//!
//!  \f$ \xi=r\sqrt{a_{1\eta}\cos^2{\theta}+a_{2\eta}\sin^2{\theta}}\f$
//!      (\f$ \xi_\mathrm{E}\f$ is the elliptical variable \f$ \xi(r,\theta)\f$ for \f$ r=R_{\mathrm{E}} \f$
//!
//!  \f$ \dfrac{df_0}{d\theta}= \left[\dfrac{d \psi_E(r,\theta)}{d\theta}\right]_{r=R_{\mathrm{E}}} = \dfrac{R^2_{\mathrm{E}}}{2}\dfrac{\alpha(\xi_{\mathrm{E}})}{\xi_{\mathrm{E}}}\mathcal{A}(\eta)\sin{2\theta},\quad \mathcal{A}(\eta)=a_{2\eta}-a_{1\eta} \f$
//!
/*!
  \param theta : angular coordinate in the lens plane
  \param pot_params[] : PNFW lens parameters
  \return \f$ \dfrac{df_0}{d\theta}\f$ */

double Df0Dtheta_nfw( double theta, double pot_params[]){
    double eta=pot_params[2];
    double a_eta1=1.0-eta,a_eta2=1.0/(1.0-eta); // If you want work with the standar parameterization
//     double a_eta1=1.0-eta,a_eta2=1.0+eta;    // // If you want work with the Angle Deflection Method
    double A_eta=a_eta2-a_eta1;
    double xie=_r_e*sqrt(a_eta1*pow(cos(theta),2)+a_eta2*pow(sin(theta),2));
    double cal_g=pow(_r_e,2)*A_eta*sin(2*theta);
   double df0dte=0.5*(alpha_nfw(xie,pot_params)/xie)*cal_g;
   return df0dte; 
}
/*********************************************************************************************************************/
//!  Solution for the PNFW model as perturbation,useful for critical and caustic lines
//!
//!  \f$ \psi_E(r)=\phi_E(\xi)-\phi_0(r) \f$
//!
//!  \f$ \xi=r\sqrt{a_{1\eta}\cos^2{\theta}+a_{2\eta}\sin^2{\theta}}\f$
//!                  (\f$ \xi_\mathrm{E}\f$ is the elliptical variable \f$ \xi(r,\theta)\f$ for \f$ r=R_{\mathrm{E}})\f$ 
//!
//!  \f$ \dfrac{d^2f_0}{d\theta^2}= \left[\dfrac{d^2 \psi_E(r,\theta)}{d\theta^2}\right]_{r=R_{\mathrm{E}}}= \mathcal{G}(\eta)\dfrac{\alpha(\xi_{\mathrm{E}})}{\xi_{\mathrm{E}}}-\mathcal{H}(\eta)\gamma(\xi_{\mathrm{E}})\f$
//! 
//!  \f$ \mathcal{G}(\eta)=\mathcal{A}(\eta)r^2_{\mathrm{E}}\cos{(2\theta)}\f$
//!
//!  \f$ \mathcal{H}(\eta)=\dfrac{r^2_{\mathrm{E}}}{2}\left\{\mathcal{A}(\eta)\dfrac{R_{\mathrm{E}}}{\xi_{\mathrm{E}}}\sin{(2\theta)}\right\}^2 \f$
//!
//! \f$ \gamma(\xi_{\mathrm{E}})=\kappa(\xi_{\mathrm{E}})-\dfrac{\alpha(\xi_{\mathrm{E}})}{\xi_{\mathrm{E}}} \f$
/*!
  \param theta : angular coordinate in the lens plane
  \param r : radial distance
  \param pot_params[] : PNFW lens parameters
  \return \f$ \dfrac{d^2f_0}{d\theta^2}\f$
*/

double D2f0Dtheta2_nfw( double theta, double pot_params[]){
    double eta=pot_params[2];
    double a_eta1=1.0-eta,a_eta2=1.0/(1.0-eta); // If you are interested in work with the standard parametrization
//     double a_eta1=1.0-eta,a_eta2=1.0+eta;    // If you are interested in work with the Angle Deflection Model
    double A_eta=a_eta2-a_eta1;
    double xie=_r_e*sqrt(a_eta1*pow(cos(theta),2)+a_eta2*pow(sin(theta),2));
    double cal_g=pow(_r_e,2)*A_eta*sin(2*theta);
    double cal_gsq=pow(cal_g,2)/_r_e;
    double dcal_g=2.0*pow(_r_e,2)*A_eta*cos(2*theta);
    
    double d2f0dte2=0.5*dcal_g*(alpha_nfw(xie,pot_params)/xie)-0.5*shear_nfw(xie,pot_params)*cal_gsq;
    return d2f0dte2;
}
/*********************************************************************************************************************/
//!  Function  useful to find the Einstein Radius (i.e. solving the equation \f$ \alpha(r)-r=0  \f$)
//!
//! \f$ f(r)=1-\dfrac{\alpha(r)}{r} \f$
//!
/*!
  \param r : radial distance
  \param pot_params[] : PNFW lens parameters
  \return \f$ f(r)\f$
*/
double lambda_t_nfw(double r, double pot_params[]) {return 1.0 - alpha_nfw(r, pot_params)/r;}

/*********************************************************************************************************************/
//!  Function useful to find the range  \f$ z_{\mathrm{min}} < R_{\mathrm{E}} < z_{\mathrm{max}}  \f$
//!
/*!
  \param f : function defining in lambda_t_nfw(double r, double pot_params[])
  \param pot_params[] : PNFW lens parameters
  \return \f$ z_{\mathrm{min}},z_{\mathrm{min}}\f$
*/

double bracketing_lambda_t(double f(double r, double params[]), double params[], double out[], double z=1.0, double dz=0.001){
  double z_min = 0.0;
  double z_max = 0.0;


  if(f(z, params) < 0){
    while( f(z, params) < 0){
      z += dz;
    }
    z_min = z-dz;
    z_max = z;
  } else {
    while(f(z, params) > 0){
      z -= dz;
    }
    z_min = z;
    z_max = z+dz;
  }

  out[0] = z_min;
  out[1] = z_max;

  return 0.0;
}

// Corrected by Habib Dumet, (I change conv by shear)
//!  Function  using the Newton-Raphson method to Einstein Radius \f$ R_{\mathrm{E}} \f$
//!
//! We define some guess value \f$ r_0 \f$ as \f$ r_0=\left(\dfrac{z_{\mathrm{min}}+z_{\mathrm{min}}}{2} \right)\f$. Therefore, following the Newton-Raphson method the Einstein Radius is defined by
//!
//! \f$ R_{\mathrm{E}}=r_0-\dfrac{r_0-\alpha(r_0)}{2\gamma(r_0)}\f$, when \f$ f(R_{\mathrm{E}})<\textrm{tol}_1=10^{-7} \f$ or \f$ |R_{\mathrm{E}}-r_0|<\textrm{tol}_2=10^{-6}\f$  (the function \f$ f(r) \f$ is defined in lambda_t_nfw(double r, double pot_params[]))  
//! 
/*!
  \param alpha : the angle deflection of the NFW model
   \param gamma : the shear of the NFW model
  \param pot_params[] : PNFW lens parameters
  \param z_min : minimum value for the range, defined at bracketing_lambda_t
  \param z_max : maximum value for the range defined at bracketing_lambda_t 
  \return \f$ R_{\mathrm{E}}\f$
*/

double root_find(double alpha(double r, double params[]), double gamma(double r, double params[]), double params[], double z_min, double z_max){

  double z0 = (z_min + z_max)/2.0;

  double z1 = z0 - (z0 - alpha(z0,params))/(2.0*gamma(z0,params));

  while( fabs(z0-z1) > 1E-6 || fabs((1.0 - alpha_nfw(z1, params)/z1)) > 1E-7 ){
    z0=z1;
    z1 = z0 - (z0 - alpha(z0,params))/(2.0*gamma(z0,params));
  }

  return z1;
}

//!  Einstein Radius of the NFW lens model \f$ R_{\mathrm{E}} \f$
//!
//! This function basically call the bracketing_lambda_t and root_find functions.
//!
/*! \param pot_params[] : PNFW lens parameters, \return \f$ R_{\mathrm{E}} \f$ */

double r_e_nfw(double pot_params[]){

  double *out = (double*)malloc(2.0*sizeof(double));
  bracketing_lambda_t( lambda_t_nfw, pot_params, out);

  return root_find(alpha_nfw, conv_nfw, pot_params,out[0] ,out[1]);
}



#endif
