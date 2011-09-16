/** @file
* Package to calcule quantities related to the NFW model
*
* - Lensing Function (Angle Deflection, Convergence, Shear)
*
* - Einstein Radius (Use bracketing_lambda_t and gsl functions)
*
* - Function useful to the Perturbative Approach \f$ \kappa_2 \f$
*
* - Function related with the third derivative of the lensing potential \f$ c_3 \f$ 
* 
*  This module need 
* - pot_params[]={\f$\kappa_s, r_s \f$}
*/

/** @package circular_models
*  Package to compute quantities related to the circular NFW model.
*
*
*/
#ifndef NFW_CIRCULAR_MODEL_H
#define NFW_CIRCULAR_MODEL_H

#include <math.h>
#include <cstdlib>

#include "../numerical_methods/root_find.h"
//#include "general_methods.h"

//! Function used to define the convergence and angle deflection of the NFW model.
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

//! Derivative of F function
/*!
  \param X radial coordinate in units of r_s (r_s is the NFW profile parameter)
  \return \f$ \frac{dF(x)}{dx} = \frac{1-x^2 F(x)}{x(x^2-1)}\f$
  \sa F()
*/
double Flinha(double X)
{
    double x2;
    double arg;
    double flinha;
    
    x2 = X*X;
    arg = x2 - 1.0;
    
    flinha = (1.0 - x2*F(X))/(X*arg);
    
    
    return flinha;
}

//pot_params[0] = kappas
//pot_params[1] = rs
//! Convergence of the NFW model.
//!
//! \f$ \kappa(r)=2\,\kappa_s\dfrac{(1-\mathcal{F}(x))}{x^2-1}\f$, 
//!
//! where \f$ \kappa_s \f$ is the characteristic convergence and \f$ \mathcal{F}(x)\f$ is defined in F(double x)
/*!   \param r : radial distance,
      \param pot_params[] : NFW lens parameters,
      \return \f$ \kappa(r)\f$
*/
double conv_nfw_circ(double r, double pot_params[])
{
    double K, ks=pot_params[0];
    double rs=pot_params[1];
    double X=r/rs;
    
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

//! Derivative of conv_nfw_circ(
/*!
  \param r radial distance
  \param pot_params[0] : NFW lens parameters kappa_s,
  \param pot_params[1] : NFW lens parameters r_s,
  \return \f$ \frac{dk(x)}{d\left(x^2\right)}\f$
  \sa F()
*/
double conv_nfw_circ_prime(double r, double pot_params[])
{

  double ks=pot_params[0];
  double rs=pot_params[1];
  double X=r/rs;
  double X2=X*X;
  double UMX2=X2-1.0;
	
	if( fabs(X-1.0) < 1e-8) return -0.4*ks;
    return ks*( 2.0*(F(fabs(X))-1.0) - UMX2*Flinha(fabs(X))/X )/(UMX2*UMX2);
}

//pot_params[0] = kappas
//pot_params[1] = rs
//! Angle deflection of the NFW model.
//!
//! \f$ \alpha(r)=   2\,\kappa_s r_s \dfrac{\log{x/2}+\mathcal{F}(x)}{x}\f$, 
//!  where \f$ \kappa_s \f$ is the characteristic convergence, \f$ r_s \f$ is the scale radius, \f$ x=r/r_s \f$  and \f$ \mathcal{F}(x)\f$ is defined in F(double x)
/*!
  \param r : radial distance,\param pot_params[] : NFW lens parameters 
  \return \f$ \alpha(r)\f$
*/
double alpha_nfw_circ(double r, double pot_params[])
{
    double alpha,rs=pot_params[1],ks=pot_params[0];
    double X=r/rs;
      alpha = 4.0*ks*rs*(log(X/2.0)+F(X))/X;
      return alpha;
}

//pot_params[0] = kappas
//pot_params[1] = rs
//! Shear of the NFW model.
//!
//! \f$ \gamma(r)= \dfrac{\alpha(r)}{r}-\kappa(r) \f$, 
//!
//!  where \f$ \alpha(r) \f$ is angle deflection and \f$ \kappa(r) \f$ is the convergence.
//!
/*!   \param r : radial distance,   \param pot_params[] : NFW lens parameters, \return \f$ \gamma(r)\f$ */
      double shear_nfw_circ(double r, double pot_params[]){
      double gamma;
      gamma=alpha_nfw_circ(r, pot_params)/r-conv_nfw_circ(r, pot_params);
      return gamma;
}


//pot_params[0] = kappas
//pot_params[1] = rs
//! Solution for the Circular NFW Lens 
//!
//! \f$ \kappa_2=2-2\kappa(R_{\mathrm{E}}) \f$,
//!
//! where \f$ \kappa(R_{\mathrm{E}})\f$ is the convergence for evaluated at the Einstein Radius
//!
/*!   \param r_e : Einstein Radius, \param pot_params[] : NFW lens parameters,\return \f$ \kappa_2 \f$ */
double kappa2_nfw(double pot_params[],double _r_e_nfw){
  double k2,r=_r_e_nfw;

  k2=2.0-2.0*conv_nfw_circ(r,pot_params);

  return k2;
}
//! Third coefficient of the Taylor Expansion of the NFW potential
//!
//! \f$ C_3(r)= \frac{4}{r_s(X^2_{\mathrm{E}}-1)}\left[\frac{1}{X_{\mathrm{E}}}-\frac{3}{2}\frac{X_{\mathrm{E}}\kappa(R_{\mathrm{E}})}{2\kappa_s}  \right] +2 \frac{\gamma(R_{\mathrm{E}})}{R_{\mathrm{E}}} \f$
//!
/*! \param r_e_nfw : Einstein Radius of the NFW model, \param pot_params[] : NFW lens parameters,\return \f$ C_3 \f$ */

double c3_nfw(double pot_params[],double _r_e_nfw){
  double ks=pot_params[0],rs=pot_params[1];
  double Re=_r_e_nfw, Xe=Re/rs;
  double dkdr=(2.*ks/(rs*(pow(Xe,2)-1.0)))*((1.0/Xe)-(3./(2.*ks))*Xe*conv_nfw_circ(Re,pot_params));
//   printf("the value of the derivative of kappa at RE is %f\n", dkdr);
  double gor=shear_nfw_circ(Re, pot_params)/Re;
//   printf("the value of the shear/r is %f\n", gor);
  return 2.0*(dkdr+gor);
}


//!  Function  useful to find the Einstein Radius (i.e. solving the equation \f$ \alpha(r)-r=0  \f$)
//!
//! \f$ f(r)=1-\dfrac{\alpha(r)}{r} \f$
//!
/*!
  \param r : radial distance
  \param pot_params[] : NFW lens parameters
  \return \f$ f(r)\f$
*/
double lambda_t_nfw_circ(double r, double pot_params[]) {return 1.0 - alpha_nfw_circ(r, pot_params)/r;}



//pot_params[0] = kappas
//pot_params[1] = rs
double re_find_func_nfw(double r, void *p){

  double *params = (double *)p;

  double rs = params[1];//(params->rs);
  double ks = params[0];//(params->ks);
  //printf("%E %E\n",rs,ks);
  double X = r/rs;

  return r - 4.0*ks*rs/X * ( F(X) + log(X/2.0) );
}
double re_find_func_nfw(double r, double params[]){

  double rs = params[1];//(params->rs);
  double ks = params[0];//(params->ks);
  //printf("%E %E\n",rs,ks);
  double X = r/rs;

  return r - 4.0*ks*rs/X * ( F(X) + log(X/2.0) );
}

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

  /*if(z_min <1.E-8){*/
  z_min+=1E-8;

  out[0] = z_min;
  out[1] = z_max;
  

  return 0.0;
}

//pot_params[0] = kappas
//pot_params[1] = rs
double r_e_nfw_find(double pot_params[], double *est_err_out=NULL, double x_lo = 1E-4, double x_hi = 10.0, int max_iter = 100, double relative_error = 1E-6, int v=1){
  double *params = (double*) malloc(2.0*sizeof(double));
  params[0] = pot_params[0];
  params[1] = pot_params[1];
  
  double out[2];
  double v_in,dz;
  double ks=pot_params[0],rs=pot_params[1];
  if(ks<=0.11){v_in=5.E-3*rs;}
  if((ks>0.11) && (ks<=0.5)){v_in=2.5E-1*rs;}
  if((ks>0.15) && (ks<=1.0)){v_in=0.5*rs;}
  if(ks>1.0){v_in=1.0*rs;}
  dz=v_in*1.E-3;
  bracketing_lambda_t(re_find_func_nfw, params, out,v_in,dz);
  if(v) printf("the range where the root is [%f , %f]\n", out[0],out[1]);
    
  double re = root_find(re_find_func_nfw, params,est_err_out, out[0],out[1],max_iter,relative_error,v);
  return re;
}



#endif
