/** @file
* Package useful to calculate some function related to the NFW lensing model
*
* Lensing Function (Angle Deflection, Convergence, Shear)
*
* Einstein Radius (Use bracketing_lambda_t and gsl functions)
*
* Function useful to the Perturbative Approach \f$ \kappa_2 \f$
*
*/

/** @package nfw_circular_model
*  Package to compute quantities related to the circular NFW model.
*
*  Detailed descrition FIXME
*
*/
#ifndef NFW_CIRCULAR_MODEL_H
#define NFW_CIRCULAR_MODEL_H

#include <math.h>
#include <cstdlib>

#include "../numerical_methods/general_methods.h"
//#include "general_methods.h"

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


//pot_params[0] = kappas
//pot_params[1] = rs
//! Convergence of the NFW model.
//!
//! \f$ \kappa(r)=2\,\kappa_s\dfrac{(1-\mathcal{F}(x))}{x^2-1}\f$, 
//!
//! where \f$ \kappa_s \f$ is the characteristic convergence and \f$ \mathcal{F}(x)\f$ is defined in F(double x)
/*!   \param r : radial distance,   \param pot_params[] : NFW lens parameters,   \return \f$ \kappa(r)\f$ */
double conv_nfw_circ(double r, double pot_params[])
{
    double K, ks=pot_params[0];
    double rs=pot_params[1],X=r/rs;
    
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

  out[0] = z_min;
  out[1] = z_max;

  return 0.0;
}




//pot_params[0] = kappas
//pot_params[1] = rs
double r_e_nfw_find(double pot_params[], double *est_err_out=NULL, double x_lo = 1E-4, double x_hi = 10.0, int max_iter = 100, double relative_error = 1E-4, int v=1){
  double *params = (double*) malloc(2.0*sizeof(double));
  params[0] = pot_params[0];
  params[1] = pot_params[1];
  
  double out[2];
  
  bracketing_lambda_t(re_find_func_nfw, params, out);
  
  
  double re = root_find(re_find_func_nfw, params,est_err_out, out[0],out[1],max_iter,relative_error,v);
  return re;
}



#endif