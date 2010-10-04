/** @file
* Example of doxygen documentation for C functions FIXME. 
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

double _r_e_nfw = 1.0;

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



//! Convergence of the NFW model.
//!
//! \f$ \kappa(r)=2\,\kappa_s\dfrac{(1-\mathcal{F}(x))}{x^2-1}\f$, 
//!
//! where \f$ \kappa_s \f$ is the characteristic convergence and \f$ \mathcal{F}(x)\f$ is defined in F(double x)
/*!   \param r : radial distance,   \param pot_params[] : NFW lens parameters,   \return \f$ \kappa(r)\f$ */
double conv_nfw_circ(double r, double pot_params[])
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
    double alpha,rs=pot_params[0],ks=pot_params[1];
    double X=r/rs;
      alpha = 4.0*ks*rs*(log(X/2.0)+F(X))/X;
      return alpha;
}


//! Shear of the NFW model.
//!
//! \f$ \gamma(r)= \dfrac{\alpha(r)}{r}-\kappa(r) \f$, 
//!
//!  where \f$ \alpha(r) \f$ is angle deflection and \f$ \kappa(r) \f$ is the convergence.
//!
/*!   \param r : radial distance,   \param pot_params[] : NFW lens parameters, \return \f$ \gamma(r)\f$ */
double shear_nfw_circ(double r, double pot_params[]){
      double gamma;
      //double rs=pot_params[0];
      //double X=r/rs;
      gamma=alpha_nfw_circ(r, pot_params)/r-conv_nfw_circ(r, pot_params);
      return gamma;
}


//! Solution for the Circular NFW Lens 
//!
//! \f$ \kappa_2=2-2\kappa(R_{\mathrm{E}}) \f$,
//!
//! where \f$ \kappa(R_{\mathrm{E}})\f$ is the convergence for evaluated at the Einstein Radius
//!
/*!   \param r_e : Einstein Radius, \param pot_params[] : NFW lens parameters,\return \f$ \kappa_2 \f$ */
double kappa2_nfw(double pot_params[]){
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
  \param pot_params[] : PNFW lens parameters
  \return \f$ f(r)\f$
*/
double lambda_t_nfw_circ(double r, double pot_params[]) {return 1.0 - alpha_nfw_circ(r, pot_params)/r;}


//pot_params[0] = rs
//pot_params[1] = kappas
double re_find_func_nfw(double r, void *p){

  double *params = (double *)p;

  double rs = params[0];//(params->rs);
  double ks = params[1];//(params->ks);
  //printf("%E %E\n",rs,ks);
  double X = r/rs;

  return r - 4.0*ks*rs/X * ( F(X) + log(X/2.0) );

}







#endif
