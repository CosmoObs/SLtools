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

double conv_nfw(double X, double pot_params[]);
double alpha_nfw(double X, double pot_params[]);

//to compile: g++ -Wall nfw_pert.cpp
/*********************************************************************************************************************/
//!  Function defined in PerturbedPNFW.pdf 
//!
//!  \f$ f_{n}(\theta) \equiv \frac{1}{n!} \left[\frac{\partial\psi(r,\theta)}{\partial r}\right]_{r=1} \f$
/*!
  \param theta angular coordinate in the lens plane
  \param pot_params[] potential parameters
  \return 0.0
*/
double f1_nfw(double theta, double pot_params[]){
    double eta=pot_params[2];
    double a_eta1=1-eta,a_eta2=1/(1-eta);
    double reta=_r_e*sqrt(a_eta1*cos(theta)*cos(theta)+a_eta2*sin(theta)*sin(theta));
    double f1=(reta/_r_e)*alpha_nfw(reta,pot_params)-alpha_nfw(_r_e,pot_params);
    return f1;
}

double Df0Dtheta_nfw( double theta, double pot_params[]){
    double eta=pot_params[2];
    double a_eta1=1-eta,a_eta2=1/(1-eta);
    double A_eta=a_eta2-a_eta1;
    double reta=_r_e*sqrt(a_eta1*cos(theta)*cos(theta)+a_eta2*sin(theta)*sin(theta));

   double df0dte=(_r_e/2)*alpha_nfw(reta,pot_params)*(_r_e/reta)*A_eta*sin(2*theta);
   return (df0dte/_r_e);
}


// pot_paramas[0]=r_s, pot_params[1]=ks,pot_params[2]=elipticity, pot_params[3]

double kappa2_nfw(double pot_params[]){
  double k2,X=_r_e/pot_params[0];

  k2=2-2*conv_nfw(X,pot_params);

  return k2;
}


/*******************************************************************************
 Lensing Functions related to the NFW model
*******************************************************************************/
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

/******************************************************************************/
double conv_nfw(double X, double pot_params[])
{
    double K, ks=pot_params[1];
    
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


double alpha_nfw(double r, double pot_params[])
{
    double X = r/pot_params[0];
    double alpha;
    //X=r_e/pot_params[0]

      alpha = 4*pot_params[0]*pot_params[1]*(log(X/2)+F(X))/X;

      return alpha;
}




#endif
