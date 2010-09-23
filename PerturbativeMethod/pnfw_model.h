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
// Added by Habib Dumet
double shear_nfw(double X, double pot_params[]);

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
    double a_eta1=1.0-eta,a_eta2=1.0/(1.0-eta);
//     double a_eta1=1.0-eta,a_eta2=1.0+eta;
    double reta=_r_e*sqrt(a_eta1*cos(theta)*cos(theta)+a_eta2*sin(theta)*sin(theta));
    double f1=(reta/_r_e)*alpha_nfw(reta,pot_params)-alpha_nfw(_r_e,pot_params);
    return f1;
}

double Df0Dtheta_nfw( double theta, double pot_params[]){
    double eta=pot_params[2];
    double a_eta1=1.0-eta,a_eta2=1.0/(1.0-eta);
//     double a_eta1=1.0-eta,a_eta2=1.0+eta;
    double A_eta=a_eta2-a_eta1;
    double reta=_r_e*sqrt(a_eta1*cos(theta)*cos(theta)+a_eta2*sin(theta)*sin(theta));

   double df0dte=(_r_e/2)*alpha_nfw(reta,pot_params)*(_r_e/reta)*A_eta*sin(2*theta);
//    return (df0dte/_r_e);
   return df0dte; 
}

double D2f0Dtheta2_nfw( double theta, double pot_params[]){
    double eta=pot_params[2];
    double a_eta1=1.0-eta,a_eta2=1.0/(1.0-eta);
//     double a_eta1=1.0-eta,a_eta2=1.0+eta;
    double A_eta=a_eta2-a_eta1;
    double reta=_r_e*sqrt(a_eta1*cos(theta)*cos(theta)+a_eta2*sin(theta)*sin(theta));
    double cal_g=A_eta*pow(_r_e,2)*cos(2.0*theta);
    double cal_h=cal_g-(pow(_r_e,2)/2.0)*pow(A_eta*(_r_e/reta)*sin(2.0*theta),2);
    
    double d2f0dte2=cal_g*conv_nfw(reta,pot_params)+cal_h*shear_nfw(reta,pot_params);
    return d2f0dte2;
}

// pot_paramas[0]=r_s, pot_params[1]=ks,pot_params[2]=elipticity, pot_params[3]

double kappa2_nfw(double pot_params[]){
  double k2,X=_r_e;

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


double alpha_nfw(double X, double pot_params[])
{
    //double X = r/pot_params[0];
    double alpha;
    //X=r_e/pot_params[0]

      alpha = 4*pot_params[1]*(log(X/2)+F(X))/X;

      return alpha;
}

// Added by Habib Dumet-Montoya
double shear_nfw(double X, double pot_params[]){
      double gamma;

      gamma=alpha_nfw(X, pot_params)/X-conv_nfw(X, pot_params);

      return gamma;
}

double lambda_t_nfw(double X, double pot_params[]) {return 1.0 - alpha_nfw(X, pot_params)/X;}

double bracketing_lambda_t(double f(double X, double params[]), double params[], double out[], double z=1.0, double dz=0.001){
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
  //printf("%f %f\n",z_min, f(z_min, params));
  //printf("%f %f\n",z_max, f(z_max, params));

  return 0.0;
}

// Corrected by Habib Dumet, (I change conv by shear)
double root_find(double alpha(double X, double params[]), double gamma(double X, double params[]), double params[], double z_min, double z_max){

  double z0 = (z_min + z_max)/2.0;

  double z1 = z0 - (z0 - alpha(z0,params))/(2.0*gamma(z0,params));

  while( fabs(z0-z1) > 1E-6 || fabs((1.0 - alpha_nfw(z1, params)/z1)) > 1E-7 ){
    z0=z1;
    z1 = z0 - (z0 - alpha(z0,params))/(2.0*gamma(z0,params));
  }
  //printf("%.10f %E\n",z1,1.0 - alpha_nfw(z1, params)/z1);


  return z1*params[0];
}

double r_e_nfw(double pot_params[]){

  double *out = (double*)malloc(2.0*sizeof(double));
  bracketing_lambda_t( lambda_t_nfw, pot_params, out);

  return root_find(alpha_nfw, conv_nfw, pot_params,out[0] ,out[1]);
}



#endif
