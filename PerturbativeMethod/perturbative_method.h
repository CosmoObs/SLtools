#ifndef PERTURBATIVEMETHOD_H
#define PERTURBATIVEMETHOD_H


#include <cmath> 

/*********************************************************************************************************************/
typedef double (*phi_type) (double r, double pot_params[]);
typedef double (*psi_type) (double r, double theta, double pot_params[]);
typedef double (*f_type)   (double theta, double pot_params[]);
/*********************************************************************************************************************/

typedef struct {
  double x0 ;     //position at x axis
  double y0 ;     //position at y axis
  double R0 ;     //characteristic size
  double eta0 ;   //parameter relater to the ellipticity ( e = sqrt(2*eta0) )
  double theta0 ; //inclination to the main axis
} elliptical_source;




double f1_bar(f_type f1_in, elliptical_source source_in, double theta, double pot_params[]){
  return f1_in(theta, pot_params) + source_in.x0*cos(theta) + source_in.y0*sin(theta);
}

double Df0Dtheta_bar(f_type Df0Dtheta_in, elliptical_source source_in, double theta, double pot_params[]){
  return Df0Dtheta_in(theta, pot_params) - source_in.x0*sin(theta) +  source_in.y0*cos(theta);
}

double arg_sqrt(f_type Df0Dtheta_in, elliptical_source source_in, double theta, double pot_params[]){
  double Df0Dtheta_bar_tmp = Df0Dtheta_bar(Df0Dtheta_in, source_in, theta, pot_params);
  double Df0Dtheta_bar_tmp2 = Df0Dtheta_bar_tmp*Df0Dtheta_bar_tmp;
  double theta_bar = theta - source_in.theta0;
  double S = 1.0 - source_in.eta0*cos(2.0*theta_bar);

  double R0_tmp2 = source_in.R0*source_in.R0;
  double eta0_tmp2 = source_in.eta0*source_in.eta0;

  return R0_tmp2*S - (1.0-eta0_tmp2)*Df0Dtheta_bar_tmp2;
}

double dr_plus(f_type f1_in, f_type Df0Dtheta_in, double kappa2, elliptical_source source_in, double theta, double pot_params[]){
  double Df0Dtheta_bar_tmp = Df0Dtheta_bar(Df0Dtheta_in, source_in, theta, pot_params);
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

double dr_minus(f_type f1_in, f_type Df0Dtheta_in, double kappa2, elliptical_source source_in, double theta, double pot_params[]){
  double Df0Dtheta_bar_tmp = Df0Dtheta_bar(Df0Dtheta_in, source_in, theta, pot_params);
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


#endif
