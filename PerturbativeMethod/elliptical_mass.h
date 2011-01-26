#ifndef ELLIPTICAL_MASS_H
#define ELLIPTICAL_MASS_H

#include <cmath>

#include "pot_derivatives.h"

// g++ -Wall ellip_mass.cpp `pkg-config gsl --cflags --libs`

double f_1_ellip_mass(conv_type conv_in,double theta, double conv_params[], double r_e, double a, double b){
  double X1 = r_e*cos(theta);
  double X2 = r_e*sin(theta);

  double J0 = J(0, conv_in, X1, X2, conv_params, a, b);
  double J1 = J(1, conv_in, X1, X2, conv_params, a, b);

  double out = r_e* ( a*b*( J0*pow(cos(theta),2) + J1*pow(sin(theta),2)) );

  return out;
}


double Df0Dtheta_ellip_mass(conv_type conv_in,double theta, double conv_params[], double r_e, double a, double b){
  double X1 = r_e*cos(theta);
  double X2 = r_e*sin(theta);

  double J0 = J(0, conv_in, X1, X2, conv_params, a, b);
  double J1 = J(1, conv_in, X1, X2, conv_params, a, b);

  double out = (a*b*r_e*r_e/2.0) * (J1 - J0)*sin(2.0*theta);

  return out;
}


double D2f0Dtheta2_ellip_mass(conv_type conv_in,double theta, double conv_params[], double r_e, double a, double b){
  double X1 = r_e*cos(theta);
  double X2 = r_e*sin(theta);

  double J0 = J(0, conv_in, X1, X2, conv_params, a, b);
  double J1 = J(1, conv_in, X1, X2, conv_params, a, b);

  double K0 = K(0, conv_in, X1, X2, conv_params, a, b);
  double K1 = K(1, conv_in, X1, X2, conv_params, a, b);
  double K2 = K(2, conv_in, X1, X2, conv_params, a, b);

  double out1 = a*b*r_e*r_e;
  double out2 = (J1 - J0)*cos(2.0*theta);
  double out3 = r_e*r_e*(K0 - 2.0*K1 + K2)*pow(sin(2.0*theta),2)/2.0;

  return out1*(out2 + out3);
}
#endif
