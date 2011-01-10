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


#endif
