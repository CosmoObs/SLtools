#ifndef ELLIPTICAL_MASS_H
#define ELLIPTICAL_MASS_H

#include <cmath>

#include "pot_derivatives.h"

// g++ -Wall ellip_mass.cpp `pkg-config gsl --cflags --libs`

<<<<<<< HEAD
double f_1_ellip_mass(conv_type conv_in,double X1, double X2, double conv_params[], double r_e, double a, double b){
  double theta = atan(X2/X1);
  //double r = sqrt(X1*X1 + X2*X2);
=======
double f_1_ellip_mass(conv_type conv_in,double theta, double conv_params[], double r_e, double a, double b){
  double X1 = r_e*cos(theta);
  double X2 = r_e*sin(theta);
>>>>>>> 3a00247c50eb340b350d432b579539c5c6d10cb8

  double J0 = J(0, conv_in, X1, X2, conv_params, a, b);
  double J1 = J(1, conv_in, X1, X2, conv_params, a, b);

  double out = r_e* ( a*b*( J0*pow(cos(theta),2) + J1*pow(sin(theta),2)) );

  return out;
}


#endif
