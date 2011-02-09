/** @file
* Example of doxygen documentation for C functions FIXME. 
*/

/** @package elliptical_mass
*  Package to compute quantities related to the elliptical mass models
*
*  Detailed descrition FIXME
*
*/

#ifndef ELLIPTICAL_MASS_H
#define ELLIPTICAL_MASS_H

#include <cmath>

#include "pot_derivatives.h"
#include "nfw_circular_model.h"

// g++ -Wall ellip_mass.cpp `pkg-config gsl --cflags --libs`

/**  function to compute the \f$ f_1 \f$ for elliptical models
*
*  \param conv_type conv_in pointer to the convergencie for the circular model
*  \param theta angular coordinate in the lens plane
*  \param conv_params[] vector with the parameters of the circular convergence
*  \param r_e is the einstein radius (\f$ R_{_{\mathrm{E}}}\f$ )
*  \param a parameter related to ellipticity (see eq 3.30 from report, section General results for elliptical models)
*  \param b parameter related to ellipticity (see eq 3.30 from report, section General results for elliptical models)
*  \return \f$ f_1 \f$ for a elliptical model
*
*/
double f_1_ellip_mass(conv_type conv_in,double theta, double conv_params[], double r_e, double a, double b){
  double X1 = r_e*cos(theta);
  double X2 = r_e*sin(theta);


  double J0 = J(0, conv_in, X1, X2, conv_params, a, b);
  double J1 = J(1, conv_in, X1, X2, conv_params, a, b);

  double out = r_e* ( a*b*( J0*pow(cos(theta),2) + J1*pow(sin(theta),2)) );

  return out;
}

/**  function to compute the \f$ \frac{df_0}{d\theta} \f$ for elliptical models
*
*  \param conv_type conv_in pointer to the convergencie for the circular model
*  \param theta angular coordinate in the lens plane
*  \param conv_params[] vector with the parameters of the circular convergence
*  \param r_e is the einstein radius (\f$ R_{_{\mathrm{E}}}\f$ )
*  \param a parameter related to ellipticity (see eq 3.30 from report, section General results for elliptical models)
*  \param b parameter related to ellipticity (see eq 3.30 from report, section General results for elliptical models)
*  \return \f$ \frac{df_0}{d\theta} \f$ for a elliptical model
*
*/
double Df0Dtheta_ellip_mass(conv_type conv_in,double theta, double conv_params[], double r_e, double a, double b){
  double X1 = r_e*cos(theta);
  double X2 = r_e*sin(theta);

  double J0 = J(0, conv_in, X1, X2, conv_params, a, b);
  double J1 = J(1, conv_in, X1, X2, conv_params, a, b);

  double out = (a*b*r_e*r_e/2.0) * (J1 - J0)*sin(2.0*theta);

  return out;
}

/**  function to compute the \f$ \frac{d^2f_0}{d\theta^2} \f$ for elliptical models
*
*  \param conv_type conv_in pointer to the convergencie for the circular model
*  \param theta angular coordinate in the lens plane
*  \param conv_params[] vector with the parameters of the circular convergence
*  \param r_e is the einstein radius (\f$ R_{_{\mathrm{E}}}\f$ )
*  \param a parameter related to ellipticity (see eq 3.30 from report, section General results for elliptical models)
*  \param b parameter related to ellipticity (see eq 3.30 from report, section General results for elliptical models)
*  \return \f$ \frac{d^2f_0}{d\theta^2} \f$ for a elliptical model
*
*/
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




double nfw_ellip_f1(double theta, double pert_params[],double r_e){
  double conv_params[]={pert_params[2],pert_params[3]};
  double eps=pert_params[0], a_1eps,a_2eps;
  int iflag=pert_params[1];

  if(iflag==1){
    // printf("You chose the Angle Deflection Parameterization");
    a_1eps=1.-eps, a_2eps=1.+eps;
  }

  if(iflag==2){
    // printf("You chose the standard parameterization");
    a_1eps=1.-eps,a_2eps=1./a_1eps;
  }

  if(iflag==3){
    // printf("You chose the Keeton's parameterization");
    a_1eps=1.0, a_2eps=1./pow((1.-eps),2);
  }

  return  f_1_ellip_mass(conv_nfw_circ,theta, conv_params, r_e, a_1eps, a_2eps);

}

#endif
