/** @file
* Example of doxygen documentation for C functions FIXME. 
*/

/** @package enfw_model
*  Package to compute quantities related to ENFW model.
*
*  Detailed descrition FIXME
*
*/
#ifndef ENFW_MODEL_H
#define ENFW_MODEL_H

#include "elliptical_mass.h"

//!  First perturbative field for the ENFW model as Perturbation.
//!
//!  
/*!
  \param theta : angular coordinate in the lens plane,
  \param pot_params[0]: ellipticity 
  \param pot_params[1]: flag control: 1)  \f$ a=1/\sqrt{1-\varepsilon} , b=1/\sqrt{1+\varepsilon}\f$ , 2)\f$ a=1/\sqrt{1-\varepsilon},b=\sqrt{1-\varepsilon}\f$, 3) \f$ a=1 , b=1-\varepsilon\f$
  \param pot_params[2]: \f$\kappa_s\f$
  \param pot_params[3]: \f$r_s\f$
  \return \f$ f_1(\theta)\f$
*/
double f1_enfw(double theta, double pert_params[],double r_e){
  double conv_params[]={pert_params[2],pert_params[3]};
  double eps=pert_params[0];
  double a_1eps, a_2eps;
  int iflag=pert_params[1];

  if(iflag==1){
    // printf("You chose the Angle Deflection Parameterization: a_1eps=sqrt(1.-eps), a_2eps=sqrt(1.+eps)");
    a_1eps=1.0/sqrt(1.0-eps);
    a_2eps=1.0/sqrt(1.0+eps);
  }

  if(iflag==2){
    // printf("You chose the standard parameterization: a_1eps=1.0/sqrt(1.0-eps),a_2eps=sqrt(1.0-eps)");
    a_1eps=1.0/sqrt(1.0-eps);
    a_2eps=sqrt(1.0-eps);
  }

  if(iflag==3){
    // printf("You chose the Keeton's parameterization: a_1eps=1.0, a_2eps=1./sqrt(1.-eps)");
    a_1eps=1.0;
    a_2eps=(1.0-eps);
  }

  return  f_1_ellip_mass(conv_nfw_circ,theta, conv_params, r_e, a_1eps, a_2eps);
}

//!  Second perturbative field for the ENFW model as Perturbation.
//!
//!  
/*!
  \param theta : angular coordinate in the lens plane,
  \param pot_params[0]: ellipticity 
  \param pot_params[1]: flag control: 1)  \f$ a=1/\sqrt{1-\varepsilon} , b=1/\sqrt{1+\varepsilon}\f$ , 2)\f$ a=1/\sqrt{1-\varepsilon},b=\sqrt{1-\varepsilon}\f$, 3) \f$ a=1 , b=1-\varepsilon\f$
  \param pot_params[2]: \f$\kappa_s\f$
  \param pot_params[3]: \f$r_s\f$
  \return \f$ \frac{df_0}{d\theta}\f$
*/
double Df0Dtheta_enfw(double theta, double pert_params[],double r_e){
  double conv_params[]={pert_params[2],pert_params[3]};
  double eps=pert_params[0];
  double a_1eps, a_2eps;
  int iflag=pert_params[1];

  if(iflag==1){
    // printf("You chose the Angle Deflection Parameterization: a_1eps=1.-eps, a_2eps=1.+eps");
    a_1eps=sqrt(1.-eps), a_2eps=sqrt(1.+eps);
  }

  if(iflag==2){
    // printf("You chose the standard parameterization: a_1eps=1.-eps,a_2eps=1./a_1eps");
    a_1eps=1.0/sqrt(1.0-eps);
    a_2eps=sqrt(1.0-eps);
  }

  if(iflag==3){
    // printf("You chose the Keeton's parameterization: a_1eps=1.0, a_2eps=1./((1.-eps)^2)");
    a_1eps=1.0, a_2eps=1./(1.-eps);
  }

  return  Df0Dtheta_ellip_mass(conv_nfw_circ,theta, conv_params, r_e, a_1eps, a_2eps);
}

//!  Solution for the ENFW model as perturbation,useful for critical and caustic lines
//!
//!  
/*!
  \param theta : angular coordinate in the lens plane,
  \param pot_params[0]: ellipticity 
  \param pot_params[1]: flag control: 1)  \f$ a=1/\sqrt{1-\varepsilon} , b=1/\sqrt{1+\varepsilon}\f$ , 2)\f$ a=1/\sqrt{1-\varepsilon},b=\sqrt{1-\varepsilon}\f$, 3) \f$ a=1 , b=1-\varepsilon\f$
  \param pot_params[2]: \f$\kappa_s\f$
  \param pot_params[3]: \f$r_s\f$
  \return  \f$ \frac{d^2f_0}{d\theta^2}\f$
*/
double D2f0Dtheta2_enfw(double theta, double pert_params[],double r_e){
  double conv_params[]={pert_params[2],pert_params[3]};
  double eps=pert_params[0];
  double a_1eps, a_2eps;
  int iflag=pert_params[1];

  if(iflag==1){
    // printf("You chose the Angle Deflection Parameterization: a_1eps=1.-eps, a_2eps=1.+eps");
    a_1eps=sqrt(1.-eps), a_2eps=sqrt(1.+eps);
  }

  if(iflag==2){
    // printf("You chose the standard parameterization: a_1eps=1.-eps,a_2eps=1./a_1eps");
    a_1eps=1.0/sqrt(1.0-eps);
    a_2eps=sqrt(1.0-eps);
  }

  if(iflag==3){
    // printf("You chose the Keeton's parameterization: a_1eps=1.0, a_2eps=1./((1.-eps)^2)");
    a_1eps=1.0, a_2eps=1./(1.-eps);
  }

  return  D2f0Dtheta2_ellip_mass(conv_nfw_circ, conv_nfw_circ_prime, theta, conv_params, r_e, a_1eps, a_2eps);
}



#endif
