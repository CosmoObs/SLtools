/** @file
* Example of doxygen documentation for C function FIXME
*/
/** @package sie_model
*  Package to compute quantities related to Pseudo-Elliptical SIS model.
*
*  Detailed descrition FIXME
*
*/

#include "sis_model.h"
#include "elliptical_mass.h"

#include "cmath"


#ifndef SIE_MODEL_H
#define SIE_MODEL_H

double J0_sie_analytic(double r, double theta, double bsis, double a, double b);
double J1_sie_analytic(double r, double theta, double bsis, double a, double b);

double K0_sie_analytic(double r, double theta, double bsis, double a, double b);
double K1_sie_analytic(double r, double theta, double bsis, double a, double b);
double K2_sie_analytic(double r, double theta, double bsis, double a, double b);


double f1_sie_analytic(double theta, double pert_params[],double r_e){
  double conv_params[]={pert_params[2]};
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
  double X1 = r_e*cos(theta);
  double X2 = r_e*sin(theta);

  double J0 = J(0, conv_sis_circ, X1, X2, conv_params, a_1eps, a_2eps);
  double J1 = J(1, conv_sis_circ, X1, X2, conv_params, a_1eps, a_2eps);

  double out = r_e* ( a_1eps*a_2eps*( J0*pow(cos(theta),2) + J1*pow(sin(theta),2)) - 1);
  return out;
}


double Df0Dtheta_sie_analytic(double theta, double pert_params[],double r_e){
  double conv_params[]={pert_params[2]};
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
  double X1 = r_e*cos(theta);
  double X2 = r_e*sin(theta);

  double J0 = J(0, conv_sis_circ, X1, X2, conv_params, a_1eps, a_2eps);
  double J1 = J(1, conv_sis_circ, X1, X2, conv_params, a_1eps, a_2eps);

  double out = (a_1eps*a_2eps*r_e*r_e/2.0) * (J1 - J0)*sin(2.0*theta);

  return out;
}


double D2f0Dtheta2_sie_analytic(double theta, double pert_params[],double r_e){
  double conv_params[]={pert_params[2]};
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
  double X1 = r_e*cos(theta);
  double X2 = r_e*sin(theta);

  double J0 = J(0, conv_sis_circ, X1, X2, conv_params, a_1eps, a_2eps);
  double J1 = J(1, conv_sis_circ, X1, X2, conv_params, a_1eps, a_2eps);

  double K0 = K(0, conv_sis_circ_prime, X1, X2, conv_params, a_1eps, a_2eps);
  double K1 = K(1, conv_sis_circ_prime, X1, X2, conv_params, a_1eps, a_2eps);
  double K2 = K(2, conv_sis_circ_prime, X1, X2, conv_params, a_1eps, a_2eps);

  double out1 = a_1eps*a_2eps*r_e*r_e;
  double out2 = (J1 - J0)*cos(2.0*theta);
  double out3 = r_e*r_e*(K0 - 2.0*K1 + K2)*pow(sin(2.0*theta),2)/2.0;

  return out1*(out2 + out3);
}


double f1_sie(double theta, double pert_params[],double r_e){
  double conv_params[]={pert_params[2]};
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

  return  f_1_ellip_mass(conv_sis_circ,theta, conv_params, r_e, a_1eps, a_2eps);
}


double Df0Dtheta_sie(double theta, double pert_params[],double r_e){
  double conv_params[]={pert_params[2]};
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

  return  Df0Dtheta_ellip_mass(conv_sis_circ,theta, conv_params, r_e, a_1eps, a_2eps);
}



double D2f0Dtheta2_sie(double theta, double pert_params[],double r_e){
  double conv_params[]={pert_params[2]};
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

  return  D2f0Dtheta2_ellip_mass(conv_sis_circ, conv_sis_circ_prime, theta, conv_params, r_e, a_1eps, a_2eps);
}




double J0_sie_analytic(double r, double theta, double bsis, double a, double b){
  if(a*a > b*b){
    return (bsis*atan(sqrt(((fabs(-pow(a,2) + pow(b,2)))*pow(sin(theta),2))/(pow(b,2) + (pow(a,2) - pow(b,2))*pow(cos(theta),2)))))/(r*sqrt((fabs(-pow(a,2) + pow(b,2)))*pow(sin(theta),2)));
  }else{
    return (bsis*atanh(sqrt(((-pow(a,2) + pow(b,2))*pow(sin(theta),2))/(pow(b,2) + (pow(a,2) - pow(b,2))*pow(cos(theta),2)))))/(r*sqrt((-pow(a,2) + pow(b,2))*pow(sin(theta),2)));
  }
}

double J1_sie_analytic(double r, double theta, double bsis, double a, double b){
  if(a*a > b*b){
    return (bsis*atanh(sqrt(((pow(a,2) - pow(b,2))*pow(cos(theta),2))/(pow(b,2) + (pow(a,2) - pow(b,2))*pow(cos(theta),2)))))/(r*sqrt((pow(a,2) - pow(b,2))*pow(cos(theta),2)));
  }else{
    return (bsis*atan(sqrt(((fabs(pow(a,2) - pow(b,2)))*pow(cos(theta),2))/(pow(b,2) + (pow(a,2) - pow(b,2))*pow(cos(theta),2)))))/(r*sqrt((fabs(pow(a,2) - pow(b,2)))*pow(cos(theta),2)));
  }
  
}


double K0_sie_analytic(double r, double theta, double bsis, double a, double b){
  if(a*a > b*b){
    return (bsis*(pow(1.0/tan(theta),2)/
        sqrt(pow(b,2) + 
          (pow(a,2) - pow(b,2))*pow(cos(theta),2)) - 
       (atan(sqrt(((fabs(-pow(a,2) + pow(b,2)))*
               pow(sin(theta),2))/
             (pow(b,2) + 
               (pow(a,2) - pow(b,2))*pow(cos(theta),2))))
           *pow(1.0/sin(theta),2))/
        sqrt((fabs(-pow(a,2) + pow(b,2)))*pow(sin(theta),2))))/
   (2.*r);
  }else{
    return (bsis*(pow(1.0/tan(theta),2)/
        sqrt(pow(b,2) + 
          (pow(a,2) - pow(b,2))*pow(cos(theta),2)) - 
       (atanh(sqrt(((fabs(-pow(a,2) + pow(b,2)))*
               pow(sin(theta),2))/
             (pow(b,2) + 
               (pow(a,2) - pow(b,2))*pow(cos(theta),2))))
           *pow(1.0/sin(theta),2))/
        sqrt((fabs(-pow(a,2) + pow(b,2)))*pow(sin(theta),2))))/
   (2.*r);
  }
}



double K1_sie_analytic(double r, double theta, double bsis, double a, double b){
  return -bsis/(2.*r*sqrt(pow(b,2) + 
       (pow(a,2) - pow(b,2))*pow(cos(theta),2)));
}




double K2_sie_analytic(double r, double theta, double bsis, double a, double b){

  if(a*a > b*b){
    return (bsis*(-((atanh(sqrt(((fabs(pow(a,2) - pow(b,2)))*
                 pow(cos(theta),2))/
               (pow(b,2) + 
                 (pow(a,2) - pow(b,2))*pow(cos(theta),2))
               ))*pow(1.0/cos(theta),2))/
          sqrt(fabs((pow(a,2) - pow(b,2)))*pow(cos(theta),2)))\
        + pow(tan(theta),2)/
        sqrt(pow(b,2) + 
          (pow(a,2) - pow(b,2))*pow(cos(theta),2))))/
   (2.*r);
  }else{
    return (bsis*(-((atan(sqrt((fabs((pow(a,2) - pow(b,2)))*
                 pow(cos(theta),2))/
               (pow(b,2) + 
                 (pow(a,2) - pow(b,2))*pow(cos(theta),2))
               ))*pow(1.0/cos(theta),2))/
          sqrt((fabs(pow(a,2) - pow(b,2)))*pow(cos(theta),2)))\
        + pow(tan(theta),2)/
        sqrt(pow(b,2) + 
          (pow(a,2) - pow(b,2))*pow(cos(theta),2))))/
   (2.*r);
  }
}

#endif
