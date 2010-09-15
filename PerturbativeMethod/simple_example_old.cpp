#include <cstdio>
#include <cmath> 

#include "perturbative_method.h"

//to compile: g++ -Wall simple_example.cpp


//Quantities related to the potential
/*********************************************************************************************************************/
//!  Function defined in http://arxiv.org/abs/0804.4277, for SIS lens model
//!
//!  \f$ f_{n}(\theta) \equiv \frac{1}{n!} \left[\frac{\partial\psi(r,\theta)}{\partial r}\right]_{r=1} \f$
/*!
  \param theta angular coordinate in the lens plane
  \param pot_params[] potential parameters
  \return 0.0
*/
double f1_SIS(double theta, double pot_params[]){
  return 0.0;
}

double Df0Dtheta_SIS(double theta, double pot_params[]){
  return 0.0;
}

/*********************************************************************************************************************/


//pot_params[0] = mp
//pot_params[1] = rp
//pot_params[2] = thetap
double f1_TEST(double theta, double pot_params[]){
  double tmp1 = pot_params[0] * (1.0 - pot_params[1]*cos(theta-pot_params[2]) );
  double tmp2 = sqrt ( 1.0 - 2.0*pot_params[1]*cos(theta-pot_params[2]) + pot_params[1]*pot_params[1]);

  return tmp1/tmp2 ;
}

double Df0Dtheta_TEST(double theta, double pot_params[]){
  double tmp1 = pot_params[0] * (pot_params[1]*sin(theta-pot_params[2]) );
  double tmp2 = sqrt ( 1.0 - 2.0*pot_params[1]*cos(theta-pot_params[2]) + pot_params[1]*pot_params[1]);
  return tmp1/tmp2;
}



/*********************************************************************************************************************/


void print_arcs_SIS(elliptical_source source_in, int npts=200){
  
  double theta = 0.0;
  double r_p = 0.0;
  double r_m = 0.0;

  for(int i=0;i<=npts;i++){

    if(arg_sqrt(Df0Dtheta_SIS, source_in, theta, NULL)>0.0){
      r_p = 1.0+dr_plus(f1_SIS, Df0Dtheta_SIS, 1.0, source_in, theta, NULL);
      r_m = 1.0+dr_minus(f1_SIS, Df0Dtheta_SIS, 1.0, source_in, theta, NULL);
      printf("%E %E\n",r_p*cos(theta),r_p*sin(theta));
      printf("%E %E\n",r_m*cos(theta),r_m*sin(theta));
    }
    theta+= 6.283185308/npts;
  }

}


void print_arcs_TEST(elliptical_source source_in, int npts=200){
  
  double theta = 0.0;
  double r_p = 0.0;
  double r_m = 0.0;
  double pot_params[3];
  pot_params[0] = 0.01;
  pot_params[1] = 1.3;
  pot_params[2] = 0.0;


  for(int i=0;i<=npts;i++){

    if(arg_sqrt(Df0Dtheta_TEST, source_in, theta, pot_params)>0.0){
      r_p = 1.0+dr_plus(f1_TEST, Df0Dtheta_TEST, 1.0, source_in, theta, pot_params);
      r_m = 1.0+dr_minus(f1_TEST, Df0Dtheta_TEST, 1.0, source_in, theta, pot_params);
      printf("%E %E\n",r_p*cos(theta),r_p*sin(theta));
      printf("%E %E\n",r_m*cos(theta),r_m*sin(theta));
    }
    theta+= 6.283185308/npts;
  }

}

int main(){
  elliptical_source source;
  source.x0 = -0.002;
  source.y0 = -0.0375;
  source.R0 = 0.0252;
  source.eta0 = 0.0;
  source.theta0 =  2.48814;


  print_arcs_TEST(source);

  return 0;
}
