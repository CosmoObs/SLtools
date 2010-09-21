#include <cstdio>
#include <cmath> 
#include <cstdlib>
#include <iostream>

#include "../perturbative_method.h"

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


void print_arcs_SIS(elliptical_source source_in, int npts=20){
  
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


void print_arcs_TEST(elliptical_source source_in, int npts=20){
  
  double theta = 0.0;
  double r_p = 0.0;
  double r_m = 0.0;
  double pot_params[3];
  pot_params[0] = 0.0;
  pot_params[1] = 2.;
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

int main(int argc, char *argv[]){

  int npts;

  elliptical_source source;

  if (argc < 2) {
    // printf("No input args, using defaults; To input, use this format: \n ");
    // printf("Num pts   Source parameters: x0   y0  R0 (size)  eta (ellip)   theta (ang. position)\n");

    npts = 100;
    
    source.x0 = 1;
    source.y0 = 0.0;
    source.R0 = 0.9;
    source.eta0 = 0.0;
    source.theta0 =  0.0;
    
  }
  else{
    npts = atoi(argv[1]);

    source.x0 = atof(argv[2]);
    source.y0 = atof(argv[3]);
    source.R0 = atof(argv[4]);
    source.eta0 = atof(argv[5]);
    source.theta0 =  atof(argv[6]);

  }

  /*   printf(" Using: npts %i  \n", npts);
   printf(" x0,y0 = (%f,%f) \n", source.x0,source.y0);
   printf(" R0, eta0, theta0 = (%f,%f,%f) \n",source.R0, source.eta0, source.theta0 );
  */

  print_arcs_TEST(source, npts);



  return 0;
}
