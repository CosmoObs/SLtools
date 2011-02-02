#ifndef INTEGRALS_H
#define INTEGRALS_H

#include <cstddef>

#include "general_methods.h"

typedef double (*conv_type) (double r, double conv_params[]);

/*---introducing elliptical symmetry---*****************************************/
double CSI2(double X1, double X2, double a, double b){
  return sqrt( pow(X1/a,2) + pow(X2/b,2) );
}


double CSI2U(double X1, double X2, double a, double b, double u){
  return sqrt(u*( X1*X1/( 1.0 - (1.0-a*a)*u ) + X2*X2/( 1.0 - (1.0-b*b)*u ) ));
}


double kCSI2(conv_type conv_in,double X1, double X2, double conv_params[], double a, double b){
  double radial_coord = CSI2(X1,X2,a,b);
  return conv_in(radial_coord,conv_params);
}


double kCSI2U(conv_type conv_in,double X1, double X2, double conv_params[], double u, double a, double b){
  double radial_coord = CSI2U(X1,X2,a,b,u);
  return conv_in(radial_coord,conv_params);
}
/******************************************************************************/



/*---defining the arguments of the integrals---**********************************/
double AJ(int n, conv_type conv_in,double X1, double X2, double conv_params[], double u, double a, double b){
  double tmp_a = pow( 1.0 - (1.0-a*a)*u , 1.5 - n);
  double tmp_b = pow( 1.0 - (1.0-b*b)*u , 0.5 + n);
  return kCSI2U(conv_in,X1,X2,conv_params,u,a,b)/(tmp_a*tmp_b);
}

double AK(int n, conv_type conv_prime_in,double X1, double X2, double conv_params[], double u, double a, double b){
  double tmp_a = pow( 1.0 - (1.0-a*a)*u , 2.5 - n);
  double tmp_b = pow( 1.0 - (1.0-b*b)*u , 0.5 + n);
  return u*kCSI2U(conv_prime_in,X1,X2,conv_params,u,a,b)/(tmp_a*tmp_b);
}

//!  J integral argument in the form double(double u, double par[])
//!
/*!
  \param u mute variable
  \param par[0] n, from J_n
  \param par[1] X1
  \param par[2] X2
  \param par[3] a
  \param par[4] b
  \param par[5 - ] lens parameters
  \return 
*/
double AJ_final(double u, double par[], conv_type conv_in){
  int n = int(par[0]);
  double X1 = par[1];
  double X2 = par[2];
  double  a = par[3];
  double  b = par[4];

  //printf("%f %f %f %f %f ",par[0],par[1],par[2],par[3] ,par[4]);

  double conv_params[10];
  //int i=5;
  for(int i=5;i<15;i++){
    conv_params[i-5] = par[i];
  }
  //printf("%f %f %f\n",par[5],par[6],par[7]);

  return AJ(n, conv_in, X1,  X2,  conv_params, u, a, b);
}

//!  K integral argument in the form double(double u, double par[])
//!
/*!
  \param u mute variable
  \param par[0] n, from J_n
  \param par[1] X1
  \param par[2] X2
  \param par[3] a
  \param par[4] b
  \param par[5 - ] lens parameters
  \return 
*/
double AK_final(double u, double par[], conv_type conv_in){
  int n = int(par[0]);
  double X1 = par[1];
  double X2 = par[2];
  double  a = par[3];
  double  b = par[4];

  //printf("%f %f %f %f %f ",par[0],par[1],par[2],par[3] ,par[4]);

  double conv_params[10];
  //int i=5;
  for(int i=5;i<15;i++){
    conv_params[i-5] = par[i];
  }
  //printf("%f %f %f\n",par[5],par[6],par[7]);

  return AK(n, conv_in, X1,  X2,  conv_params, u, a, b);
}
/******************************************************************************/


double J(int n, conv_type conv_in,double X1, double X2, double conv_params[], double a, double b){

  double int_params[15];
  int_params[0] = n;
  int_params[1] = X1;
  int_params[2] = X2;
  int_params[3] = a;
  int_params[4] = b;

  //printf("%f %f %f %f %f ",int_params[0],int_params[1],int_params[2],int_params[3] ,int_params[4]);

  for(int i=5;i<15;i++){
    int_params[i] = conv_params[i-5];
  }
  //printf("%f %f %f\n",int_params[5],int_params[6],int_params[7]);

  return IntGaulegSub50_elip_model(AJ_final, int_params, 0.0, 1.0, conv_in);
}


double K(int n, conv_type conv_in,double X1, double X2, double conv_params[], double a, double b){

  double int_params[15];
  int_params[0] = n;
  int_params[1] = X1;
  int_params[2] = X2;
  int_params[3] = a;
  int_params[4] = b;

  //printf("%f %f %f %f %f ",int_params[0],int_params[1],int_params[2],int_params[3] ,int_params[4]);

  for(int i=5;i<15;i++){
    int_params[i] = conv_params[i-5];
  }
  //printf("%f %f %f\n",int_params[5],int_params[6],int_params[7]);

  return IntGaulegSub50_elip_model(AK_final, int_params, 0.0, 1.0, conv_in);
}


#endif
