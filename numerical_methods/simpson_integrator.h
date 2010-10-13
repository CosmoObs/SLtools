/** @file
* Example of doxygen documentation for C functions FIXME. 
*/

/** @package general_methods
*  Package with general functions, root finders, integrators etc.
*
*  Detailed descrition FIXME
*
*/
#include <cstdlib>
#include <cstdio>
#include <fstream>        //funcoes de entrada e saida para arquivos
#include <iostream>

#ifndef SYMPSON_INT_H
#define SYMPSON_INT_H

#include <cmath> 

// typedef double (*f_type)  (double r); 

typedef double (*f_type)  (double theta, double pot_params[], double src_params[], double _r_e); 

// double simpson_int(f_type fx_in, double x_lim[2], int npts=2001 ){
double simpson_int(f_type fx_in, double pot_params[], double src_params[], double _r_e, double x_lim[2], int npts=2001){
  double x[npts], ftemp[npts], result;
  double dx=(x_lim[1]-x_lim[0])/(npts-1);
  for (int i=1; i<=npts;i++){
      x[i]=x_lim[0]+(i-1)*dx;
//     ftemp[i]=fx_in(x[i]);
    ftemp[i]=fx_in(x[i], pot_params, src_params, _r_e);
    }
  int mpts=(npts-1)/2;

  double sum_1=0.0, sum_2=0.0, sum_3=0.0;

  for (int k=1; k<=mpts; k++){
      int k1=2*k-1, k2=2*k,k3=2*k+1;
      sum_1+=ftemp[k1];
      sum_2+=ftemp[k2];
      sum_3+=ftemp[k3];
    }
  printf("the sums are %f %f and %f \n", sum_1, sum_2, sum_3);
  result=(dx/3.0)*(sum_1+4.0*sum_2+sum_3);
  printf("The result of the integral is %f\n",result);
  return result;  
}

#endif