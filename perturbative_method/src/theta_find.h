/** @file
* Example of doxygen documentation for C functions FIXME. 
*/

/** @package theta_find
*  Short description FIXME 
*   
*
*  Detailed description FIXME
*
*/

#ifndef THETA_FIND_H
#define THETA_FIND_H

#include <cstdio>
#include <cstdlib>

#include "perturbative_method.h"

void theta_find(f_type Df0Dtheta_in, elliptical_source source_in,  double pert_params[], double _r_e, int N=1000){

  double theta[N];
  double arg[N];
  int N_pair=0;

  for(int i=0;i<N;i++) theta[i] = double(i)*2.0*3.141592654/double(N);
  for(int i=0;i<N;i++) arg[i] = arg_sqrt(Df0Dtheta_in, source_in, theta[i], pert_params,_r_e);


  for(int i=0;i<N;i++){
    if(arg[i]*arg[i+1]  < 0){
      //printf("%f  %f  %f\n",theta[i], arg[i],arg[i+1]);
      N_pair++;
    }
  }
  double theta_N[N_pair+1];
  int count = 0;
  for(int i=0;i<N;i++){
    if(arg[i]*arg[i+1]  < 0){
      theta_N[count] = theta[i];
      count++;
    }
  }

  theta_N[N_pair] =  theta_N[0];
  
  double test = 0;
  test=arg_sqrt(Df0Dtheta_in, source_in, (theta_N[0]+theta_N[1])/2.0, pert_params,_r_e);


  
  FILE * pFile = fopen ("theta_out.txt","w");
  fprintf(pFile,"%i\n",N_pair/2);

  if(test>0){
    for(int i=0;i<N_pair;i+=2){
      fprintf(pFile,"%f %f\n",theta_N[i],theta_N[i+1]);
    }
  } else {
    for(int i=0;i<N_pair;i+=2){
      fprintf(pFile,"%f %f\n",theta_N[i+1],theta_N[i+2]);
    }
  }

  fclose(pFile);

}

void theta_find(f_type Df0Dtheta_in, elliptical_source source_in,  double pert_params[], double _r_e, double out[10][2], int N=1000){

  double theta[N];
  double arg[N];
  int N_pair=0;

  for(int i=0;i<N;i++) theta[i] = double(i)*2.0*3.141592654/double(N);
  for(int i=0;i<N;i++) arg[i] = arg_sqrt(Df0Dtheta_in, source_in, theta[i], pert_params,_r_e);


  for(int i=0;i<N;i++){
    if(arg[i]*arg[i+1]  < 0){
      //printf("%f  %f  %f\n",theta[i], arg[i],arg[i+1]);
      N_pair++;
    }
  }
  double theta_N[N_pair+1];
  int count = 0;
  for(int i=0;i<N;i++){
    if(arg[i]*arg[i+1]  < 0){
      theta_N[count] = theta[i];
      count++;
    }
  }

  theta_N[N_pair] =  theta_N[0];
  
  double test = 0;
  test=arg_sqrt(Df0Dtheta_in, source_in, (theta_N[0]+theta_N[1])/2.0, pert_params,_r_e);


  
//  FILE * pFile = fopen ("theta_out.txt","w");
//  fprintf(pFile,"%i\n",N_pair/2);
 double d_theta = 2.0*3.141592654/double(N);
 out[0][0] = int(N_pair/2);
 out[0][1] = d_theta;

  int n_pos = 1;
  if(test>0){
    for(int i=0;i<N_pair;i+=2){
      //fprintf(pFile,"%f %f\n",theta_N[i],theta_N[i+1]);
      double test_inf = arg_sqrt(Df0Dtheta_in, source_in, theta_N[i], pert_params,_r_e);
      while(test_inf<0){
	theta_N[i]+=d_theta;
	test_inf = arg_sqrt(Df0Dtheta_in, source_in, theta_N[i], pert_params,_r_e);
      }
      out[n_pos][0] = theta_N[i];
      double test_sup = arg_sqrt(Df0Dtheta_in, source_in, theta_N[i+1], pert_params,_r_e);
      while(test_sup<0){
	theta_N[i+1]-=d_theta;
	test_sup = arg_sqrt(Df0Dtheta_in, source_in, theta_N[i+1], pert_params,_r_e);
      }      
      out[n_pos][1] = theta_N[i+1];
      n_pos++;
    }
  } else {
    for(int i=0;i<N_pair;i+=2){
      //fprintf(pFile,"%f %f\n",theta_N[i+1],theta_N[i+2]);
      double test_inf = arg_sqrt(Df0Dtheta_in, source_in, theta_N[i+1], pert_params,_r_e);
      while(test_inf<0){
	theta_N[i+1]+=d_theta;
	test_inf = arg_sqrt(Df0Dtheta_in, source_in, theta_N[i+1], pert_params,_r_e);
      }
      out[n_pos][0] = theta_N[i+1];
      double test_sup = arg_sqrt(Df0Dtheta_in, source_in, theta_N[i+2], pert_params,_r_e);
      while(test_sup<0){
	theta_N[i+2]-=d_theta;
	test_sup = arg_sqrt(Df0Dtheta_in, source_in, theta_N[i+2], pert_params,_r_e);
      }
      out[n_pos][1] = theta_N[i+2];
      n_pos++;
    }
  }

//  Corrected values

}
#endif
