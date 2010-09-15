#include <cstdio>
#include <cstdlib>

#include "perturbative_method.h"

void theta(f_type Df0Dtheta_in, elliptical_source source_in,int N=100){

  double theta[N];
  double arg[N];
  int N_pair=0;

  for(int i=0;i<N;i++) theta[i] = double(i)*2.0*3.14159/double(N);
  for(int i=0;i<N;i++) arg[i] = arg_sqrt(Df0Dtheta_in, source_in, theta[i], NULL);


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
  test=arg_sqrt(Df0Dtheta_in, source_in, (theta_N[0]+theta_N[1])/2.0, NULL);




  
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

}
