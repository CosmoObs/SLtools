#include "../theta_find.h"
#include  "../perturbative_method.h" 

#include <cstdio>
#include <cstdlib>


double f1_SIS(double theta, double pot_params[]){
  return 0.0;
}
double Df0Dtheta_SIS(double theta, double pot_params[]){
  return 0.0;
}

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


int main(){

  elliptical_source source;
  source.x0     = -0.1; 
  source.y0     = 0.0;
  source.R0     = 0.0308;
  source.eta0   = 0.0;
  source.theta0 = 0.0;

  theta_find(Df0Dtheta_SIS,source,NULL,2000);

  print_arcs_SIS(source,2000);

  return 0;
}
