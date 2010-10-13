#include <cstdlib>
#include <cstdio>
#include <fstream>        //funcoes de entrada e saida para arquivos
//#include <cstdio>
#include <iostream>

#include "../perturbative_method.h"
#include "../pnfw_model.h"
#include "../theta_find.h"
#include "../generate_arcs.h"
#include "../generate_curves.h"
#include "../arc_properties.h"

#define TAM_MAX 20

// to compile make: g++ -Wall pnfw_test.cpp `pkg-config gsl --cflags --libs`

int main(){

  FILE *outarc = fopen ("arcs_pnfw.dat" , "w");
  FILE *outls = fopen ("lensing_data.dat" , "w");  
  FILE *outtc = fopen ("tang_crit.dat" , "w");
  FILE *outcau = fopen ("tang_caust.dat" , "w");
  FILE *outsrc = fopen ("src_plot.dat" , "w");
 

  int npts = 4000;

  double twpi= 6.283185308;  
  double pot_params[] = {1.5,0.65,0.25};
//   double theta_in[TAM_MAX];
//   double theta_sup[TAM_MAX];

  double _r_e_nfw= r_e_nfw_find(pot_params);
  printf("Einstein Radius = %f\n",_r_e_nfw);

  elliptical_source source;
/*
  source.x0 = 0.04;
  source.y0 = 0.04;
  source.R0 = 0.01;
  source.eta0 = 0.0;
  source.theta0 = 0.0;
*/
  source.x0 = (_r_e_nfw/4.);
  source.y0 = source.x0;
  printf(" Source Position %f %f\n",source.x0,source.y0);
  source.R0 = (sqrt(2.0)/25.0)*_r_e_nfw;
  source.eta0 = 0.8;
  source.theta0 = twpi/10.0;
//   source.theta0 = 0.0;

  fprintf(outls," PNFW Model Parameters: \n");  
  fprintf(outls," Scale Radius = %f\n",pot_params[0]);
  fprintf(outls," Characteristic Convergence = %f\n",pot_params[1]);
  fprintf(outls," Ellipticity of the Lensing Potential =  %f\n\n",pot_params[2]);
  fprintf(outls," Source Data: \n");  
  fprintf(outls," Source Position %f %f\n",source.x0,source.y0);
  fprintf(outls," Source Radius %f\n",source.R0);
  fprintf(outls," Source Ellipticity and Inclination %f %f\n\n",source.eta0,source.theta0);

  
//   theta_find(Df0Dtheta_pnfw, source,pot_params,_r_e_nfw, theta_in, theta_sup);
  theta_find(Df0Dtheta_pnfw, source,pot_params,_r_e_nfw);
//   printf("algun angle\%f\t\%f\n", theta_in[0],theta_sup[0])
  

// 
  plot_arcs(source, f1_pnfw, Df0Dtheta_pnfw, pot_params, kappa2_nfw(pot_params,_r_e_nfw),_r_e_nfw,npts, outarc);
  plot_curves(source, f1_pnfw, Df0Dtheta_pnfw, D2f0Dtheta2_pnfw, pot_params, kappa2_nfw(pot_params,_r_e_nfw),_r_e_nfw, npts, outtc, outcau, outsrc);
}


}
