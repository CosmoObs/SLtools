#include <cstdlib>
#include <cstdio>
#include <fstream>        // Functions for input and output files
//#include <cstdio>
#include <iostream>

#include "../perturbative_method.h"
#include "../psis_model.h"
#include "../sis_model.h"
#include "../theta_find.h"
#include "../generate_arcs.h"
#include "../generate_curves.h"
#include "../arc_properties.h" 
#define TAM_MAX 20


// To compile: g++ -o exec_name pe_sis_test.cpp `pkg-config gsl --cflags --libs`
// To running: exec_name


int main(){
  
  FILE *outarc = fopen ("arcs_pesis.dat" , "w");
  FILE *outls = fopen ("lensing_data_pesis.txt" , "w");  
  FILE *outtc = fopen ("tang_crit_pesis.dat" , "w");
  FILE *outcau = fopen ("tang_caust_pesis.dat" , "w");
  FILE *outsrc = fopen ("src_plot_pesis.dat" , "w");

   int npts = 4000;

  double twpi= 6.283185308;  
//   double pert_params[]={0.64,1.5, 0.1};
//   double pot_params[]={pert_params[0],pert_params[1],pert_params[2]};
  double pert_params[]={0.1,3,0.64};
  double pot_params[] = {pert_params[2]};
  
  double _r_e_sis= pot_params[0];
  printf("Einstein Radius = %f\n",_r_e_sis);
  double kappa_2 = kappa2_sis(pot_params,_r_e_sis);

  elliptical_source source;

  source.x0 = (_r_e_sis/10.);
  source.y0 = source.x0;
  printf(" Source Position %f %f\n",source.x0, source.y0);
  source.R0 = (sqrt(2.0)/25.0)*_r_e_sis;
  source.eta0 = 0.0;
  source.theta0 = twpi/8.0;

  fprintf(outls," PE_SIS Model Parameters: \n");  
  fprintf(outls," Mass = %f\n",pot_params[0]);
  fprintf(outls," Ellipticity of the Lensing Potential =  %f\n\n",pert_params[0]);
  fprintf(outls," Source Data: \n");  
  fprintf(outls," Source Position %f %f\n",source.x0,source.y0);
  fprintf(outls," Source Radius %f\n",source.R0);
  fprintf(outls," Source Ellipticity and Inclination %f %f\n\n",source.eta0,source.theta0);

  plot_arcs(source, f1_psis, Df0Dtheta_psis, pert_params, kappa_2,_r_e_sis,npts, outarc);
  plot_curves(f1_psis, Df0Dtheta_psis, D2f0Dtheta2_psis, pert_params, kappa_2,_r_e_sis, npts, outtc, outcau);
  arc_measures(f1_psis, Df0Dtheta_psis, kappa_2, source, pert_params, _r_e_sis, npts,outls);   
}
