#include <cstdlib>
#include <cstdio>
#include <fstream>        //funcoes de entrada e saida para arquivos
#include <iostream>

#include "../perturbative_method.h"
#include "../pnfw_model2.h"
#include "../nfw_circular_model.h"
#include "../theta_find.h"
#include "../generate_arcs.h"
#include "../generate_curves.h"
#include "../arc_properties.h"
//#include "gof_per_method.h"

#define TAM_MAX1 20

// to compile make: g++ -Wall -o exec_name pnfw_test_new.cpp `pkg-config gsl --cflags --libs`

int main(){
  
//   double pert_params[4];
  FILE *entrada1;
  entrada1=fopen("input_lens.txt","r");
  float lens_par[TAM_MAX1]; 
  int i=0;
  
  if (entrada1== NULL) 
    {
      printf("Nao foi possivel abrir o arquivo.\n");
      lens_par[1]=0.8,lens_par[2]=1.5,lens_par[3]=0. ;
//       break;
//        exit();
     }else{
      while (! feof(entrada1))
      {
        i++;
        fscanf(entrada1," %f",&lens_par[i]);
	
      }}
//       int nlp = i++;
  
  int npts = 4001;
//   double twpi= 6.283185308;  
  double pert_params[]={lens_par[3],1,lens_par[1],lens_par[2]};
  double pot_params[] = {pert_params[2],pert_params[3]};
  
  double _r_e_nfw= r_e_nfw_find(pot_params);
  printf("Einstein Radius = %f\n",_r_e_nfw);

  double kappa_2=kappa2_nfw(pot_params,_r_e_nfw);

  elliptical_source source;

  source.x0 = (_r_e_nfw/15.5);
  source.y0 = source.x0;
//   source.y0 = 0.0;
  printf(" Source Position %f %f\n",source.x0, source.y0);
  source.R0 = (sqrt(2.0)/60.0)*_r_e_nfw;
  source.eta0 = 0.0;
//   source.theta0 = twpi/8.0;
  source.theta0 = 0.0;

  FILE *outarc = fopen ("arcs_pnfw.dat" , "w");
  FILE *outls = fopen ("lensing_data.txt" , "w");  
  FILE *outtc = fopen ("tang_crit.dat" , "w");
  FILE *outcau = fopen ("tang_caust.dat" , "w");
  FILE *outsrc = fopen ("src_plot.dat" , "w");
  
  plot_arcs(source, f1_pnfw, Df0Dtheta_pnfw, pert_params, kappa_2,_r_e_nfw,npts, outarc);
  plot_curves(source, f1_pnfw, Df0Dtheta_pnfw, D2f0Dtheta2_pnfw, pert_params, kappa_2,_r_e_nfw, npts, outtc, outcau, outsrc);
  

  fprintf(outls,"\t\t PNFW Model Parameters: \n");  
  fprintf(outls," Scale Radius = %f\n",pot_params[1]);
  fprintf(outls," Characteristic Convergence = %f\n",pot_params[0]);
  fprintf(outls," Ellipticity of the Lensing Potential =  %f\n\n",pert_params[0]);
  fprintf(outls,"\t\t Source Data: \n");  
  fprintf(outls," Source Position %f %f\n",source.x0,source.y0);
  fprintf(outls," Source Radius %f\n",source.R0);
  fprintf(outls," Source Ellipticity and Inclination %f %f\n\n",source.eta0,source.theta0);

  arc_measures(f1_pnfw, Df0Dtheta_pnfw, kappa_2, source, pert_params, _r_e_nfw, npts,outls);  

  FILE *outls2 = fopen ("input_src.txt","w"); 
  int step=1, nang=400;

  fprintf(outls2, "%f\n", source.x0);
  fprintf(outls2, "%f\n", source.y0);
  fprintf(outls2, "%f\n", source.eta0);
  fprintf(outls2, "%f\n", source.theta0);  
  fprintf(outls2, "%f\n", source.R0);
  fprintf(outls2, "%f\n", source.R0);
  fprintf(outls2, "%i\n", step);
  fprintf(outls2, "%i\n", nang);

 // double gof_crit= gof_per_method(f1_pnfw, Df0Dtheta_pnfw,D2f0Dtheta2_pnfw,kappa_2, pert_params, _r_e_nfw, 1);
 // printf(" The value of the gof, for kappa_s %f and ellipticity %f, is %f\n", pert_params[2],pert_params[0],gof_crit);
  
 // double gof_caust= gof_per_method(f1_pnfw, Df0Dtheta_pnfw,D2f0Dtheta2_pnfw,kappa_2, pert_params, _r_e_nfw, 2);
 // printf(" The value of the gof, for kappa_s %f and ellipticity %f, is %f\n", pert_params[2],pert_params[0],gof_caust);  

}
