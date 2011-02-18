#include <cstdlib>
#include <cstdio>
#include <stdio.h>
#include <string.h>
#include <fstream>        //funcoes de entrada e saida para arquivos
#include <iostream>

#include "../perturbative_method.h"
#include "../pnfw_model2.h"
#include "../nfw_circular_model.h"
#include "../theta_find.h"
#include "../generate_arcs.h"
#include "../generate_curves.h"
#include "../arc_properties.h"
#include "../arc_cross_section.h"
//#include "gof_per_method.h"

#define TAM_MAX1 20

// to compile make: g++ -Wall -o exec_pnfw pnfw_test_new.cpp `pkg-config gsl --cflags --libs`

int main(){
  
//   double pert_params[4];
  FILE *entrada1;
  entrada1=fopen("input_lens.txt","r");
  float argv[3];
  float lens_par[TAM_MAX1];
  int iargv=1;
  
  

  int i=0;
  if (entrada1== NULL) 
    {
      printf("Nao foi possivel abrir o arquivo... then \n");
      printf("Enter the convergence, scale radius,  ellipticity: ");
      scanf("%f %f %f",&argv[0],&argv[1], &argv[2]);
      lens_par[1]=argv[0];
      lens_par[2]=argv[1];
      lens_par[3]=argv[2];

      printf("enter the flag control to choice the parameterization: ");
      scanf("%i",&iargv);
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
  double twpi= 6.283185308;  
  double pert_params[]={lens_par[3],iargv,lens_par[1],lens_par[2]};
  double pot_params[] = {pert_params[2],pert_params[3]};

  double _r_e_nfw= r_e_nfw_find(pot_params);   //func
  printf("Einstein Radius = %f\n",_r_e_nfw);

  double kappa_2=kappa2_nfw(pot_params,_r_e_nfw);   //func

  elliptical_source source;

  source.x0 = (_r_e_nfw/2);
//   source.y0 = source.x0;
  source.y0 = 0.0;
  printf(" Source Position %f %f\n",source.x0, source.y0);
  source.R0 = (sqrt(2.0)/60.0)*_r_e_nfw;
  source.eta0 = 0.0;
   source.theta0 = twpi/8.0;
  //source.theta0 = 0.0;

  FILE *outarc = fopen ("arcs_pnfw.dat" , "w");
  FILE *outls = fopen ("lensing_data.txt" , "w");  
  FILE *outtc = fopen ("tang_curves.dat" , "w");
  FILE *outsrc = fopen ("src_plot.dat" , "w");
  FILE *outrthp = fopen ("curves_r_th_pos.dat" , "w");
  FILE *outrthn = fopen ("curves_r_th_neg.dat" , "w");
  
  plot_arcs(source, f1_pnfw, Df0Dtheta_pnfw, pert_params, kappa_2,_r_e_nfw,npts, outarc);   //func
  plot_curves(f1_pnfw, Df0Dtheta_pnfw, D2f0Dtheta2_pnfw, pert_params, kappa_2,_r_e_nfw, npts, outtc);   //func
  plot_rth_curves(f1_pnfw, Df0Dtheta_pnfw, D2f0Dtheta2_pnfw, pert_params, kappa_2, _r_e_nfw, 10., npts, outrthp, outrthn);   //func
  plot_sources(source,npts,outsrc);   //func

  fprintf(outls,"\t\t PNFW Model Parameters: \n");  
  fprintf(outls," Scale Radius = %f\n",pot_params[1]);
  fprintf(outls," Characteristic Convergence = %f\n",pot_params[0]);
  fprintf(outls," Ellipticity of the Lensing Potential =  %f\n\n",pert_params[0]);
  fprintf(outls,"\t\t Source Data: \n");  
  fprintf(outls," Source Position %f %f\n",source.x0,source.y0);
  fprintf(outls," Source Radius %f\n",source.R0);
  fprintf(outls," Source Ellipticity and Inclination %f %f\n\n",source.eta0,source.theta0);

  arc_measures(f1_pnfw, Df0Dtheta_pnfw, kappa_2, source, pert_params, _r_e_nfw, npts,outls);   //func

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

   fclose(outls2);
   fclose(outarc);
   fclose(outls); 
   fclose(outtc); 
   fclose(outsrc); 
   fclose(outrthp); 
   fclose(outrthn);
   
  system("xmgrace -view 0.15 0.15 0.85 0.85 -block curves_r_th_pos.dat -bxy 3:4  -block tang_curves.dat -bxy 3:4 -block curves_r_th_neg.dat -bxy 3:4 src_plot.dat");
  system("xmgrace -view 0.15 0.15 0.85 0.85 curves_r_th_pos.dat tang_curves.dat curves_r_th_neg.dat  arcs_pnfw.dat");
//   system("xmgrace -view 0.15 0.15 0.85 0.85 tang_crit.dat arcs_pnfw.dat");
 
  double sig_test = sigma_defarcs(f1_pnfw, D2f0Dtheta2_pnfw, pert_params, kappa_2,_r_e_nfw, 10);
    printf("cross-sec=%E\n",sig_test);
}
