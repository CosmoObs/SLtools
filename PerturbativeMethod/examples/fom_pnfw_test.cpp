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
// #include "../arc_properties.h"
#include "../gof_per_method.h"

// to compile make: g++  -o exec_fom fom_pnfw_test.cpp `pkg-config gsl --cflags --libs`
int main(){
  

  FILE *out2;
  out2=fopen("gof_pnfw.dat","w");

  system(" echo  compiling the fortran code, please wait...  ");
//   system(" g77-3.4 -Wall -o c_pnfw   fortran_srcs/*.f ");
  system(" gfortran-4.4 -o c_pnfw   fortran_srcs/*.f ");  
  system(" echo  ending the compilation  ");

  double ks=1.0, rs=2.0;
  int npt=251;
//   double elp=0.1;
  double elp=0.0;
  double s_el=0.005;

  for(int i=1;i<=199;i++){
      FILE *out1;
      out1=fopen("in_pnfw_par.txt","w");
  
      elp=i*s_el;
      
      fprintf(out1,"%f %f %f %i\n",ks,rs,elp,npt);
      fclose(out1);    
      system("./c_pnfw");
      system("cp fort.61 cc.dat");
      system("cp fort.71 ca.dat");
//       system("rm in_pnfw_par.txt");
  
      double pert_params[]={elp,1,ks,rs};
      double pot_params[] = {pert_params[2],pert_params[3]};

//       printf("counter %i, k_s %f, r_s=%f, ellipticity %f\n", i,pert_params[2],pert_params[3],pert_params[0]);
  
      double _r_e_nfw= r_e_nfw_find(pot_params);
      printf("Einstein Radius = %f\n",_r_e_nfw);

      double kappa_2=kappa2_nfw(pot_params,_r_e_nfw);

      elliptical_source source;

      source.x0 = (_r_e_nfw/15.5);
      source.y0 = source.x0;
      printf(" Source Position %f %f\n",source.x0, source.y0);
      source.R0 = (sqrt(2.0)/60.0)*_r_e_nfw;
      source.eta0 = 0.0;
      source.theta0 = 0.0;

      FILE *outtc = fopen ("tang_crit.dat" , "w");
      FILE *outcau = fopen ("tang_caust.dat" , "w");
      FILE *outsrc = fopen ("src_plot.dat" , "w");
      
      

      plot_curves(f1_pnfw, Df0Dtheta_pnfw, D2f0Dtheta2_pnfw, pert_params, kappa_2,_r_e_nfw, npt, outtc, outcau, outsrc);
      printf(" terminei as curvas");
      
      double gof_crit= gof_per_method(f1_pnfw, Df0Dtheta_pnfw,D2f0Dtheta2_pnfw,kappa_2, pert_params, _r_e_nfw, 1);
      printf(" The value of the gof, for kappa_s %f and ellipticity %f, is %f\n", pert_params[2],pert_params[0],gof_crit);
  
      double gof_caust= gof_per_method(f1_pnfw, Df0Dtheta_pnfw,D2f0Dtheta2_pnfw,kappa_2, pert_params, _r_e_nfw, 2);
      printf(" The value of the gof, for kappa_s %f and ellipticity %f, is %f\n", pert_params[2],pert_params[0],gof_caust);
      fprintf(out2," %f %f %f \n",elp,gof_crit,gof_caust);
       }
      
      fclose(out2);

  return 0.0;
}



/*
//   
// //   plot_arcs(source, f1_pnfw, Df0Dtheta_pnfw, pert_params, kappa_2,_r_e_nfw,npts, outarc);
//    
// 
//     system("xmgrace tang_crit.dat tang_caust.dat fort.65 fort.75 ");
*/
