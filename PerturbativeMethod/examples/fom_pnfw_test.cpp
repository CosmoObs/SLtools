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
#include "../comparison_fns.h"
// to compile make: g++  -o exec_fom fom_pnfw_test.cpp `pkg-config gsl --cflags --libs`
int main(){
  

  FILE *out2;
  out2=fopen("gof_pnfw.dat","w");

  system(" echo  compiling the fortran code, please wait...  ");
  system(" g77-3.4 -Wall -o c_pnfw ../../lensing_models/pnfw/pnfw_contours.f ../../lensing_models/pnfw/root_finding.f ../../lensing_models/pnfw/lensing_functions.f ");
//   system(" gfortran-4.4 -o c_pnfw   fortran_srcs/*.f ");  
//   system(" gfortran-4.4 -Wall -o c_pnfw ../../lensing_models/pnfw/pnfw_contours.f ../../lensing_models/pnfw/root_finding.f ../../lensing_models/pnfw/lensing_functions.f ");

  system(" echo  ending the compilation  ");

   float argv[3];
   int iargv=1, iflagc=1;

  printf("Enter the convergence, scale radius,  ellipticity: ");
  scanf("%f %f %f",&argv[0],&argv[1], &argv[2]);
//     printf("Enter the convergence, scale radius: ");
//     scanf("%f %f",&argv[0],&argv[1]);
    double ks=argv[0];
    double rs=argv[1];
    double elp=argv[2];
//   printf("enter the flag control to choice the parameterization: ");
//   scanf("%i",&iargv);
    double lw_u=10.0;
    int npt=251;
//     double elp=0.0;
    double s_el=0.005;
  
//   for(int i=1;i<=199;i++){
      FILE *out1;
      out1=fopen("in_pnfw_par.txt","w");
  
//       elp=i*s_el;
      
      fprintf(out1,"%f %f %f %f %i %i\n",ks,rs,elp,lw_u,iflagc,npt);
      fclose(out1);    

      system("./c_pnfw");
      system("cp fort.61 cc.dat");
      system("cp fort.71 ca.dat");
  
      double pert_params[]={elp,1,ks,rs};
      double pot_params[] = {pert_params[2],pert_params[3]};

  
      double _r_e_nfw= r_e_nfw_find(pot_params);
      printf("Einstein Radius = %f\n",_r_e_nfw);

      double kappa_2=kappa2_nfw(pot_params,_r_e_nfw);
      double c3_model= c3_nfw(pot_params,_r_e_nfw);
      printf("the value of C3 is %f\n",c3_model);
//       plot_curves(f1_pnfw, Df0Dtheta_pnfw, D2f0Dtheta2_pnfw, pert_params, kappa_2,_r_e_nfw, npt, outtc, outcau, outsrc);
//       printf(" terminei as curvas");
      
      double mwsrfd_crit= d2_crit(f1_pnfw, D2f0Dtheta2_pnfw,kappa_2, pert_params, _r_e_nfw);
      printf(" The value of the D2, for kappa_s %f and ellipticity %f, is %f\n", pert_params[2],pert_params[0],mwsrfd_crit);
      double mwsrfd_caust= d2_caust(f1_pnfw, Df0Dtheta_pnfw, D2f0Dtheta2_pnfw,kappa_2, pert_params, _r_e_nfw);
      double dmax=lin_grad(f1_pnfw,D2f0Dtheta2_pnfw,c3_model, kappa_2,pert_params, _r_e_nfw,200);
      printf(" The value of the D2, for kappa_s %f and ellipticity %f, is %f\n", pert_params[2],pert_params[0],mwsrfd_caust);
      fprintf(out2," %f %f  %f %f \n",elp,mwsrfd_crit,mwsrfd_caust,dmax);  
//        }
      
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
