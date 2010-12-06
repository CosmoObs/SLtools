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
// to compile make: g++  -o exec_gen fom_ks-e_space_pnfw.cpp `pkg-config gsl --cflags --libs`
int main(){
  

  FILE *out2;
  out2=fopen("cut1_crit.dat","w");
//   out2=fopen("cut1_caust.dat","w");
  
  FILE *out3;
  out3=fopen("cut2_crit.dat","w");
//   out3=fopen("cut2_caustt.dat","w");

  FILE *out4;
  out4=fopen("cut3_crit.dat","w");
//   out4=fopen("cut3_caust.dat","w");
    
  FILE *out5;
  out5=fopen("cut4_crit.dat","w");
//   out5=fopen("cut4_caust.dat","w");



  system(" echo  compiling the fortran code, please wait...  ");
  system(" g77-3.4 -Wall -o c_pnfw ../../lensing_models/pnfw/pnfw_contours.f ../../lensing_models/pnfw/root_finding.f ../../lensing_models/pnfw/lensing_functions.f ");
//   system(" gfortran-4.4 -Wall -o c_pnfw ../../lensing_models/pnfw/pnfw_contours.f ../../lensing_models/pnfw/root_finding.f ../../lensing_models/pnfw/lensing_functions.f ");

  system(" echo  ending the compilation  ");

    float argv[3];
    int iargv=1, iflagc=1;
    double lw_u=10.0;
    double rs=1.0;
    int npt=251;
    double elp=0.0;
    double ks;
    double s_el=0.005;
    double s_ks=0.01;
    double cut1_l=0.0009, cut2_l=0.0019, cut3_l=0.0039, cut4_l=0.0059;
    double cut1_u=0.0011, cut2_u=0.0021, cut3_u=0.0041, cut4_u=0.0061;

  for(int j=0; j<=50; j++){
          ks=0.1+j*s_ks;
          printf("the value of ks is %f",ks);
      for(int i=1;i<=199;i++){
          FILE *out1;
          out1=fopen("in_pnfw_par.txt","w");
          elp=i*s_el;
          fprintf(out1,"%f %f %f %f %i %i\n",ks,rs,elp,lw_u,iflagc,npt);
          fclose(out1);    
          system("./c_pnfw");
          system("cp fort.61 cc.dat");
          system("cp fort.71 ca.dat");
          double pert_params[]={elp,1,ks,rs};
          double pot_params[] = {pert_params[2],pert_params[3]};
          double _r_e_nfw= r_e_nfw_find(pot_params);
          double kappa_2=kappa2_nfw(pot_params,_r_e_nfw);
          double d2cc= d2_crit(f1_pnfw, D2f0Dtheta2_pnfw,kappa_2, pert_params, _r_e_nfw);
          if((d2cc > cut1_l) && (d2cc < cut1_u)){fprintf(out2," %f %f %f \n",ks,elp,d2cc);}    
          if((d2cc > cut2_l) && (d2cc < cut2_u)){fprintf(out3," %f %f %f \n",ks,elp,d2cc);}    
          if((d2cc > cut3_l) && (d2cc < cut3_u)){fprintf(out4," %f %f %f \n",ks,elp,d2cc);}    
          if((d2cc > cut4_l) && (d2cc < cut4_u)){fprintf(out5," %f %f %f \n",ks,elp,d2cc);}    
        }
  }

      fclose(out2);
      fclose(out3);
      fclose(out4);
      fclose(out5);  
      system("rm fort.*");
      system("xmgrace -view  0.15 0.15 0.85 0.85 cut1_crit.dat cut2_crit.dat cut3_crit.dat cut4_crit.dat -saveall e_max_crit.agr");
      return 0.0;
}
