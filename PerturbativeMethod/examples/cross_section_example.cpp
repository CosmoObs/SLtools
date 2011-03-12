#include <cstdlib>
#include <cstdio>
#include <fstream>        // Functions for input and output files
//#include <cstdio>
#include <iostream>


#include "../perturbative_method.h"
//
#include "../siep_model.h"
#include "../sis_model.h"
//
#include "../pnfw_model.h"
#include "../nfw_circular_model.h"
//
#include "../arc_cross_section.h"

// to compile make: g++ -Wall -o exec_dcs cross_section_example.cpp `pkg-config gsl --cflags --libs`


// Falta definir os arquivos da saida para secao de choque ambos modelos
// Escolher que parametros desejamos mudar, colocar um flag para isso!! par_flag[0]=ellipticity, par_flag[1]=mass, par_flag[2]=rth, par_flag[3]=rs
//Documentar estes código

int main(){

// Defining the lenght-to-width ratio threshold

double r_th=10;

// Defining the parameters for the PNFW model  

float argv[3];
int iargv;

      printf("Enter the PNFW parameter as convergence, scale radius,  ellipticity: ");
      scanf("%f %f %f",&argv[0],&argv[1], &argv[2]);

      printf("enter the flag control to choice the parameterization: ");
      scanf("%i",&iargv);

      double pnfw_params[]={argv[2],iargv,argv[0],argv[1]};
      double nfw_params[]={argv[0],argv[1]};
//    
      double _r_e_nfw= r_e_nfw_find(nfw_params);
      printf("Einstein Radius = %f\n",_r_e_nfw);
//
      double kappa_2_nfw=kappa2_nfw(nfw_params,_r_e_nfw); 
//
//       FILE *outcspnfw = fopen ("dcs_pnfw_vs_ks_pm.dat" , "w");
//       FILE *outcspnfw = fopen ("dcs_pnfw_vs_e_pm.dat" , "w");
      FILE *outcspnfw = fopen ("dcs_pnfw_vs_rth_pm.dat" , "w");  
      double sig_pnfw,ks_nfw,epsi;

//           for(int j=0; j<=160; j++){
//             ks_nfw=0.1+0.01*j;
//             pnfw_params[2]=ks_nfw;
//             nfw_params[0]=ks_nfw;
//             _r_e_nfw= r_e_nfw_find(nfw_params);
//             kappa_2_nfw=kappa2_nfw(nfw_params,_r_e_nfw); 
//             sig_pnfw = sigma_defarcs(f1_pnfw, D2f0Dtheta2_pnfw, pnfw_params, kappa_2_nfw,_r_e_nfw, r_th);
//             fprintf(outcspnfw,"%f %f \n", ks_nfw, sig_pnfw);
//           }

//           for(int j=0; j<100; j++){
//             epsi=0.+0.01*j;
//             pnfw_params[0]=epsi;
//             sig_pnfw = sigma_defarcs(f1_pnfw, D2f0Dtheta2_pnfw, pnfw_params, kappa_2_nfw,_r_e_nfw, r_th);
//             fprintf(outcspnfw,"%f %f \n", epsi, sig_pnfw);
//           }

          for(int j=0; j<=100; j++){
            r_th=1.1+(24./100.)*j;
            sig_pnfw = sigma_defarcs(f1_pnfw, D2f0Dtheta2_pnfw, pnfw_params, kappa_2_nfw,_r_e_nfw, r_th);
            fprintf(outcspnfw,"%f %f \n", r_th, sig_pnfw);
          }


// Defining the parameters for the SIEP model

/*    float argu[3];
    int iargu;
      printf("Enter the SIEP parameter as the ellipticity, Einstein Radius: ");
      scanf("%f %f",&argu[0],&argu[1]);
      printf("enter the flag control to choice the parameterization: ");
      scanf("%i",&iargu);
//
      double siep_params[]={argu[0],iargu,argu[1]};
      double sis_params[]={argu[1]};
//
       double _r_e_sis= sis_params[0];
       printf("Einstein Radius = %f\n",_r_e_sis);
//
       double kappa_2_sis = kappa2_sis(sis_params,_r_e_sis);*/
//
//       double sig_siep = sigma_defarcs(f1_siep, D2f0Dtheta2_siep, siep_params, kappa_2_sis,_r_e_sis, rth);
//       printf("cross-sec=%f\n",sig_siep);

//       FILE *outcsda = fopen ("dcs_siep_vs_re_pm.dat" , "w");      
//       FILE *outcsda = fopen ("dcs_siep_vs_e_pm.dat" , "w");
//      FILE *outcsda = fopen ("dcs_siep_vs_rth_pm.dat" , "w");
//         double sig_siep,re_sis;
//         double sis_params[1];
/*
      for(int j=0; j<=1000; j++){
          re_sis=0.01+0.01*j;
          siep_params[2]=re_sis;
          sis_params[0]=re_sis;
          kappa_2_sis = kappa2_sis(sis_params,re_sis);
          sig_siep = sigma_defarcs(f1_siep, D2f0Dtheta2_siep, siep_params, kappa_2_sis,re_sis, r_th);
          fprintf(outcsda,"%f %f \n", re_sis, sig_siep);
      }
*/
/*
      for(int j=0; j<100; j++){
        epsi=0.+0.01*j;
        siep_params[0]=epsi;
        sig_siep = sigma_defarcs(f1_siep, D2f0Dtheta2_siep, siep_params, kappa_2_sis,_r_e_sis, r_th);
        fprintf(outcsda,"%f %f \n", epsi, sig_siep);
      }
            
*/

// /*      double r_th;
//       for(int j=0; j<=2000; j++){
//           r_th=1.1+0.01*j;
//           sig_siep = sigma_defarcs(f1_siep, D2f0Dtheta2_siep, siep_params, kappa_2_sis,_r_e_sis, r_th);
//           fprintf(outcsda,"%f %f \n", r_th, sig_siep);
//       }

 //   fclose(outcsda);

}
