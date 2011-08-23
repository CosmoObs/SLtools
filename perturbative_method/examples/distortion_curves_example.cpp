#include <cstdlib>
#include <cstdio>
#include <fstream>        // Functions for input and output files
//#include <cstdio>
#include <iostream>


#include "../perturbative_method.h"

#include "../siep_model.h"
#include "../sis_model.h"

#include "../pnfw_model.h"
#include "../nfw_circular_model.h"

#include "../generate_curves.h"


// #include "perturbative_method.h"
// //
// #include "siep_model.h"
// #include "sis_model.h"
// //
// #include "pnfw_model.h"
// #include "nfw_circular_model.h"
// //
// #include "generate_curves.h"

// to compile make: g++ -Wall -o exec_cdc distortion_curves_example.cpp `pkg-config gsl --cflags --libs`


// Falta definir os arquivos da saida para ambos modelos
// usar o batch file del xmgrace para fazer os gráficos de forma automática
//Documentar estes código
int main(){

// Defining the lenght-to-width ratio threshold

double rth=5;
int npts=1000;

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
      FILE *outpnfw1 = fopen ("pnfw-tang_curves.dat" , "w");
      FILE *outpnfw2 = fopen ("pnfw-curves_r_th_pos.dat" , "w");
      FILE *outpnfw3 = fopen ("pnfw-curves_r_th_neg.dat" , "w");

      plot_curves(f1_pnfw, Df0Dtheta_pnfw, D2f0Dtheta2_pnfw, pnfw_params, kappa_2_nfw,_r_e_nfw, npts,outpnfw1);   
      plot_rth_curves(f1_pnfw, Df0Dtheta_pnfw, D2f0Dtheta2_pnfw, pnfw_params, kappa_2_nfw, _r_e_nfw, 10., npts,outpnfw2,outpnfw3);   

// Defining the parameters for the SIEP model
    float argu[3];
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
      double kappa_2_sis = kappa2_sis(sis_params,_r_e_sis);
//
    
      plot_curves(f1_siep, Df0Dtheta_siep, D2f0Dtheta2_siep, siep_params, kappa_2_sis,_r_e_sis, npts);
      plot_rth_curves(f1_siep, Df0Dtheta_siep, D2f0Dtheta2_siep, siep_params, kappa_2_sis,_r_e_sis, rth, npts);   

}
