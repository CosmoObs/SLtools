#include<cstdlib>
#include<cstdio>
#include <fstream>        //funcoes de entrada e saida para arquivos
#include <iostream>

#include "../perturbative_method.h"
//
#include "../siep_model.h"
#include "../sis_model.h"
//
#include "../comparison_fns.h"
//
#include "../../lensing_models/siep/siep_curves.h"

// #define TAM_MAX 2000

// to compile g++ -Wall -o e_fom_siep fom_siep_test.cpp


int main(){
              int npts=251;
              double elp=0.0;
              double s_el=0.005;
//
              double x1tan[2000], x2tan[2000];
//
              FILE *outfom=fopen("fom_siep.dat","w");
//
              float argu[1];
              int iargu=1;
              printf("Enter Einstein Radius: ");
              scanf("%f", &argu[1]);
      
              for(int i=1;i<=199;i++){
	 	 FILE *out=fopen("cc_siep.dat","w");
		  FILE *out1=fopen("cc.dat","w");
		  FILE *out2=fopen("ca.dat","w");
		  elp=i*s_el;
		  double siep_params[]={elp,iargu,argu[1]};
		  double sis_params[]={argu[1]};
// Calculating the caustic and critical curves in the tangential direction for the SIEP	
		  cc_tan_siep(siep_params, x1tan,x2tan,npts,out,out1,out2);
//      
		  fclose(out);
	  	fclose(out1);
		  fclose(out2);
//	  
		  double _r_e_sis= sis_params[0];
	  	printf("Einstein Radius = %f\n",_r_e_sis);
		  double kappa_2_sis = kappa2_sis(sis_params,_r_e_sis);
		  double mwsrfd_cc= d2_crit(f1_siep,D2f0Dtheta2_siep,kappa_2_sis, siep_params, _r_e_sis);
		  double mwsrfd_ca= d2_caust(f1_siep,Df0Dtheta_siep,D2f0Dtheta2_siep, kappa_2_sis, siep_params,_r_e_sis);
		  fprintf(outfom,"%lg %lg %lg \n",elp,mwsrfd_cc,mwsrfd_ca);
              }

           }


