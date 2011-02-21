// File: siep_example.cpp
// Author: Habib Dumet-Montoya
// Date: January, 2011
/** @file
* This code compute:
*
* Critical Curves and Caustics for
*
* Tangential Direction, it uses "cc_tan_siep"
*
* Radial Direction, it uses it uses "cc_tan_siep"  
*
* Curves of constant distortion in both planes. It uses "cc_lw_pos_siep" and "cc_lw_neg_siep"
*
* to compile do: g++ -Wall -o exec_siep siep_example.cpp
*
* to run the program, there is two options:
*
* the first option: we need a text file "input_lens.txt" containing 
*
* iflag (integer to choice the parameterization, can be 1 or 2)
*
* ellipticity value
*
* Einstein Radius value
*
* In the second option, in terminal to put the values of iflag, ellipticity value, Einstein Radius Value 
*
*/
/** @package siep_model
*
**/


//! The main program
/*!
\param argv[] : SIEP parameters 
  \return  vector containing cartesian coordinates of the curves, show the plot in the screen and the data files of the curves if requested by the person.
*/

#include <cstdlib>
#include <cstdio>
#include <fstream>        
#include <iostream>
#include "siep_curves.h" 
#define PAR_MAX1 5
#define TAM_MAX 2000



int main(){
  double x1tan[TAM_MAX], x2tan[TAM_MAX];
  double x1lwp[TAM_MAX], x2lwp[TAM_MAX];
  double x1lwn[TAM_MAX], x2lwn[TAM_MAX];
  double x1rad[TAM_MAX], x2rad[TAM_MAX];  
  
  FILE *entrada1;
  entrada1=fopen("input_lens.txt","r");
  float argv[3];
  float lens_par[PAR_MAX1];
//   int iargv=1;

int i=0;
  if (entrada1== NULL) 
    {
      printf("Nao foi possivel abrir o arquivo... then \n");
      printf("Enter the flag, ellipticity and Einstein radius: ");
      scanf("%f %f %f",&argv[0],&argv[1], &argv[2]);
      lens_par[1]=argv[0];
      lens_par[2]=argv[1];
      lens_par[3]=argv[2];

//       printf("enter the flag control to choice the parameterization: ");
//       scanf("%i",&iargv);
//       break;
//        exit();
     }else{
      while (! feof(entrada1))
      {
        i++;
        fscanf(entrada1," %f",&lens_par[i]);
      }}

      double siep_params[]={lens_par[1],lens_par[2],lens_par[3]};
      int npts = 251;
      double ratio=10.0;

      FILE *outtc = fopen ("tang_curves.dat" , "w");
      
      FILE *outrc = fopen ("rad_curves.dat" , "w");


      FILE *ollwp = fopen ("lw_pos_curves.dat" , "w");


      FILE *ollwn = fopen ("lw_neg_curves.dat" , "w");
      
//       cc_tan_siep(siep_params, npts, ratio, outtc, outcau);

      cc_tan_siep(siep_params, x1tan,x2tan,npts, outtc);
      cc_rad_siep(siep_params,x1rad,x2rad, npts,  outrc);
      cc_lw_pos_siep(siep_params,ratio, x1lwp,x2lwp,npts, ollwp);
      cc_lw_neg_siep(siep_params,ratio,x1lwn,x2lwn,npts,ollwn);
      fclose(outtc); 

      fclose(outrc); 

      fclose(ollwp);

      fclose(ollwn);
      
//       system("xmgrace -view 0.15 0.15 0.85 0.85 tang_curves.dat rad_curves.dat lw_pos_curves.dat lw_neg_curves.dat -saveall curvas.agr ");
//       system("xmgrace -view 0.15 0.15 0.85 0.85 -block tang_curves.dat -bxy 3:4 -block rad_curves.dat -bxy 3:4  -block lw_pos_curves.dat -bxy 3:4 -block  lw_neg_curves.dat -bxy 3:4 -saveall curvas2.agr ");
      
//       FILE *outcsda = fopen ("dcs_siep_vs_e.dat" , "w");
//       FILE *outcsda = fopen ("dcs_siep_vs_re.dat" , "w");
      FILE *outcsda = fopen ("dcs_siep_vs_rth.dat" , "w");
      double sig1,re_sis;
 /*
      for(int j=0; j<=1000; j++){
          re_sis=0.+0.01*j;
          siep_params[2]=re_sis;
          sig1=cross_section_siep(siep_params,ratio);
          fprintf(outcsda,"%f %f \n", re_sis, sig1);
      }

      for(int j=0; j<100; j++){
        re_sis=0.+0.01*j;
        siep_params[1]=re_sis;
        sig1=cross_section_siep(siep_params,ratio);
        fprintf(outcsda,"%f %f \n", re_sis, sig1);
      }
            
*/	
      double r_th;
      for(int j=0; j<=2000; j++){
          r_th=1.1+0.01*j;
          sig1=cross_section_siep(siep_params,r_th);
          fprintf(outcsda,"%f %f \n", r_th, sig1);
      }

    fclose(outcsda);
    
}