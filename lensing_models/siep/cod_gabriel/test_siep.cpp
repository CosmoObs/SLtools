#include <cstdlib>
#include <cstdio>
#include <fstream>        // Functions for input and output files
//#include <cstdio>
#include <iostream>

// #include "siep_lens_funct.h" 
#include "siep_curves.h" 

#define PAR_MAX1 5

// to compile make: g++ -Wall -o exec_siep test_siep.cpp

int main(){
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
      double ratio=5.0;

      FILE *outtc = fopen ("tang_crit.dat" , "w");
      FILE *outcau = fopen ("tang_caust.dat" , "w");
      
      FILE *outrc = fopen ("rad_crit.dat" , "w");
      FILE *outcar = fopen ("rad_caust.dat" , "w");


      FILE *ollwp = fopen ("lw_pos_lens.dat" , "w");
      FILE *oslwp = fopen ("lw_pos_src.dat" , "w");


      FILE *ollwn = fopen ("lw_neg_lens.dat" , "w");
      FILE *oslwn = fopen ("lw_neg_src.dat" , "w");
      
//       cc_tan_siep(siep_params, npts, ratio, outtc, outcau);

      cc_tan_siep(siep_params, npts, outtc,outcau);
      cc_rad_siep(siep_params, npts, outrc,outcar);
      cc_lw_pos_siep(siep_params,ratio,npts,ollwp,oslwp);
      cc_lw_neg_siep(siep_params,ratio,npts,ollwn,oslwn);
      fclose(outtc); 
      fclose(outcau); 

      fclose(outrc); 
      fclose(outcar); 

      fclose(ollwp);
      fclose(oslwp);

      fclose(ollwn);
      fclose(oslwn);
      
//       system("xmgrace -view 0.15 0.15 0.85 0.85 tang_crit.dat  -saveall curvas.agr ");
      system("xmgrace -view 0.15 0.15 0.85 0.85 tang_crit.dat rad_crit.dat lw_pos_lens.dat lw_neg_lens.dat -saveall curvas.agr ");
      system("xmgrace -view 0.15 0.15 0.85 0.85 tang_caust.dat rad_caust.dat lw_pos_src.dat lw_neg_src.dat -saveall curvas2.agr ");

}