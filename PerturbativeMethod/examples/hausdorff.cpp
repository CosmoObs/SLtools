#include <cstdio>
#include <cmath>
#include <algorithm>
#include <cstdlib>

#include "../../numerical_methods/hausdorff_distance.h"




int main(){

  FILE * h_in = fopen ("../../lensing_models/pnfw/in_pnfw_par.txt","w");
  double ks = 1.2;
  double rs = 1.0;
  double elp= 0.1;
  double lw_u = 1.0;
  int iflag=1;
  int npt = 20;
  fprintf(h_in,"%f %f %f %f %i %i",ks, rs, elp, lw_u, iflag, npt); 
  fclose(h_in);

  system("sh hausdorff.sh");


//Read the output from fortran program and put in cc_x[] and cc_y[]
  FILE *cc_file = fopen("../../lensing_models/pnfw/fort.65","r");
  FILE *ca_file = fopen("../../lensing_models/pnfw/fort.75","r");
  int number_columns = 500;  //number of coluns (characteres)
  char buffer[number_columns]; //variavel a receber as linhas do arquivo
  double cc_x[4*npt], cc_y[4*npt];
  double ca_x[4*npt], ca_y[4*npt];


  int i=0;
  while (!feof(cc_file)){
    fgets (buffer , number_columns , cc_file);
    sscanf(buffer ,"%lg %lg", &cc_x[i] ,&cc_y[i]);
    fgets (buffer , number_columns , ca_file);
    sscanf(buffer ,"%lg %lg", &ca_x[i] ,&ca_y[i]);
    printf("%i %f %f %f %f \n",i,cc_x[i],cc_y[i],ca_x[i],ca_y[i]);
    i++;
  }

  fclose(cc_file);
  fclose(ca_file);





/*  int npts1 = 1000; //number of poits for ellipse 1
  int npts2 = 1000; //number of poits for ellipse 2
  double c1x[npts1],c1y[npts1];//vector for curve 1
  double c2x[npts2],c2y[npts2];//vector for curve 2
  double a1,e1; //major semi axis and ellipticity for ellipse 1
  double a2,e2; //major semi axis and ellipticity for ellipse 2
  double x02,y02;//center of the ellipse 2

  a1=2.0;
  e1 = 0.;

  a2 =2.25;
  e2 = 0.5;
  x02 = 0.1;
  y02 = 0.2;

  double t=0;
  for(int i=0;i<npts1;i++){//defining curve 1
    c1x[i] = a1*cos(t);
    c1y[i] = a1*sqrt(1.0-e1)*sin(t);
    t+=6.28319/double(npts1);
  }
  t=0.0;
  for(int i=0;i<npts2;i++){//defining curve 2
    c2x[i] = x02 + a2*cos(t);
    c2y[i] = y02 + a2*sqrt(1.0-e2)*sin(t);
    t+=6.28319/double(npts2);
  }

*/
/********************************************************************************************************************
  printf("Distancies from curve A to curve B\n");
  printf("  d1    : %f\n", dist2curves(c1x,c1y,npts1,c2x,c2y,npts2,1) );
  printf("  d2    : %f\n", dist2curves(c1x,c1y,npts1,c2x,c2y,npts2,2) );
  printf("  d3    : %f\n", dist2curves(c1x,c1y,npts1,c2x,c2y,npts2,0) );
  printf("  d4    : %f\n", dist2curves_max(c1x,c1y,npts1,c2x,c2y,npts2) );
  printf("  d3/d2 : %f\n", dist2curves(c1x,c1y,npts1,c2x,c2y,npts2,0)/dist2curves(c1x,c1y,npts1,c2x,c2y,npts2,2) );
  printf("  d3/d4 : %f\n", dist2curves(c1x,c1y,npts1,c2x,c2y,npts2,0)/dist2curves_max(c1x,c1y,npts1,c2x,c2y,npts2) );

  printf("Distancies from curve B to curve A\n");
  printf("  d1    : %f\n", dist2curves(c2x,c2y,npts2,c1x,c1y,npts1,1) );
  printf("  d2    : %f\n", dist2curves(c2x,c2y,npts2,c1x,c1y,npts1,2) );
  printf("  d3    : %f\n", dist2curves(c2x,c2y,npts2,c1x,c1y,npts1,0) );
  printf("  d4    : %f\n", dist2curves_max(c1x,c1y,npts1,c2x,c2y,npts2) );
  printf("  d3/d2 : %f\n", dist2curves(c1x,c1y,npts1,c2x,c2y,npts2,0)/dist2curves(c1x,c1y,npts1,c2x,c2y,npts2,2) );
  printf("  d3/d4 : %f\n", dist2curves(c1x,c1y,npts1,c2x,c2y,npts2,0)/dist2curves_max(c1x,c1y,npts1,c2x,c2y,npts2) );



  FILE * c1File = fopen ("c1.dat","w");
  for(int i=0;i<npts1;i++){ fprintf(c1File,"%f  %f\n",c1x[i], c1y[i]); }
  FILE * c2File = fopen ("c2.dat","w");
  for(int i=0;i<npts2;i++){ fprintf(c2File,"%f  %f\n",c2x[i], c2y[i]); }
  fclose(c1File);
  fclose(c2File);
  system("xmgrace c1.dat c2.dat &");
********************************************************************************************************************/






/********************************************************************************************************************
  double factor = 2.0;
  printf("\n\nEnlarging curves by a factor %f\n",factor);
  for(int i=0;i<npts1;i++){ c1x[i]=factor*c1x[i]; c1y[i]=factor*c1y[i]; }
  for(int i=0;i<npts2;i++){ c2x[i]=factor*c2x[i]; c2y[i]=factor*c2y[i]; }
  printf("Distancies from curve A to curve B\n");
  printf("  d1    : %f\n", dist2curves(c1x,c1y,npts1,c2x,c2y,npts2,1) );
  printf("  d2    : %f\n", dist2curves(c1x,c1y,npts1,c2x,c2y,npts2,2) );
  printf("  d3    : %f\n", dist2curves(c1x,c1y,npts1,c2x,c2y,npts2,0) );
  printf("  d4    : %f\n", dist2curves_max(c1x,c1y,npts1,c2x,c2y,npts2) );
  printf("  d3/d2 : %f\n", dist2curves(c1x,c1y,npts1,c2x,c2y,npts2,0)/dist2curves(c1x,c1y,npts1,c2x,c2y,npts2,2) );
  printf("  d3/d4 : %f\n", dist2curves(c1x,c1y,npts1,c2x,c2y,npts2,0)/dist2curves_max(c1x,c1y,npts1,c2x,c2y,npts2) );

  printf("Distancies from curve B to curve A\n");
  printf("  d1    : %f\n", dist2curves(c2x,c2y,npts2,c1x,c1y,npts1,1) );
  printf("  d2    : %f\n", dist2curves(c2x,c2y,npts2,c1x,c1y,npts1,2) );
  printf("  d3    : %f\n", dist2curves(c2x,c2y,npts2,c1x,c1y,npts1,0) );
  printf("  d4    : %f\n", dist2curves_max(c1x,c1y,npts1,c2x,c2y,npts2) );
  printf("  d3/d2 : %f\n", dist2curves(c1x,c1y,npts1,c2x,c2y,npts2,0)/dist2curves(c1x,c1y,npts1,c2x,c2y,npts2,2) );
  printf("  d3/d4 : %f\n", dist2curves(c1x,c1y,npts1,c2x,c2y,npts2,0)/dist2curves_max(c1x,c1y,npts1,c2x,c2y,npts2) );
  FILE * c12File = fopen ("c1.dat","w");
  for(int i=0;i<npts1;i++){ fprintf(c12File,"%f  %f\n",c1x[i], c1y[i]); }
  FILE * c22File = fopen ("c2.dat","w");
  for(int i=0;i<npts2;i++){ fprintf(c22File,"%f  %f\n",c2x[i], c2y[i]); }
  fclose(c12File);
  fclose(c22File);
  system("xmgrace c1.dat c2.dat &");
********************************************************************************************************************/

  return 0;
}
