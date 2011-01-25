/** @file
* Example of doxygen documentation for C function FIXME
*/
/** @package siep_model
*  Package to compute critical curves and caustics for the  SIEP (Singular Isothermal Elliptical Potential) model 
*  Notice: This model have 3 main parameters.
*   siep_params[0]: Integer quantity to choose the parameterization for the ellipticity
*  
*   siep_params[1]: Ellipticiy value (must be parameter or ellipticity by itself whether the angle deflection is used)
*
*   siep_params[2]: Mass of the SIS model, related wit the velocity dispersion.
*
*/
#ifndef SIEP_CURVES_H
#define SIEP_CURVES_H

#include <math.h>
#include <cstdlib>
#include "siep_lens_funct.h" 

#define PI 3.1415926535897932384626433832795 
#define TAM_MAX 2000

double r_cct_siep(double theta,double siep_params[]){
int pflag=siep_params[0];
double elp,b,a1e,a2e,theta_e,Ae,Be;
elp=siep_params[1];
b=siep_params[2];
if(pflag==1){a1e=1.-elp;a2e=1.+elp;}
if(pflag==2){a1e=1.;a2e=pow(1.-elp,-2);}
theta_e=atan(sqrt(a2e/a1e)*tan(theta));  
Ae=0.5*(a1e+a2e), Be=0.5*(a1e-a2e);
double r_crit=Ae-Be*cos(2.*theta_e);
return r_crit;
}


double r_th_siep(double theta, double siep_params[], double rth){
double qth=1.-1/rth;
double r_th_c=r_cct_siep(theta,siep_params)*(1./qth);
return r_th_c;
}


double r_ccr_siep(double theta,double siep_params[]){
return 0.0;
}




/** A functions to print the tangential curves ( critical and  caustic) for the SIEP
*/

void cc_tan_siep(double siep_params[],double xtcc[], double ytcc[], int npts=1000,  FILE *file_in1=NULL){


double theta = 0.0, theta_e=0.0;
// double x1_cct[TAM_MAX],x2_cct[TAM_MAX];
// double y1_cct[TAM_MAX],y2_cct[TAM_MAX];
double x1_cct[npts],x2_cct[npts];
double y1_cct[npts],y2_cct[npts];
// double x1_cct,x2_cct,y1_cct,y2_cct;
// double x1_rn, x2_rn, x1_rp, x2_rp;
int pflag=siep_params[0];
double elp=siep_params[1];
double b=siep_params[2];
double a1e,a2e;
if(pflag==1){a1e=1.-elp;a2e=1.+elp;}
if(pflag==2){a1e=1.; a2e=1./pow(1.-elp,2);}


for(int i=0;i<npts;i++){

  theta_e=atan(sqrt(a2e/a1e)*tan(theta));  

  x1_cct[i]=r_cct_siep(theta,siep_params)*cos(theta_e)/sqrt(a1e);
  x2_cct[i]=r_cct_siep(theta,siep_params)*sin(theta_e)/sqrt(a2e);
  xtcc[i]=x1_cct[i];
  ytcc[i]=x2_cct[i];

  y1_cct[i]= fabs(y1_siep(x1_cct[i]/b, x2_cct[i]/b,siep_params));  
  y2_cct[i]= fabs(y2_siep(x1_cct[i]/b, x2_cct[i]/b,siep_params)); 

  if(file_in1==NULL){
      printf("%lg %lg %lg %lg \n",x1_cct[i]*b, x2_cct[i]*b,y1_cct[i]*b, y2_cct[i]*b);
   } else {
     fprintf(file_in1,"%lg %lg %lg %lg \n", x1_cct[i]*b, x2_cct[i]*b,y1_cct[i]*b, y2_cct[i]*b);}
      
   theta+= PI/(2.*npts); 
}
// writting in the second quadrant
for(int i=npts-1;i>0;i--){
     if(file_in1==NULL){
      printf("%lg %lg %lg %lg \n",-x1_cct[i]*b, x2_cct[i]*b,-y1_cct[i]*b, y2_cct[i]*b);
   } else {
    fprintf(file_in1,"%lg %lg %lg %lg \n", -x1_cct[i]*b, x2_cct[i]*b,-y1_cct[i]*b, y2_cct[i]*b);}
}
for(int i=0; i<npts; i++){
    if(file_in1==NULL){
      printf("%lg %lg %lg %lg \n",-x1_cct[i]*b, -x2_cct[i]*b, -y1_cct[i]*b, -y2_cct[i]*b);
   } else {
    fprintf(file_in1," %lg %lg %lg %lg \n", -x1_cct[i]*b, -x2_cct[i]*b, -y1_cct[i]*b, -y2_cct[i]*b);}
    
}
for(int i=npts-1;i>0;i--){
    if(file_in1==NULL){
      printf("%lg %lg %lg %lg \n",x1_cct[i]/b, -x2_cct[i]/b,y1_cct[i]/b, -y2_cct[i]/b);
   } else {
    fprintf(file_in1,"%lg %lg %lg %lg \n", x1_cct[i]*b, -x2_cct[i]*b, y1_cct[i]*b, -y2_cct[i]*b);}
}

}

void cc_lw_pos_siep(double siep_params[], double lw_th, double xlwp[], double ylwp[], int npts=1000, FILE *file_in1=NULL){


double theta = 0.0, theta_e=0.0;
double x1_cct[npts],x2_cct[npts];
double y1_cct[npts],y2_cct[npts];
int pflag=siep_params[0];
double elp=siep_params[1];
double b=siep_params[2];
double a1e,a2e;
if(pflag==1){a1e=1.-elp;a2e=1.+elp;}
if(pflag==2){a1e=1.; a2e=1./pow(1.-elp,2);}


for(int i=0;i<npts;i++){

  theta_e=atan(sqrt(a2e/a1e)*tan(theta));  

  x1_cct[i]=r_th_siep(theta,siep_params,lw_th)*cos(theta_e)/sqrt(a1e);
  x2_cct[i]=r_th_siep(theta,siep_params,lw_th)*sin(theta_e)/sqrt(a2e);
  xlwp[i]=x1_cct[i], ylwp[i]=x2_cct[i];

  y1_cct[i]= fabs(y1_siep(x1_cct[i]/b, x2_cct[i]/b,siep_params));  
  y2_cct[i]= fabs(y2_siep(x1_cct[i]/b, x2_cct[i]/b,siep_params)); 

  if(file_in1==NULL){
      printf("%lg %lg %lg %lg \n",x1_cct[i]*b, x2_cct[i]*b, y1_cct[i]*b, y2_cct[i]*b);
   } else {
     fprintf(file_in1,"%lg %lg %lg %lg \n", x1_cct[i]*b, x2_cct[i]*b,y1_cct[i]*b, y2_cct[i]*b);}
      
   theta+= PI/(2.*npts); 
}
// writting in the second quadrant
for(int i=npts-1;i>0;i--){
    if(file_in1==NULL){
      printf("%lg %lg %lg %lg \n", -x1_cct[i]*b, x2_cct[i]*b, -y1_cct[i]*b, y2_cct[i]*b);
   } else {
    fprintf(file_in1,"%lg %lg %lg %lg \n", -x1_cct[i]*b, x2_cct[i]*b, -y1_cct[i]*b, y2_cct[i]*b);}
}
for(int i=0; i<npts; i++){
    if(file_in1==NULL){
      printf("%lg %lg %lg %lg \n", -x1_cct[i]*b, -x2_cct[i]*b, -y1_cct[i]*b, -y2_cct[i]*b);
   } else {
    fprintf(file_in1,"%lg %lg %lg %lg \n", -x1_cct[i]*b, -x2_cct[i]*b, -y1_cct[i]*b, -y2_cct[i]*b);}
}
for(int i=npts-1;i>0;i--){
    if(file_in1==NULL){
      printf("%lg %lg %lg %lg \n", x1_cct[i]*b, -x2_cct[i]*b, y1_cct[i]*b, -y2_cct[i]*b);
   } else {
    fprintf(file_in1,"%lg %lg %lg %lg \n", x1_cct[i]*b, -x2_cct[i]*b, y1_cct[i]*b, -y2_cct[i]*b);}
}

}

void cc_lw_neg_siep(double siep_params[], double lw_th, double xlwn[], double ylwn[], int npts=1000, FILE *file_in1=NULL){


double theta = 0.0, theta_e=0.0;
double x1_cct[npts],x2_cct[npts];
double y1_cct[npts],y2_cct[npts];
int pflag=siep_params[0];
double elp=siep_params[1];
double b=siep_params[2];
double a1e,a2e;
if(pflag==1){a1e=1.-elp;a2e=1.+elp;}
if(pflag==2){a1e=1.; a2e=1./pow(1.-elp,2);}


for(int i=0;i<npts;i++){

  theta_e=atan(sqrt(a2e/a1e)*tan(theta));  

  x1_cct[i]=r_th_siep(theta,siep_params,-lw_th)*cos(theta_e)/sqrt(a1e);
  x2_cct[i]=r_th_siep(theta,siep_params,-lw_th)*sin(theta_e)/sqrt(a2e);
  xlwn[i]=x1_cct[i], ylwn[i]=x2_cct[i];

  y1_cct[i]= fabs(y1_siep(x1_cct[i]/b, x2_cct[i]/b,siep_params));  
  y2_cct[i]= fabs(y2_siep(x1_cct[i]/b, x2_cct[i]/b,siep_params)); 

  if(file_in1==NULL){
      printf("%lg %lg %lg %lg \n",x1_cct[i]*b, x2_cct[i]*b, y1_cct[i]*b, y2_cct[i]*b);
   } else {
     fprintf(file_in1," %lg %lg %lg %lg \n", x1_cct[i]*b, x2_cct[i]*b, y1_cct[i]*b, y2_cct[i]*b);}
      
   theta+= PI/(2.*npts); 
}
// writting in the second quadrant
for(int i=npts-1;i>0;i--){
    if(file_in1==NULL){
      printf("%lg %lg %lg %lg \n",-x1_cct[i]*b, x2_cct[i]*b,-y1_cct[i]*b, y2_cct[i]*b);
   } else {
    fprintf(file_in1,"%lg %lg %lg %lg  \n", -x1_cct[i]*b, x2_cct[i]*b,-y1_cct[i]*b, y2_cct[i]*b);}
}
for(int i=0; i<npts; i++){
   if(file_in1==NULL){
   printf("%lg %lg %lg %lg \n",-x1_cct[i]*b, -x2_cct[i]*b, -y1_cct[i]*b, -y2_cct[i]*b);
   } else {
    fprintf(file_in1,"%lg %lg %lg %lg \n", -x1_cct[i]*b, -x2_cct[i]*b, -y1_cct[i]*b, -y2_cct[i]*b);}
}
for(int i=npts-1;i>0;i--){
    if(file_in1==NULL){printf("%lg %lg %lg %lg \n",x1_cct[i]*b, -x2_cct[i]*b,y1_cct[i]*b, -y2_cct[i]*b);
   } else {fprintf(file_in1,"%lg %lg %lg %lg \n", x1_cct[i]*b, -x2_cct[i]*b,y1_cct[i]*b, -y2_cct[i]*b);}
}

}



void cc_rad_siep(double siep_params[], double xrad[], double yrad[], int npts=1000, FILE *file_in1=NULL){

double theta = 0.0, theta_e=0.0;
double x1_ccr[npts],x2_ccr[npts];
double y1_ccr[npts],y2_ccr[npts];
int pflag=siep_params[0];
double elp=siep_params[1];
double b=siep_params[2];
double a1e,a2e;
if(pflag==1){a1e=1.-elp;a2e=1.+elp;}
if(pflag==2){a1e=1.; a2e=1./pow(1.-elp,2);}


for(int i=0;i<npts;i++){

  theta_e=atan(sqrt(a2e/a1e)*tan(theta));  

  x1_ccr[i]=r_ccr_siep(theta,siep_params)*cos(theta_e)/sqrt(a1e);
  x2_ccr[i]=r_ccr_siep(theta,siep_params)*sin(theta_e)/sqrt(a2e);
  xrad[i]=x1_ccr[i],yrad[i]=x2_ccr[i];

  y1_ccr[i]= fabs(x1_ccr[i]/b-sqrt(a1e)*cos(theta_e));  
  y2_ccr[i]= fabs(x2_ccr[i]/b-sqrt(a2e)*sin(theta_e)); 


  if(file_in1==NULL){
    printf("%lg %lg %lg %lg \n",x1_ccr[i]*b, x2_ccr[i]*b, y1_ccr[i]*b, y2_ccr[i]*b);
   } else {
  fprintf(file_in1,"%lg %lg %lg %lg \n", x1_ccr[i]*b, x2_ccr[i]*b, y1_ccr[i]*b, y2_ccr[i]*b);}
      
   theta+= PI/(2.*npts); 
}
// writting in the second quadrant
for(int i=npts-1;i>=0;i--){
    if(file_in1==NULL){
    printf("%lg %lg %lg %lg \n", -x1_ccr[i]*b, x2_ccr[i]*b, -y1_ccr[i]*b, y2_ccr[i]*b);
   } else {
    fprintf(file_in1,"%lg %lg %lg %lg \n", -x1_ccr[i]*b, x2_ccr[i]*b, -y1_ccr[i]*b, y2_ccr[i]*b);}
}
//writting in the third quadrant
for(int i=0; i<npts; i++){
    if(file_in1==NULL){
    printf("%lg %lg %lg %lg \n",-x1_ccr[i]*b, -x2_ccr[i]*b, -y1_ccr[i]*b, -y2_ccr[i]*b);
   } else {
    fprintf(file_in1,"%lg %lg %lg %lg \n", -x1_ccr[i]*b, -x2_ccr[i]*b, -y1_ccr[i]*b, -y2_ccr[i]*b);}
}
//writting in the fourth quadrant
for(int i=npts-1;i>=0;i--){
    if(file_in1==NULL){
    printf("%lg %lg %lg %lg \n",x1_ccr[i]*b, -x2_ccr[i]*b, y1_ccr[i]*b, -y2_ccr[i]*b);
   } else {
    fprintf(file_in1,"%lg %lg %lg %lg \n", x1_ccr[i]*b, -x2_ccr[i]*b, y1_ccr[i]*b, -y2_ccr[i]*b);}

}
}

double cross_section_siep(double siep_params[], double lw_ratio, int npts=1000){


return 0.0;
}

#endif