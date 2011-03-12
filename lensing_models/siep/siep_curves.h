// File: siep_curves.h
// Author: Habib Dumet-Montoya
// Date: January, 2011
/** @file
* In this file, we find the functions that are useful to compute:
*
* Tangential and Radial Critical Curves
*
* Tangential and Radial Caustics
*
* Curves of constant distortion in both planes.
*
* Cross section for deformation arcs (Analytical Expression)
*
* Notice: This model have 3 main parameters.
*
* siep_params[0]: Ellipticity value (must be ellipticity parameter or ellipticity by itself whether the angle deflection model is used)
*
* siep_params[1]: Integer useful to choose the flag for the parameterization for the ellipticity
*
* If siep_params[1]=1 : \f$ a_{1\varepsilon}=1-\varepsilon \f$ and \f$ a_{2\varepsilon}=1+\varepsilon \f$ (Angle Deflection Model)
*
* If siep_params[1]=2 : \f$ a_{1\varepsilon}=1-\varepsilon \f$ and \f$ a_{2\varepsilon}=\frac{1}{1-\varepsilon} \f$ (Standard parametrization)
*
* If siep_params[1]=3 : \f$ a_{1\varepsilon}=1 \f$ and and \f$ a_{2\varepsilon}=\frac{1}{(1-\varepsilon)^{2}} \f$ (Keeton's Parametrization)
*
* siep_params[2]: Mass of the SIS model, related wit the velocity dispersion.*
*
*/

/** @package siep_model
*
*  Package to compute the main lensing functions, critical curves, caustics, constant distortion curves  and cross section for deformation arcs for the SIEP (Singular Isothermal Elliptical Potential) model 
*
* See Chap.4 in ...sltools/PerturbativeMethod/writeups/Report_on_Perturbative_Method.pdf
*
*   Notice: This model have 3 main parameters.
*
* siep_params[0]: Ellipticity value (must be ellipticity parameter or ellipticity by itself whether the angle deflection model is used)
*
* siep_params[1]: Integer useful to choose the flag for the parameterization for the ellipticity
*
* If siep_params[1]=1 : \f$ a_{1\varepsilon}=1-\varepsilon \f$ and \f$ a_{2\varepsilon}=1+\varepsilon \f$ (Angle Deflection Model)
*
* If siep_params[1]=2 : \f$ a_{1\varepsilon}=1-\varepsilon \f$ and \f$ a_{2\varepsilon}=\frac{1}{1-\varepsilon} \f$ (Standard parametrization)
*
* If siep_params[1]=3 : \f$ a_{1\varepsilon}=1 \f$ and  \f$ a_{2\varepsilon}=\frac{1}{(1-\varepsilon)^2} \f$ (Keeton's Parametrization)
*
* siep_params[2]: Mass of the SIS model, related wit the velocity dispersion.
*
*/


#ifndef SIEP_CURVES_H
#define SIEP_CURVES_H

#include <math.h>
//#include <cstdlib>
#include "siep_lens_funct.h" 

#define PI 3.1415926535897932384626433832795 
// #define TAM_MAX 2000

/** Pseudo-Elliptical Radial Coordinate of the Tangential Critical Curve of the SIEP model
*
* \f$ x_{\varepsilon,tcc}= \mathcal{A}(\varepsilon)-\mathcal{B}(\varepsilon)\cos{2\phi_\varepsilon}  \f$, 
* where \f$ \mathcal{A}(\varepsilon)\f$ and \f$ \mathcal{B}(\varepsilon)\f$ are already defined in siep_lens_funct.h
*
* \param theta : angular coordinate, and therefore \f$ \phi_\varepsilon= \arctan{\left(\sqrt{\frac{a_{2\varepsilon}}{a_{1\varepsilon}}}\tan{(\theta)}\right)} \f$
* \param siep_params[] : parameters of the SIEP model
* \return \f$ x_{\varepsilon,tcc}(x)\f$
*/
double r_cct_siep(double theta,double siep_params[]){
double elp=siep_params[0];
int pflag=siep_params[1];
// double b=siep_params[2];
double a1e,a2e,theta_e,Ae,Be;
if(pflag==1){a1e=1.-elp;a2e=1.+elp;}
if(pflag==2){a1e=1.-elp;a2e=1./a1e;}
if(pflag==3){a1e=1.;a2e=1./pow(1.-elp,2);}
theta_e=atan(sqrt(a2e/a1e)*tan(theta));  
Ae=0.5*(a1e+a2e), Be=0.5*(a1e-a2e);
double r_crit=Ae-Be*cos(2.*theta_e);
return r_crit;
}

/** Pseudo-Elliptical Radial Coordinate of the Constant Distortion curve.
*
* \f$ x_{\varepsilon,R_{\rm th}}=x_{\varepsilon,tcc}\times\left\{\begin{array}{cc} \dfrac{|R_{\rm th}|}{|R_{\rm th}|-1}, & R_{\rm th}>0 \\ & \\ \dfrac{|R_{\rm th}|}{|R_{\rm th}|+1}, & R_{\rm th}<0 \end{array}\right.\f$, 
*
*
*\return \f$ x_{\varepsilon,R_{\rm th}}(x)\f$
*\see r_cct_siep
*/

double r_th_siep(double theta, double siep_params[], double rth){
double qth=1.-1/rth;
double r_th_c=r_cct_siep(theta,siep_params)*(1./qth);
return r_th_c;
}

/**  Radial Coordinate of the Radial Critical Curve
*
* \f$ x_{\varepsilon,rcc}=0 \f$, 
*
*
* \return \f$ x_{\varepsilon,rcc}(x)\f$
*/

double r_ccr_siep(double theta,double siep_params[]){
return 0.0;
}


/** A functions to print the tangential curves ( critical and  caustic) for the SIEP
*
* For the tangential critical curve, this function plot the following parametric equation:
* 
* \f$ x_{1,tcc}=\dfrac{x_{\varepsilon,tcc}\cos{\phi_\varepsilon}}{\sqrt{a_{1\varepsilon}}}\f$  
*
* \f$ x_{2,tcc}=\dfrac{x_{\varepsilon,tcc}\sin{\phi_\varepsilon}}{\sqrt{a_{2\varepsilon}}}\f$
* 
* For the tangential caustic, this function plot the following parametric equation:
*
* \f$ y_{1,tca} =x_{1,tcc}-\sqrt{a_{1\varepsilon}}\cos{\phi_\varepsilon}\f$ 
*
* \f$ y_{2,tca} =x_{2,tcc}-\sqrt{a_{2\varepsilon}}\sin{\phi_\varepsilon} \f$
*
*
*\param siep_params[] : SIEP parameters
*\param xtcc[] : Vector containing the first coordinate of the tangential critical curve  
*\param ytcc[] : Vector containing the second coordinate of the tangential critical curve  
*\param npts : number of point to divide the range in values of theta, by default=1000
*\return by default the vector (xtcc,ytcc) or if it is indicated  a data file containing such vector 
*\see r_cct_siep
*/

void cc_tan_siep(double siep_params[],double xtcc[], double ytcc[], int npts=1000,  FILE *file_in1=NULL, FILE *file_in2=NULL ,FILE *file_in3=NULL){

double theta = 0.0, theta_e=0.0;
double x1_cct[npts],x2_cct[npts];
double y1_cct[npts],y2_cct[npts];
double elp=siep_params[0];
int pflag=siep_params[1];
double b=siep_params[2];
double a1e,a2e;
if(pflag==1){a1e=1.-elp;a2e=1.+elp;}
if(pflag==2){a1e=1.-elp;a2e=1./a1e;}
if(pflag==3){a1e=1.;a2e=1./pow(1.-elp,2);}
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
    
    if((file_in2 != NULL) && (file_in3 != NULL)){
      fprintf(file_in2,"%lg %lg \n", x1_cct[i]*b, x2_cct[i]*b);
      fprintf(file_in3,"%lg %lg\n", y1_cct[i]*b, y2_cct[i]*b);  
   }
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

/** A functions to print the constant distortion curve for \f$ R_{\rm th}>0 \f$ for the SIEP
*
* In the lens plane  this function plot the following parametric equation:
* 
* \f$ x_{1,lwp}=\dfrac{x_{\varepsilon,R_{\rm th}}\cos{\phi_\varepsilon}}{\sqrt{a_{1\varepsilon}}}\f$  
*
* \f$ x_{2,lwp}=\dfrac{x_{\varepsilon,R_{\rm th}}\sin{\phi_\varepsilon}}{\sqrt{a_{2\varepsilon}}}\f$
* 
* In the source plane, this function plot the following parametric equation:
*
* \f$ y_{1,lwp} =x_{1,lwp}-\sqrt{a_{1\varepsilon}}\cos{\phi_\varepsilon}\f$ 
*
* \f$ y_{2,lwp} =x_{2,lwp}-\sqrt{a_{2\varepsilon}}\sin{\phi_\varepsilon} \f$
*
*
*\param siep_params[] : SIEP parameters
*\param xlwp[] : Vector containing the first coordinate of the  constant distortion curve in the lens plane for \f$ R_{\rm th}>0 \f$
*\param ylwp[] : Vector containing the second coordinate of the constant distortion curve in the lens plane for \f$ R_{\rm th}>0 \f$
*\param npts : number of point to divide the range in values of theta, by default=1000
*\return by default the vector (xlwp,ylwp) or if it is indicated  a data file containing such vector 
*\see r_th_siep
*/
void cc_lw_pos_siep(double siep_params[], double lw_th, double xlwp[], double ylwp[], int npts=1000, FILE *file_in1=NULL){

double theta = 0.0, theta_e=0.0;
double x1_cct[npts],x2_cct[npts];
double y1_cct[npts],y2_cct[npts];
double elp=siep_params[0];
int pflag=siep_params[1];
double b=siep_params[2];
double a1e,a2e;
if(pflag==1){a1e=1.-elp;a2e=1.+elp;}
if(pflag==2){a1e=1.-elp;a2e=1./a1e;}
if(pflag==3){a1e=1.;a2e=1./pow(1.-elp,2);}

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
/** A functions to print the constant distortion curve for \f$ R_{\rm th}<0 \f$ for the SIEP
*
* In the lens plane  this function plot the following parametric equation:
* 
* \f$ x_{1,lwn}=\dfrac{x_{\varepsilon,R_{\rm th}}\cos{\phi_\varepsilon}}{\sqrt{a_{1\varepsilon}}}\f$  
*
* \f$ x_{2,lwn}=\dfrac{x_{\varepsilon,R_{\rm th}}\sin{\phi_\varepsilon}}{\sqrt{a_{2\varepsilon}}}\f$
* 
* In the source plane, this function plot the following parametric equation:
*
* \f$ y_{1,lwn} =x_{1,lwn}-\sqrt{a_{1\varepsilon}}\cos{\phi_\varepsilon}\f$ 
*
* \f$ y_{2,lwn} =x_{2,lwn}-\sqrt{a_{2\varepsilon}}\sin{\phi_\varepsilon} \f$
*
*\param siep_params[] : SIEP parameters
*\param xlwn[] : Vector containing the first coordinate of the constant distortion curve in the lens plane for \f$ R_{\rm th}<0 \f$
*\param ylwn[] : Vector containing the second coordinate of the constant distortion curve in the lens plane for \f$ R_{\rm th}<0 \f$
*\param npts : number of point to divide the range in values of theta, by default=1000
*\return by default the vector (xlwn,ylwn) or if it is indicated  a data file containing such vector 
*\see r_th_siep
*/
void cc_lw_neg_siep(double siep_params[], double lw_th, double xlwn[], double ylwn[], int npts=1000, FILE *file_in1=NULL){
//
double theta = 0.0, theta_e=0.0;
double x1_cct[npts],x2_cct[npts];
double y1_cct[npts],y2_cct[npts];
double elp=siep_params[0];
int pflag=siep_params[1];
double b=siep_params[2];
double a1e,a2e;
if(pflag==1){a1e=1.-elp;a2e=1.+elp;}
if(pflag==2){a1e=1.-elp;a2e=1./a1e;}
if(pflag==3){a1e=1.;a2e=1./pow(1.-elp,2);}
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
/** A functions to print the radial curves ( critical and  caustic) for the SIEP
*
* For the radial critical curve, this function plot the following parametric equation:
* 
* \f$ x_{1,rcc}=0\f$  
*
* \f$ x_{2,rcc}=0\f$
* 
* For the radial caustic, this function plot the following parametric equation:
*
* \f$ y_{1,rca} =-\sqrt{a_{1\varepsilon}}\cos{\phi_\varepsilon}\f$ 
*
* \f$ y_{2,rca} =-\sqrt{a_{2\varepsilon}}\sin{\phi_\varepsilon} \f$
*
*
* \param siep_params[] : SIEP parameters
* \param xrad[] : Vector containing the first coordinate of the radial critical curve  
* \param yrad[] : Vector containing the second coordinate of the radial critical curve  
* \param npts : number of point to divide the range in values of theta, by default=1000
* \return by default the vector (xrad,yrad) or if it is indicated  a data file containing such vector 
*\see r_ccr_siep
*/

void cc_rad_siep(double siep_params[], double xrad[], double yrad[], int npts=1000, FILE *file_in1=NULL){

double theta = 0.0, theta_e=0.0;
double x1_ccr[npts],x2_ccr[npts];
double y1_ccr[npts],y2_ccr[npts];
double elp=siep_params[0];
int pflag=siep_params[1];
double b=siep_params[2];
double a1e,a2e;
if(pflag==1){a1e=1.-elp;a2e=1.+elp;}
if(pflag==2){a1e=1.-elp;a2e=1./a1e;}
if(pflag==3){a1e=1.;a2e=1./pow(1.-elp,2);}

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
    fprintf(file_in1,"%lg %lg %lg %lg \n", x1_ccr[i]*b, -x2_ccr[i]*b, y1_ccr[i]*b, -y2_ccr[i]*b);}}

}
/** A functions to compute the cross section for arc formation for the SIE model
*
* 
* \f$ \sigma_{\rm siep}= R^2_{\rm E}\frac{|R_{\rm th}|^2+1}{(|R_{\rm th}|^2-1)^2}\times \left\{2\pi\left[ [\mathcal{A}(\varepsilon)-\mathcal{B}(\varepsilon)]^2+\frac{2\,a_{1\varepsilon}\mathcal{B}^2(\varepsilon)}{a_{1\varepsilon}+4\,a_{2\varepsilon}\pi^{2}} \right]+ 2[2\mathcal{A}(\varepsilon)\mathcal{B}(\varepsilon)-\mathcal{B}^2(\varepsilon)]\sqrt{\frac{a_{1\varepsilon}}{a_{2\varepsilon}}}\arctan{\left( 2\pi\sqrt{\frac{a_{2\varepsilon}}{a_{1\varepsilon}}}\right)} \right\} \f$
*
* where \f$ \mathcal{A}(\varepsilon)=\frac{1}{2} (a_{1\varepsilon}+a_{2\varepsilon})\f$ and \f$ \mathcal{B}(\varepsilon)=\frac{1}{2}(a_{1\varepsilon}-a_{2\varepsilon})\f$
*
*
* \param siep_params[] : SIEP parameters
* \param r_th: length-to-width ratio threshold 
* \return \f$ \sigma_{\rm siep} \f$ 
*/
double cross_section_siep(double siep_params[], double r_th){
  double rth_fact=(fabs(r_th)*fabs(r_th)+1.)/(pow(fabs(r_th)*fabs(r_th)-1,2));
  double elp=siep_params[0];
  int pflag=siep_params[1];
  double r_e=siep_params[2];
  double a1e,a2e;
  if(pflag==1){a1e=1.-elp;a2e=1.+elp;}
  if(pflag==2){a1e=1.-elp;a2e=1./a1e;}
  if(pflag==3){a1e=1.;a2e=1./pow(1.-elp,2);}
  double term1,term2,term3;
//   double arg1=pow(10,-8)*sqrt(a2e/a1e);
  double arg2=pow(10,8)*sqrt(a2e/a1e);
//   term1=(0.75*a1e*a1e+0.5*a1e*a2e+0.75*a2e*a2e)*(2.*atan(arg2)-atan(arg1));
//   term2=0.5*(a1e*a1e-a2e*a2e)*(sin(2.*atan(arg1))-2.*sin(2.*atan(arg2)));
//   term3=-0.0625*(a1e-a2e)*(a1e-a2e)*(sin(4.*atan(arg1))-2.*sin(4.*atan(arg2)));
  term1=(0.75*a1e*a1e+0.5*a1e*a2e+0.75*a2e*a2e)*(2.*atan(arg2));
  term2=(a1e*a1e-a2e*a2e)*sin(2.*atan(arg2));
  term3=0.0625*(a1e-a2e)*(a1e-a2e)*2.*sin(4.*atan(arg2));
  double sig_tmp=(1./sqrt(a1e*a2e))*rth_fact*r_e*r_e*(term1+term2+term3);
  return sig_tmp;
}
#endif
