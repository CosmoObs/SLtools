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
* Cross section for deformation arcs (... in progress)
*/
/** @package siep_model
*
*  Package to compute the main lensing functions, critical curves, caustics  and constant distortion curves for the  SIEP (Singular Isothermal Elliptical Potential) model 
*
* See Chap.4 in ...sltools/PerturbativeMethod/writeups/Report_on_Perturbative_Method.pdf
*
*  Notice: This model have 3 main parameters.
*
*   siep_params[0]: Integer useful to choose the parameterization for the ellipticity
*
* If siep_params[0]=1, It uses the Angle Deflection Model, i.e., \f$ a_{1\varepsilon}=1-\varepsilon \f$ and \f$ a_{2\varepsilon}=1+\varepsilon \f$
*
* If siep_params[0]=2, It uses the Keeton´s parametrization, i.e., \f$ a_{1\varepsilon}=1\f$ and \f$ a_{2\varepsilon}=(1-\varepsilon)^{-2} \f$
*  
*   siep_params[1]: Ellipticity value (must be ellipticity parameter or the ellipticity by itself whether the angle deflection is used)
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

//! Pseudo-Elliptical Radial Coordinate of the Tangential Critical Curve of the SIEP model
//!
//! \f$ x_{\varepsilon,tcc}= \mathcal{A}(\varepsilon)-\mathcal{B}(\varepsilon)\cos{2\phi_\varepsilon}  \f$, 
//! where \f$ \mathcal{A}(\varepsilon)\f$ and \f$ \mathcal{B}(\varepsilon)\f$ are already defined in siep_lens_funct.h
/*!
//*!
\param theta : angular coordinate, and therefore \f$ \phi_\varepsilon= \arctan{\left(\sqrt{\frac{a_{2\varepsilon}}{a_{1\varepsilon}}}\tan{(\theta)}\right)} \f$
\param siep_params[] : parameters of the SIEP model
  \return \f$ x_{\varepsilon,tcc}(x)\f$

*/
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

//! Pseudo-Elliptical Radial Coordinate of the Constant Distortion curve.
//!
//! \f$ x_{\varepsilon,R_{\rm th}}=x_{\varepsilon,tcc}\times\left\{\begin{array}{cc} \dfrac{|R_{\rm th}|}{|R_{\rm th}|-1}, & R_{\rm th}>0 \\ & \\ //!\dfrac{|R_{\rm th}|}{|R_{\rm th}|+1}, & R_{\rm th}<0 \end{array}\right.\f$, 
/*!
//*
\return \f$ x_{\varepsilon,R_{\rm th}}(x)\f$
\see r_cct_siep
*/
double r_th_siep(double theta, double siep_params[], double rth){
double qth=1.-1/rth;
double r_th_c=r_cct_siep(theta,siep_params)*(1./qth);
return r_th_c;
}

//!  Radial Coordinate of the Radial Critical Curve
//!
//! \f$ x_{\varepsilon,rcc}=0 \f$, 
/*!
//*
\return \f$ x_{\varepsilon,rcc}(x)\f$
*/

double r_ccr_siep(double theta,double siep_params[]){
return 0.0;
}


//! A functions to print the tangential curves ( critical and  caustic) for the SIEP
//!
//! For the tangential critical curve, this function plot the following parametric equation:
//! 
//! \f$ x_{1,tcc}=\dfrac{x_{\varepsilon,tcc}\cos{\phi_\varepsilon}}{\sqrt{a_{1\varepsilon}}}\f$  
//!
//! \f$ x_{2,tcc}=\dfrac{x_{\varepsilon,tcc}\sin{\phi_\varepsilon}}{\sqrt{a_{2\varepsilon}}}\f$
//! 
//! For the tangential caustic, this function plot the following parametric equation:
//!
//! \f$ y_{1,tca} =x_{1,tcc}-\sqrt{a_{1\varepsilon}}\cos{\phi_\varepsilon}\f$ 
//!
//! \f$ y_{2,tca} =x_{2,tcc}-\sqrt{a_{2\varepsilon}}\sin{\phi_\varepsilon} \f$
//!
/*!
\param siep_params[] : SIEP parameters
\param xtcc[] : Vector containing the first coordinate of the tangential critical curve  
\param ytcc[] : Vector containing the second coordinate of the tangential critical curve  
\param npts : number of point to divide the range in values of theta, by default=1000
\return by default the vector (xtcc,ytcc) or if it is indicated  a data file containing such vector 
\see r_cct_siep
*/

void cc_tan_siep(double siep_params[],double xtcc[], double ytcc[], int npts=1000,  FILE *file_in1=NULL){

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

//! A functions to print the constant distortion curve for \f$ R_{\rm th}>0 \f$ for the SIEP
//!
//! In the lens plane  this function plot the following parametric equation:
//! 
//! \f$ x_{1,lwp}=\dfrac{x_{\varepsilon,R_{\rm th}}\cos{\phi_\varepsilon}}{\sqrt{a_{1\varepsilon}}}\f$  
//!
//! \f$ x_{2,lwp}=\dfrac{x_{\varepsilon,R_{\rm th}}\sin{\phi_\varepsilon}}{\sqrt{a_{2\varepsilon}}}\f$
//! 
//! In the source plane, this function plot the following parametric equation:
//!
//! \f$ y_{1,lwp} =x_{1,lwp}-\sqrt{a_{1\varepsilon}}\cos{\phi_\varepsilon}\f$ 
//!
//! \f$ y_{2,lwp} =x_{2,lwp}-\sqrt{a_{2\varepsilon}}\sin{\phi_\varepsilon} \f$
//!
/*!
\param siep_params[] : SIEP parameters
\param xlwp[] : Vector containing the first coordinate of the  constant distortion curve in the lens plane for \f$ R_{\rm th}>0 \f$
\param ylwp[] : Vector containing the second coordinate of the constant distortion curve in the lens plane for \f$ R_{\rm th}>0 \f$
\param npts : number of point to divide the range in values of theta, by default=1000
\return by default the vector (xlwp,ylwp) or if it is indicated  a data file containing such vector 
\see r_th_siep
*/
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
//! A functions to print the constant distortion curve for \f$ R_{\rm th}<0 \f$ for the SIEP
//!
//! In the lens plane  this function plot the following parametric equation:
//! 
//! \f$ x_{1,lwn}=\dfrac{x_{\varepsilon,R_{\rm th}}\cos{\phi_\varepsilon}}{\sqrt{a_{1\varepsilon}}}\f$  
//!
//! \f$ x_{2,lwn}=\dfrac{x_{\varepsilon,R_{\rm th}}\sin{\phi_\varepsilon}}{\sqrt{a_{2\varepsilon}}}\f$
//! 
//! In the source plane, this function plot the following parametric equation:
//!
//! \f$ y_{1,lwn} =x_{1,lwn}-\sqrt{a_{1\varepsilon}}\cos{\phi_\varepsilon}\f$ 
//!
//! \f$ y_{2,lwn} =x_{2,lwn}-\sqrt{a_{2\varepsilon}}\sin{\phi_\varepsilon} \f$
//!
/*!
\param siep_params[] : SIEP parameters
\param xlwn[] : Vector containing the first coordinate of the constant distortion curve in the lens plane for \f$ R_{\rm th}<0 \f$
\param ylwn[] : Vector containing the second coordinate of the constant distortion curve in the lens plane for \f$ R_{\rm th}<0 \f$
\param npts : number of point to divide the range in values of theta, by default=1000
\return by default the vector (xlwn,ylwn) or if it is indicated  a data file containing such vector 
\see r_th_siep
*/
void cc_lw_neg_siep(double siep_params[], double lw_th, double xlwn[], double ylwn[], int npts=1000, FILE *file_in1=NULL){
//
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
//! A functions to print the radial curves ( critical and  caustic) for the SIEP
//!
//! For the radial critical curve, this function plot the following parametric equation:
//! 
//! \f$ x_{1,rcc}=0\f$  
//!
//! \f$ x_{2,rcc}=0\f$
//! 
//! For the radial caustic, this function plot the following parametric equation:
//!
//! \f$ y_{1,rca} =-\sqrt{a_{1\varepsilon}}\cos{\phi_\varepsilon}\f$ 
//!
//! \f$ y_{2,rca} =-\sqrt{a_{2\varepsilon}}\sin{\phi_\varepsilon} \f$
//!
/*!
\param siep_params[] : SIEP parameters
\param xrad[] : Vector containing the first coordinate of the radial critical curve  
\param yrad[] : Vector containing the second coordinate of the radial critical curve  
\param npts : number of point to divide the range in values of theta, by default=1000
\return by default the vector (xrad,yrad) or if it is indicated  a data file containing such vector 
\see r_ccr_siep
*/

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