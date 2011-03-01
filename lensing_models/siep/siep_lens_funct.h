// File: siep_lens_funct.h
// Author: Habib Dumet-Montoya
// Date: January, 2011
/** @file
* This module is useful to compute quantities related to the Singular Isothermal Elliptical Potential (SIEP) 
*
* See sltools/PerturbativeMethod/writeups/Report_on_Perturbative_Method.pdf
*
* Lensing functions to be considered: angle deflection, convergence, components of the shear.
*
* Lens Equation 
*
*  Notice: This model have 3 main parameters.
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
* siep_params[2]: Mass of the SIS model, related wit the velocity dispersion.
*
*/

/** @package siep_model
*
*
*/

#ifndef SIEP_MODEL_H
#define SIEP_MODEL_H

#include <math.h>
#include <cstdlib>


//! Module of the Angle deflection of the SIS model.
//!
//! \f$ \alpha(r)= 1 \f$, 
/*!
  \return \f$ \alpha(r)\f$
*/
double alpha_sis(double siep_params[]){
return 1.;
}

//! Convergence of the SIS model.
//!
//! \f$ \kappa(r)= \frac{1}{2x}\f$, 
//!
//! \f$ x=b/r \f$, where \f$ b \f$ is the parameter characterizing the mass of the SIS, related with the velocity dispersion
/*!   \param r : radial distance,   \param b : SIS lens parameters,   \return \f$ \kappa(r)\f$ */

double kappa_sis(double r, double siep_params[]){
// double b=siep_params[2];
return 1./(2*r);
}

//! Shear of the SIS model.
//!
//! \f$ \gamma(r)= \frac{1}{2x}\f$, 
//!
//! \f$ x=b/r \f$, where \f$ b \f$ is the parameter characterizing the mass of the SIS, related with the velocity dispersion
/*!   \param r : radial distance,   \param b : SIS lens parameters,    \return \f$ \kappa(r)\f$ */

double gamma_sis(double r, double siep_params[]){
return kappa_sis(r,siep_params);
}

// This part corresponds to pseudo elliptical functions
//! Convergence of the SIEP model
//!
//! \f$ \kappa_\varepsilon(\vec{x})=\mathcal{A}(\varepsilon)\kappa(x_\varepsilon)-\mathcal{B}(\varepsilon)\gamma(x_\varepsilon)\cos{2\phi_\varepsilon} \f$
//!
//! where \f$ \mathcal{A}(\varepsilon)= \frac{1}{2}(a_{1\varepsilon}+a_{2\varepsilon})\f$, \f$ \mathcal{B}(\varepsilon)=\frac{1}{2}(a_{1\varepsilon}-a_{2\varepsilon})\f$ and \f$\phi_\varepsilon=\arctan \left(\frac{\sqrt{a_{2\varepsilon}}x_{2}}{\sqrt{a_{1\varepsilon}}x_{1}}\right)\f$
/*!
\param  xi1,xi2 : are the cartesian coordinates
\param siep_params[] : SIEP parameters  
\return \f$ \kappa_\varepsilon(x)\f$
*/
double kappa_siep(double xi1, double xi2, double siep_params[]){
double elp=siep_params[0];
int pflag=siep_params[1];
double b=siep_params[2];
double x1,x2,a1e,a2e;
x1=b*xi1, x2=b*xi2;
if(pflag==1){a1e=1.-elp;a2e=1.+elp;}
if(pflag==2){a1e=1.-elp;a2e=1./a1e;}
if(pflag==3){a1e=1.;a2e=1./pow(1.-elp,2);}
double theta_e,xe,Ae,Be;
theta_e=atan(sqrt(a2e/a1e)*(x2/x1));  
xe=sqrt(a1e*pow(x1,2)+a2e*pow(x2,2));
Ae=0.5*(a1e+a2e);
Be=0.5*(a1e-a2e);
double kappa=Ae*kappa_sis(xe,siep_params)-Be*gamma_sis(xe,siep_params)*cos((2.*theta_e));
return kappa;
}

//! First component of the shear of the SIEP model
//!
//! \f$ \gamma_{1\varepsilon}(\vec{x}) = \mathcal{B}(\varepsilon)\kappa(x_\varepsilon)-\mathcal{A}(\varepsilon)\gamma(x_\varepsilon)\cos{2\phi_\varepsilon} \f$
//!
//! where \f$ \mathcal{A}(\varepsilon)= \frac{1}{2}(a_{1\varepsilon}+a_{2\varepsilon})\f$, \f$ \mathcal{B}(\varepsilon)=\frac{1}{2}(a_{1\varepsilon}-a_{2\varepsilon})\f$ and \f$\phi_\varepsilon=\arctan \left(\frac{\sqrt{a_{2\varepsilon}}x_{2}}{\sqrt{a_{1\varepsilon}}x_{1}}\right)\f$
/*!
\param xi1,xi2: cartesian coordinates
\param siep_params[] : SIEP parameters
\return \return \f$ \gamma_{1\varepsilon(x)}\f$
*/
double gamma1_siep(double xi1, double xi2, double siep_params[]){
double elp=siep_params[0];
int pflag=siep_params[1];
double b=siep_params[2];
double x1=b*xi1, x2=b*xi2;
double a1e, a2e;
if(pflag==1){a1e=1.-elp;a2e=1.+elp;}
if(pflag==2){a1e=1.-elp;a2e=1./a1e;}
if(pflag==3){a1e=1.;a2e=1./pow(1.-elp,2);}
double theta_e=atan(sqrt(a2e/a1e)*(x2/x1));  
double xe=sqrt(a1e*pow(x1,2)+a2e*pow(x2,2));
double Ae=0.5*(a1e+a2e);
double Be=0.5*(a1e-a2e);
double gama1=Be*kappa_sis(xe,siep_params)-Ae*gamma_sis(xe,siep_params)*cos(2.*theta_e);
return gama1;
}

//! Second component of the shear of the SIEP model (see report)
//!
//! \f$ \gamma_{2\varepsilon}(\vec{x})=-\sqrt{\mathcal{A}^2(\varepsilon)-\mathcal{B}^2(\varepsilon)}\gamma(x_\varepsilon)\sin{ 2\phi_\varepsilon} \f$
//!
//! where \f$ \mathcal{A}(\varepsilon)= \frac{1}{2}(a_{1\varepsilon}+a_{2\varepsilon})\f$, \f$ \mathcal{B}(\varepsilon)=\frac{1}{2}(a_{1\varepsilon}-a_{2\varepsilon})\f$ and \f$\phi_\varepsilon=\arctan \left(\frac{\sqrt{a_{2\varepsilon}}x_{2}}{\sqrt{a_{1\varepsilon}}x_{1}}\right)\f$
/*!
\param xi1,xi2 : cartesian coordinates
\param siep_params[] : SIEP parameters
\return \return \f$ \gamma_{2\varepsilon(x)}\f$
*/
double gamma2_siep(double xi1, double xi2, double siep_params[]){
double elp=siep_params[0];
int pflag=siep_params[1];
double b=siep_params[2];
double x1=b*xi1, x2=b*xi2;
double a1e, a2e;
if(pflag==1){a1e=1.-elp;a2e=1.+elp;}
if(pflag==2){a1e=1.-elp;a2e=1./a1e;}
if(pflag==3){a1e=1.;a2e=1./pow(1.-elp,2);}
double theta_e=atan(sqrt(a2e/a1e)*(x2/x1));  
double xe=sqrt(a1e*pow(x1,2)+a2e*pow(x2,2));
double Ae=0.5*(a1e+a2e);
double Be=0.5*(a1e-a2e);
double gama2=-sqrt(pow(Ae,2)-pow(Be,2))*gamma_sis(xe,siep_params)*sin((2.*theta_e));
return gama2;
}

//! Shear (modulus) of the SIEP model
//!
//! \f$ \gamma_\varepsilon(\vec{x})=\mathcal{A}(\varepsilon)\kappa(x_\varepsilon)-\mathcal{B}(\varepsilon)\gamma(x_\varepsilon)\cos{2\phi_\varepsilon} \f$
//!
//! where \f$ \mathcal{A}(\varepsilon)= \frac{1}{2}(a_{1\varepsilon}+a_{2\varepsilon})\f$, \f$ \mathcal{B}(\varepsilon)=\frac{1}{2}(a_{1\varepsilon}-a_{2\varepsilon})\f$ and \f$ \phi_\varepsilon=\arctan \left(\frac{\sqrt{a_{2\varepsilon}}x_{2}}{\sqrt{a_{1\varepsilon}}x_{1}}\right)\f$
/*!
\param  xi1,xi2 : are the cartesian coordinates
\param siep_params[] : SIEP parameters  
\return \f$ \gamma_\varepsilon(x)\f$
*/
double gamma_siep(double xi1,double xi2, double siep_params[]){
double gamma=kappa_siep(xi1,xi2, siep_params);
return gamma;
}
//! First component of the deflection angle of the SIEP model
//!
//! \f$ \alpha_{1\varepsilon} = \alpha(x_\varepsilon)\sqrt{a_{1\varepsilon}} \cos { \phi_ \varepsilon} \f$
//!
//! where \f$ \phi_\varepsilon=\arctan \left(\frac{\sqrt{a_{2\varepsilon}}x_{2}}{\sqrt{a_{1\varepsilon}}x_{1}}\right)\f$
//!
/*!
\param  xi1,xi2 : are the cartesian coordinates
\param siep_params[] : SIEP parameters  
\return \f$ \alpha_{1,\varepsilon}(x)\f$
*/

double alpha1_siep(double xi1, double xi2, double siep_params[]){
double elp=siep_params[0];
int pflag=siep_params[1];
double b=siep_params[2];
double x1,x2,a1e,a2e;
x1=b*xi1, x2=b*xi2;
if(pflag==1){a1e=1.-elp;a2e=1.+elp;}
if(pflag==2){a1e=1.-elp;a2e=1./a1e;}
if(pflag==3){a1e=1.;a2e=1./pow(1.-elp,2);}
double theta_e=atan(sqrt(a2e/a1e)*(x2/x1));  
double alpha1=alpha_sis(siep_params)*sqrt(a1e)*cos(theta_e);
return alpha1;

}

//! Second component of the deflection angle of the SIEP model
//!
//! \f$ \alpha_{2\varepsilon} = \alpha(x_\varepsilon)\sqrt{a_{2\varepsilon}} \sin{\phi_ \varepsilon} \f$
//!
//! where \f$ \phi_\varepsilon=\arctan \left(\frac{\sqrt{a_{2\varepsilon}}x_{2}}{\sqrt{a_{1\varepsilon}}x_{1}}\right)\f$
//!
/*!
\param  xi1,xi2 : are the cartesian coordinates
\param siep_params[] : SIEP parameters  
\return \f$ \alpha_{2,\varepsilon}(x)\f$
*/

double alpha2_siep(double xi1, double xi2, double siep_params[]){
double elp=siep_params[0];
int pflag=siep_params[1];
double b=siep_params[2];
double x1,x2,a1e,a2e;
double theta_e,alpha2;
elp=siep_params[0];
b=siep_params[2];
x1=b*xi1, x2=b*xi2;
if(pflag==1){a1e=1.-elp;a2e=1.+elp;}
if(pflag==2){a1e=1.-elp;a2e=1./a1e;}
if(pflag==3){a1e=1.;a2e=1./pow(1.-elp,2);}
theta_e=atan(sqrt(a2e/a1e)*(x2/x1));  
alpha2=alpha_sis(siep_params)*sqrt(a2e)*sin(theta_e);
return alpha2;

}
//! First component of the lens equation for the SIEP model
//!
//! \f$  y_1(x)=\frac{x_\varepsilon\cos{\phi_\varepsilon}}{\sqrt{a_{1\varepsilon}}}-\sqrt{a_{1\varepsilon}}\cos{\phi_ \varepsilon}\f$
//!
//! where \f$ \phi_\varepsilon=\arctan \left(\frac{\sqrt{a_{2\varepsilon}}x_{2}}{\sqrt{a_{1\varepsilon}}x_{1}}\right)\f$
//!
/*!
\param  xi1,xi2 : are the cartesian coordinates
\param siep_params[] : SIEP parameters  
\return \f$ y_{1}(x)\f$
*/

double y1_siep(double xi1, double xi2, double siep_params[]){
double b=siep_params[2];
double x1=b*xi1;
double y1=x1-alpha1_siep(xi1,xi2,siep_params);
return y1;
}

//! Second component of the lens equation for the SIEP model
//!
//! \f$  y_2(x)=\frac{x_\varepsilon\sin{\phi_\varepsilon}}{\sqrt{a_{2\varepsilon}}}-\sqrt{a_{2\varepsilon}}\sin{\phi_ \varepsilon}\f$
//!
//! where \f$ \phi_\varepsilon=\arctan \left(\frac{\sqrt{a_{2\varepsilon}}x_{2}}{\sqrt{a_{1\varepsilon}}x_{1}}\right)\f$
//!
//!
/*!
\param  xi1,xi2 : are the cartesian coordinates
\param siep_params[] : SIEP parameters  
\return \f$ y_{2}(x)\f$
*/
double y2_siep(double xi1, double xi2, double siep_params[]){
double b=siep_params[2];
double x2=b*xi2;
double y2=x2-alpha2_siep(xi1,xi2,siep_params);
return y2;
}

#endif
