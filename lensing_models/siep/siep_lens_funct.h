/** @file
* Example of doxygen documentation for C function FIXME
*/
/** @package siep_model
*  Package to compute quantities related to the circular SIEP (Singular Isothermal Elliptical Potential) model for.
*  dimensionles coordinates.
*  Notice: This model have 3 main parameters.
*
*   siep_params[0]: Integer quantity to choose flag for the parameterization for the ellipticity
*  
*   siep_params[1]: Ellipticiy value (must be parameter or ellipticity by itself whether the angle deflection is used)
*
*   siep_params[2]: Mass of the SIS model, related wit the velocity dispersion.
*  Detailed descrition FIXME
*
*/
#ifndef SIEP_MODEL_H
#define SIEP_MODEL_H

#include <math.h>
#include <cstdlib>


//! Angle deflection of the SIS model.
//!
//! \f$ \alpha(r)=   b \f$, 
//!  where \f$ b \f$ is defined in the function above (Note: we only considerer the magnitude of the angle deflection)
/*!
  \param b : quantity related with the mass of the SIS.
  \return \f$ \alpha(r)\f$
*/
double alpha_sis(double siep_params[]){
// double b=siep_params[2];
return 1.;
}

//! Convergence of the SIS model.
//!
//! \f$ \kappa(r)= \dfrac{b}{2r}\f$, 
//!
//! where \f$ b \f$ is the parameter characterizing the mass of the SIS, related with the velocity dispersion
/*!   \param r : radial distance,   \param b : SIS lens parameters,   \return \f$ \kappa(r)\f$ */

double kappa_sis(double r, double siep_params[]){
// double b=siep_params[2];
return 1./(2*r);
}

//! Shear of the SIS model.
//!
//! \f$ \gamma(r)= \dfrac{b}{2r}\f$, 
//!
//! where \f$ b \f$ is already defined, Note that shear= convergence in this case
/*!   \param r : radial distance,   \param pot_params[] : SIS lens parameters,   \return \f$ \kappa(r)\f$ */

double gamma_sis(double r, double siep_params[]){
return kappa_sis(r,siep_params);
}

// This part corresponds to pseudo elliptical functions
//! Convergence of the SIEP model (see my paper )
double kappa_siep(double xi1, double xi2, double siep_params[]){
int pflag=siep_params[0];
double elp,b,x1,x2,a1e,a2e;
elp=siep_params[1];
b=siep_params[2];
x1=b*xi1, x2=b*xi2;
if(pflag==1){a1e=1.-elp;a2e=1.+elp;}
if(pflag==2){a1e=1.;a2e=1./pow(1.-elp,2);}
double theta_e,xe,Ae,Be;
theta_e=atan(sqrt(a2e/a1e)*(x2/x1));  
xe=sqrt(a1e*pow(x1,2)+a2e*pow(x2,2));
Ae=0.5*(a1e+a2e);
Be=0.5*(a1e-a2e);
double kappa=Ae*kappa_sis(xe,siep_params)-Be*gamma_sis(xe,siep_params)*cos((2.*theta_e));
return kappa;
}

//! First component of the shear of the SIEP model (see my paper )
double gamma1_siep(double xi1, double xi2, double siep_params[]){
double elp=siep_params[0];
int pflag=siep_params[1];
double b=siep_params[2];
double x1=b*xi1, x2=b*xi2;
double a1e, a2e;
if(pflag==1){a1e=1.-elp;a2e=1.+elp;}
if(pflag==2){a1e=1.;a2e=1./pow(1.-elp,2);}
double theta_e=atan(sqrt(a2e/a1e)*(x2/x1));  
double xe=sqrt(a1e*pow(x1,2)+a2e*pow(x2,2));
double Ae=0.5*(a1e+a2e);
double Be=0.5*(a1e-a2e);
double gama1=Be*kappa_sis(xe,siep_params)-Ae*gamma_sis(xe,siep_params)*cos(2.*theta_e);
return gama1;
}

//! Second component of the shear of the SIEP model (see my paper )
double gamma2_siep(double xi1, double xi2, double siep_params[]){
double elp=siep_params[0];
int pflag=siep_params[1];
double b=siep_params[2];
double x1=b*xi1, x2=b*xi2;
double a1e, a2e;
if(pflag==1){a1e=1.-elp;a2e=1.+elp;}
if(pflag==2){a1e=1.;a2e=1./pow(1.-elp,2);}
double theta_e=atan(sqrt(a2e/a1e)*(x2/x1));  
double xe=sqrt(a1e*pow(x1,2)+a2e*pow(x2,2));
double Ae=0.5*(a1e+a2e);
double Be=0.5*(a1e-a2e);
double gama2=-sqrt(pow(Ae,2)-pow(Be,2))*gamma_sis(xe,siep_params)*sin((2.*theta_e));
return gama2;
}

//! Shear (modulus) of the SIEP model (see my paper )
double gamma_siep(double xi1,double xi2, double siep_params[]){
double gamma=kappa_siep(xi1,xi2, siep_params);
return gamma;
}
//! First component of the deflection angle of the SIEP model

double alpha1_siep(double xi1, double xi2, double siep_params[]){
int pflag;
double elp,b,x1,x2,a1e,a2e;
pflag=siep_params[0];
elp=siep_params[1];
b=siep_params[2];
x1=b*xi1, x2=b*xi2;
if(pflag==1){a1e=1.-elp;a2e=1.+elp;}
if(pflag==2){a1e=1.; a2e=1./pow(1.-elp,2);}
double theta_e=atan(sqrt(a2e/a1e)*(x2/x1));  
double alpha1=alpha_sis(siep_params)*sqrt(a1e)*cos(theta_e);
return alpha1;

}

//! Second component of the deflection angle of the SIEP model
double alpha2_siep(double xi1, double xi2, double siep_params[]){
int pflag;
pflag=siep_params[0];
double elp,b,x1,x2,a1e,a2e;
double theta_e,alpha2;
elp=siep_params[1];
b=siep_params[2];
x1=b*xi1, x2=b*xi2;
if(pflag==1){a1e=1.-elp;a2e=1.+elp;}
if(pflag==2){a1e=1.;a2e=1./pow(1.-elp,2);}
theta_e=atan(sqrt(a2e/a1e)*(x2/x1));  
alpha2=alpha_sis(siep_params)*sqrt(a2e)*sin(theta_e);
return alpha2;

}
//! First component of the lens equation for the SIEP model (see my paper )
double y1_siep(double xi1, double xi2, double siep_params[]){
double b=siep_params[2];
double x1=b*xi1;
double y1=x1-alpha1_siep(xi1,xi2,siep_params);
return y1;
}

//! Second component of the lens equation for the SIEP model (see my paper )
double y2_siep(double xi1, double xi2, double siep_params[]){
double b=siep_params[2];
double x2=b*xi2;
double y2=x2-alpha2_siep(xi1,xi2,siep_params);
return y2;
}




#endif