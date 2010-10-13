/** @file
* Example of doxygen documentation for C functions FIXME. 
*/

/** @package nfw_model
*  Package to compute quantities related to PNFW model.
*
*  Detailed descrition FIXME
*
*/
#ifndef PNFW_MODEL_H
#define PNFW_MODEL_H

#include "nfw_circular_model.h"

//!  First perturbative field for the PNFW model as Perturbation.
//!
//!  \f$ \psi_E(r)=\phi_E(\xi)-\phi_0(r) \f$ 
//!
//!  \f$ \xi=r\sqrt{a_{1\eta}\cos^2{\theta}+a_{2\eta}\sin^2{\theta}}\f$ 
//!        (\f$ \xi_E\f$ is the elliptical variable \f$ \xi(r,\theta)\f$ for \f$ r=R_{E}\f$) 
//!
//!  \f$ f_{1}(\theta) \equiv \left[\dfrac{\partial\psi_E(r,\theta)}{\partial r}\right]_{r=R_{E}} = \dfrac{\xi_{\mathrm{E}}}{R_{\mathrm{E}}}\alpha(\xi_{\mathrm{E}})-\alpha(R_{\mathrm{E}})\f$
/*!
  \param theta : angular coordinate in the lens plane,
  \param pot_params[] : PNFW lens parameters
  \return \f$ f_1(\theta)\f$
*/
double f1_pnfw(double theta, double pot_params[],double _r_e_nfw){
    double eta=pot_params[2];
//     double a_eta1=1.0-eta,a_eta2=1.0/(1.0-eta); // If you want work with the standar parameterization
//     double a_eta1=1.0-eta,a_eta2=1.0+eta;   // If you want work with the Angle Deflection Method
    double a_eta1=1.0,a_eta2=1.0/pow((1.0-eta),2); // If you want work with the Keeton's parameterization
    double xie=_r_e_nfw*sqrt(a_eta1*pow(cos(theta),2)+a_eta2*pow(sin(theta),2));
    double f1=(xie/_r_e_nfw)*alpha_nfw_circ(xie,pot_params)-alpha_nfw_circ(_r_e_nfw,pot_params);
    return f1;
}



//!  Second perturbative field for the PNFW model as Perturbation (see PerturbedPNFW.pdf) 
//!
//!  \f$ \psi_E(r)=\phi_E(\xi)-\phi_0(r) \f$
//!
//!  \f$ \xi=r\sqrt{a_{1\eta}\cos^2{\theta}+a_{2\eta}\sin^2{\theta}}\f$
//!      (\f$ \xi_\mathrm{E}\f$ is the elliptical variable \f$ \xi(r,\theta)\f$ for \f$ r=R_{\mathrm{E}} \f$
//!
//!  \f$ \dfrac{df_0}{d\theta}= \left[\dfrac{d \psi_E(r,\theta)}{d\theta}\right]_{r=R_{\mathrm{E}}} = \dfrac{R^2_{\mathrm{E}}}{2}\dfrac{\alpha(\xi_{\mathrm{E}})}{\xi_{\mathrm{E}}}\mathcal{A}(\eta)\sin{2\theta},\quad \mathcal{A}(\eta)=a_{2\eta}-a_{1\eta} \f$
//!
/*!
  \param theta : angular coordinate in the lens plane
  \param pot_params[] : PNFW lens parameters
  \return \f$ \dfrac{df_0}{d\theta}\f$ */

double Df0Dtheta_pnfw( double theta, double pot_params[],double _r_e_nfw){
    double eta=pot_params[2];
//     double a_eta1=1.0-eta,a_eta2=1.0/(1.0-eta); // If you want work with the standar parameterization
//     double a_eta1=1.0-eta,a_eta2=1.0+eta;    // // If you want work with the Angle Deflection Method
    double a_eta1=1.0,a_eta2=1.0/pow((1.0-eta),2); // If you want work with the Keeton's parameterization
    double A_eta=a_eta2-a_eta1;
    double xie=_r_e_nfw*sqrt(a_eta1*pow(cos(theta),2)+a_eta2*pow(sin(theta),2));
    double cal_g=pow(_r_e_nfw,2)*A_eta*sin(2*theta);
   double df0dte=0.5*(alpha_nfw_circ(xie,pot_params)/xie)*cal_g;
   return df0dte; 
}



//!  Solution for the PNFW model as perturbation,useful for critical and caustic lines
//!
//!  \f$ \psi_E(r)=\phi_E(\xi)-\phi_0(r) \f$
//!
//!  \f$ \xi=r\sqrt{a_{1\eta}\cos^2{\theta}+a_{2\eta}\sin^2{\theta}}\f$
//!                  (\f$ \xi_\mathrm{E}\f$ is the elliptical variable \f$ \xi(r,\theta)\f$ for \f$ r=R_{\mathrm{E}})\f$ 
//!
//!  \f$ \dfrac{d^2f_0}{d\theta^2}= \left[\dfrac{d^2 \psi_E(r,\theta)}{d\theta^2}\right]_{r=R_{\mathrm{E}}}= \mathcal{G}(\eta)\dfrac{\alpha(\xi_{\mathrm{E}})}{\xi_{\mathrm{E}}}-\mathcal{H}(\eta)\gamma(\xi_{\mathrm{E}})\f$
//! 
//!  \f$ \mathcal{G}(\eta)=\mathcal{A}(\eta)r^2_{\mathrm{E}}\cos{(2\theta)}\f$
//!
//!  \f$ \mathcal{H}(\eta)=\dfrac{r^2_{\mathrm{E}}}{2}\left\{\mathcal{A}(\eta)\dfrac{R_{\mathrm{E}}}{\xi_{\mathrm{E}}}\sin{(2\theta)}\right\}^2 \f$
//!
//! \f$ \gamma(\xi_{\mathrm{E}})=\kappa(\xi_{\mathrm{E}})-\dfrac{\alpha(\xi_{\mathrm{E}})}{\xi_{\mathrm{E}}} \f$
/*!
  \param theta : angular coordinate in the lens plane
  \param r : radial distance
  \param pot_params[] : PNFW lens parameters
  \return \f$ \dfrac{d^2f_0}{d\theta^2}\f$
*/

double D2f0Dtheta2_pnfw( double theta, double pot_params[],double _r_e_nfw){
    double eta=pot_params[2];
//     double a_eta1=1.0-eta,a_eta2=1.0/(1.0-eta); // If you are interested in work with the standard parametrization
//     double a_eta1=1.0-eta,a_eta2=1.0+eta;    // If you are interested in work with the Angle Deflection Model
    double a_eta1=1.0,a_eta2=1.0/pow((1.0-eta),2); // If you want work with the Keeton's parameterization
    double A_eta=a_eta2-a_eta1;
    double xie=_r_e_nfw*sqrt(a_eta1*pow(cos(theta),2)+a_eta2*pow(sin(theta),2));
    double cal_g=pow(_r_e_nfw,2)*A_eta*sin(2*theta);
    double cal_gsq=pow(cal_g,2)/_r_e_nfw;
    double dcal_g=2.0*pow(_r_e_nfw,2)*A_eta*cos(2*theta);
    
    double d2f0dte2=0.5*dcal_g*(alpha_nfw_circ(xie,pot_params)/xie)-0.5*shear_nfw_circ(xie,pot_params)*cal_gsq;
    return d2f0dte2;
}


#endif
