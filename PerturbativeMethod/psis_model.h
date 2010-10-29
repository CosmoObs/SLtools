/** @file
* Example of doxygen documentation for C function FIXME
*/
/** @package pnfw_model
*  Package to compute quantities related to Pseudo-Elliptical SIS model.
*
*  Detailed descrition FIXME
*
*/

#ifndef PSIS_MODEL_H
#define PSIS_MODEL_H

#include "sis_model.h"
#include "pseudo_elliptical_model.h"

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
double f1_psis(double theta, double pert_params[],double _r_e_sis){
    double f1=f1_pe(alpha_sis_circ, _r_e_sis, theta,pert_params);
    return f1;
}
//!  Second perturbative field for the PE-SIS model as Perturbation (see PerturbedPNFW.pdf) 
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
  \param pert_params[] : SIS lens parameters
  \return \f$ \dfrac{df_0}{d\theta}\f$ */

double Df0Dtheta_psis( double theta, double pert_params[],double _r_e_sis){
  return df0dtheta_pe(alpha_sis_circ, _r_e_sis, theta, pert_params);
}
//!  Second perturbative field for the PE-SIS model as Perturbation
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
  \param pert_params[] : PNFW lens parameters
  \return \f$ \dfrac{df_0}{d\theta}\f$ */

double D2f0Dtheta2_psis( double theta, double pert_params[],double _r_e_sis){
  return d2f0dtheta2_pe(alpha_sis_circ,shear_sis_circ, _r_e_sis, theta, pert_params);
}
#endif
