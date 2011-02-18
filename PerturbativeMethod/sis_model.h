/** @file
* Module useful to calculate some function related to the SIS lensing model
*
* - Lensing Function (Angle Deflection, Convergence, Shear)
*
* - Function useful to the Perturbative Approach \f$ \kappa_2 \f$
*
* - Function related with the third derivative of the lensing potential \f$ c_3 \f$ 
* 
*  This module need 
* - pot_params[]={\f$ R_{\rm E}\f$}
*/

*/
/** @package circular_models
*  Package to compute quantities related to the circular SIS model.
*
*
*/
#ifndef SIS_MODEL_H
#define SIS_MODEL_H

#include <math.h>
#include <cstdlib>


//! Convergence of the SIS model.
//!
//! \f$ \kappa(r)= \dfrac{b}{2r}\f$, 
//!
//! where \f$ b \f$ is the parameter characterizing the mass of the SIS, related with the velocity dispersion
/*!   \param r : radial distance,   \param pot_params[] : SIS lens parameters,   \return \f$ \kappa(r)\f$ */
double conv_sis_circ(double r, double pot_params[]){
double b=pot_params[0];
return b/(2.0*r);
}

//! Angle deflection of the SIS model.
//!
//! \f$ \alpha(r)=   b \f$, 
//!  where \f$ b \f$ is defined in the function above (Note: we only considerer the magnitude of the angle deflection)
/*!
  \param r : radial distance,\param pot_params[] : SIS lens parameters 
  \return \f$ \alpha(r)\f$
*/
double alpha_sis_circ(double r, double pot_params[]){
double b=pot_params[0];
return b;
}

//! Shear of the SIS model.
//!
//! \f$ \gamma(r)= \dfrac{b}{2r}\f$, 
//!
//! where \f$ b \f$ is already defined, Note that shear= convergence in this case
/*!   \param r : radial distance,   \param pot_params[] : SIS lens parameters,   \return \f$ \kappa(r)\f$ */
double shear_sis_circ(double r, double pot_params[]){
double b=pot_params[0];
return b/(2.0*r);
}

//! Solution for the Circular NFW Lens 
//!
//! \f$ \kappa_2=2-2\kappa(R_{\mathrm{E}}) \f$,
//!
//! where \f$ \kappa(R_{\mathrm{E}})\f$ is the convergence for evaluated at the Einstein Radius
//!
/*!   \param r_e : Einstein Radius, \param pot_params[] : SIS lens parameters,\return \f$ \kappa_2 \f$ */
double kappa2_sis(double pot_params[],double _r_e_sis){
  double k2,r=_r_e_sis;

  k2=2.0-2.0*conv_sis_circ(r,pot_params);

  return k2;
}

#endif
