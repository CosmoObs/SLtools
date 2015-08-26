/** @file
* This module computes:
* - Main fields of the PNFW in the Perturbative Approach \f$ f_{n}(\theta)\f$ and their derivatives
*
* This module needs:  
* - pert_params[]={\f$ \eta \f$ ,flag, \f$ \kappa_s \f$, \f$ r_s \f$}, where flag is a integer to choice the parametrization of the ellipticity
* - Einstein Radius of the NFW model, i.e.,  \f$ R_{\rm E}\f$ 
*/

/** @package pseudo_elliptical_models
*  Package to compute quantities related to PNFW model.
*
*/

#ifndef PNFW_MODEL_H
#define PNFW_MODEL_H

double f1_pnfw(double theta, double pert_params[],double _r_e_nfw);


double Df0Dtheta_pnfw( double theta, double pert_params[],double _r_e_nfw);

double D2f0Dtheta2_pnfw( double theta, double pert_params[],double _r_e_nfw);

#endif
