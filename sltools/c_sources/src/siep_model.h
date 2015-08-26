/** @file
* This module computes:
* - Main fields of the SIEP in the Perturbative Approach \f$ f_{n}(\theta)\f$ and their derivatives
*
* This module needs:  
* - pert_params[]={\f$ \eta \f$ ,flag, \f$ R_{\rm E} \f$}, where flag is a integer to choice the parametrization of the ellipticity
* - Einstein Radius of the SIS model, i.e.,  \f$ R_{\rm E}\f$ 
*/
/** @package pseudo_elliptical_models
*  Package to compute quantities related to Singular Isothermal Elliptical Potential (SIEP) model.
*
*/

#ifndef PSIS_MODEL_H
#define PSIS_MODEL_H


double f1_siep(double theta, double pert_params[],double _r_e_sis);

double Df0Dtheta_siep( double theta, double pert_params[],double _r_e_sis);

double D2f0Dtheta2_siep( double theta, double pert_params[],double _r_e_sis);

#endif
