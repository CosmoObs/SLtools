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


/** @package circular_models
*  Package to compute quantities related to the circular SIS model.
*
*
*/
#ifndef SIS_MODEL_H
#define SIS_MODEL_H

#include <math.h>
#include <stdlib.h>


double conv_sis_circ(double r, double pot_params[]);

double conv_sis_circ_prime(double r, double pot_params[]);

double alpha_sis_circ(double r, double pot_params[]);

double shear_sis_circ(double r, double pot_params[]);

double kappa2_sis(double pot_params[],double _r_e_sis);

#endif
