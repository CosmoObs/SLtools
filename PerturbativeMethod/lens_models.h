/** @file
* Example of doxygen documentation for C functions FIXME. 
*/

/** @package lens_models
*  Package to compute quantities related to several lens models, perturbations and radial simetric models.
*
*  Detailed descrition FIXME
*
*/



#ifndef LENSMODELS_H
#define LENSMODELS_H

/**********************************************************************************************************************/
/**  \f$ f_1(\theta) \f$ for SIS model
*
*  \param theta angular coordinate in the lens plane
*  \param pert_params[] a NULL vector
*  \return 0.0

*  \sa Df0Dtheta_SIS, f_type
*/
double f1_SIS(double theta, double pert_params[]){
  return 0.0;
}

/**  \f$ \frac{d f_0(\theta)}{d \theta} \f$ for SIS model
*
*  \param theta angular coordinate in the lens plane
*  \param pert_params[] a NULL vector
*  \return 0.0

*  \sa f1_SIS, f_type
*/
double Df0Dtheta_SIS(double theta, double pert_params[]){
  return 0.0;
}
/**********************************************************************************************************************/


/**********************************************************************************************************************/
//pot_params[0] = mp
//pot_params[1] = rp
//pot_params[2] = thetap
/**  \f$ f_1(\theta) \f$ for SIS as a perturbation
*
*  \param theta angular coordinate in the lens plane
*  \param pot_params[0] = mp
*  \param pot_params[1] = rp
*  \param pot_params[2] = thetap
*  \return \f$ f_1(\theta) \f$

*  \sa Df0Dtheta_pert_SIS, f_type
*/
double f1_pert_SIS(double theta, double pot_params[]){
  double tmp1 = pot_params[0] * (1.0 - pot_params[1]*cos(theta-pot_params[2]) );
  double tmp2 = sqrt ( 1.0 - 2.0*pot_params[1]*cos(theta-pot_params[2]) + pot_params[1]*pot_params[1]);

  return tmp1/tmp2 ;
}

/**  \f$ \frac{d f_0(\theta)}{d \theta} \f$ for SIS as a perturbation
*
*  \param theta angular coordinate in the lens plane
*  \param pot_params[0] = mp
*  \param pot_params[1] = rp
*  \param pot_params[2] = thetap
*  \return \f$ \frac{d f_0(\theta)}{d \theta} \f$

*  \sa f1_pert_SIS, f_type
*/
double Df0Dtheta_pert_SIS(double theta, double pot_params[]){
  double tmp1 = pot_params[0] * (pot_params[1]*sin(theta-pot_params[2]) );
  double tmp2 = sqrt ( 1.0 - 2.0*pot_params[1]*cos(theta-pot_params[2]) + pot_params[1]*pot_params[1]);
  return tmp1/tmp2;
}
/**********************************************************************************************************************/



#endif
