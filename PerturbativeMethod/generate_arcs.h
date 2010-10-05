/** @file
* Example of doxygen documentation for C functions FIXME. 
*/

/** @package perturbative_method
*  Package to generate the arcs
*
*  Detailed descrition FIXME
*
*/


#ifndef GENERATE_ARCS_H
#define GENERATE_ARCS_H

#include "perturbative_method.h"

/** A functions to print the arcs
*
*  \param source_in a struct that define elliptical sources
*  \param Df0Dtheta_in function related to the perturbed potential
*  \param f1_in function related to the perturbed potential
*  \param Df0Dtheta_in function related to the perturbed potential
*  \param pert_params[] a vector that contains all the perturbation parameters
*  \param kappa2 \f$  \kappa_2 = 1 - \left[\frac{d^2 \phi_0(r)}{dr^2}\right]_{r=r_e} \f$
*  \param npts number of theta divisions 
*  \param file_in a file to write the arc's contours
*  \return nothing
*
*  \sa f_type, dr_plus, dr_minus
*/
void plot_arcs(elliptical_source source_in, f_type f1_in, f_type Df0Dtheta_in, double pert_params[], double kappa2, int npts=500, FILE *file_in=NULL){
  double theta = 0.0;
  double r_p = 0.0;
  double r_m = 0.0;

  for(int i=0;i<=npts;i++){

    if(arg_sqrt(Df0Dtheta_in, source_in, theta, pert_params)>0.0){
      r_p = _r_e+dr_plus(f1_in, Df0Dtheta_in, kappa2, source_in, theta, pert_params);
      r_m = _r_e+dr_minus(f1_in, Df0Dtheta_in, kappa2, source_in, theta, pert_params);
      if(file_in==NULL){
        printf("%E %E\n",r_p*cos(theta),r_p*sin(theta));
        printf("%E %E\n",r_m*cos(theta),r_m*sin(theta));
      } else {
        fprintf(file_in,"%E %E\n",r_p*cos(theta),r_p*sin(theta));
        fprintf(file_in,"%E %E\n",r_m*cos(theta),r_m*sin(theta));
      }
    }
    theta+= 6.283185308/npts;
  }

}


#endif
