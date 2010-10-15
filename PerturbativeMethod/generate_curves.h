/** @file
* Example of doxygem documentation for C functions FIXME
*/

/** @package perturbative_method
* Package to generate the main curves
*
*
* Detailde description FIXME
*
*/
#ifndef GENERATE_CUVES_H
#define GENERATE_CURVES_H

#include "perturbative_method.h"

/** A functions to print the curves model (tangential critical curve,tangential caustic curves and source )
*
*  \param source_in a struct that define elliptical sources
*  \param Df0Dtheta_in function related to the perturbed potential
*  \param f1_in function related to the perturbed potential
*  \param D2f0Dtheta2_in function related to the perturbed potential
*  \param pert_params[] a vector that contains all the perturbation parameters
*  \param kappa2 \f$  \kappa_2 = 1 - \left[\frac{d^2 \phi_0(r)}{dr^2}\right]_{r=r_e} \f$
*  \param npts number of theta divisions 
*  \param _r_e Einstein Radius (If not indicated, the code assumes \f$ R_{_{\mathrm{E}}}=1  \f$ )
*  \param file_in1 a file to write tangential critical curve
*  \param file_in2 a file to write tangential caustic
*  \param file_in3 a file to write the source contour
*  \return nothing
*
*  \sa f_type, r_crit, caustic_y1, caustic_y2, y1_src,y2_src
*/

void plot_curves(elliptical_source source_in, f_type f1_in, f_type Df0Dtheta_in, f_type D2f0Dtheta2_in, double pert_params[], double kappa2, double _r_e, int npts=1000, FILE *file_in1=NULL, FILE *file_in2=NULL, FILE *file_in3=NULL){
  
  double theta = 0.0;
  double r_tcrit=0.0;
  double y1_caust=0.0;
  double y2_caust=0.0;
  double y1s_plot=0.0;
  double y2s_plot=0.0;

  for(int i=0;i<=npts;i++){
  r_tcrit= r_crit(f1_in,D2f0Dtheta2_in,kappa2,theta,pert_params,_r_e); 
  y1_caust=caustic_y1(Df0Dtheta_in,D2f0Dtheta2_in,theta,pert_params,_r_e);
  y2_caust=caustic_y2(Df0Dtheta_in,D2f0Dtheta2_in,theta,pert_params,_r_e);
  y1s_plot=y1_src(source_in, theta);
  y2s_plot=y2_src(source_in, theta);

  if(file_in1==NULL){
        printf("%E %E\n",r_tcrit*cos(theta),r_tcrit*sin(theta));
      } else {
        fprintf(file_in1,"%f %f\n",r_tcrit*cos(theta),r_tcrit*sin(theta));
      }

  if(file_in2==NULL){
        printf("%E %E\n",y1_caust,y2_caust);
      } else {
        fprintf(file_in2,"%f %f\n",y1_caust,y2_caust);
      }

  if(file_in3==NULL){
        printf("%E %E\n",y1s_plot,y2s_plot);
      } else {
        fprintf(file_in3,"%f %f\n",y1s_plot,y2s_plot);
      }

  theta+= 6.283185308/npts;
}

}

#endif