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
#include "theta_find.h"


/** A functions to print the arcs
*
*  \param source_in a struct that define elliptical sources
*  \param Df0Dtheta_in function related to the perturbed potential
*  \param f1_in function related to the perturbed potential
*  \param pert_params[] a vector that contains all the perturbation parameters
*  \param kappa2 \f$  \kappa_2 = 1 - \left[\frac{d^2 \phi_0(r)}{dr^2}\right]_{r=r_e} \f$
*  \param _r_e Einstein Radius (if is not write, the code assumes \f$ R_{_{\mathrm{E}}}=1 \f$ )
*  \param npts number of theta divisions 
*  \param file_in a file to write the arc's contours
*  \return nothing
*
*  \sa f_type, dr_plus, dr_minus
*/
void plot_arcs(elliptical_source source_in, f_type f1_in, f_type Df0Dtheta_in, double pert_params[], double kappa2, double _r_e, int npts=500, FILE *file_in=NULL){
  double theta = 0.0;
  double r_p = 0.0;
  double r_m = 0.0;

  double arcs[10][2];
  theta_find(Df0Dtheta_in, source_in,  pert_params, _r_e, arcs);

  for(int i=0;i<=npts;i++){

    if(arg_sqrt(Df0Dtheta_in, source_in, theta, pert_params, _r_e)>0.0){
      r_p = _r_e+dr_plus(f1_in, Df0Dtheta_in, kappa2, source_in, theta, pert_params,_r_e);
      //r_m = _r_e+dr_minus(f1_in, Df0Dtheta_in, kappa2, source_in, theta, pert_params,_r_e);
      if(file_in==NULL){
        printf("%E %E\n",r_p*cos(theta),r_p*sin(theta));
        //printf("%E %E\n",r_m*cos(theta),r_m*sin(theta));
      } else {
        fprintf(file_in,"%E %E\n",r_p*cos(theta),r_p*sin(theta));
        //fprintf(file_in,"%E %E\n",r_m*cos(theta),r_m*sin(theta));
      }
    }
    theta+= 6.283185308/npts;
  }

  for(int i=0;i<=npts+1;i++){

    if(arg_sqrt(Df0Dtheta_in, source_in, theta, pert_params, _r_e)>0.0){
      //r_p = _r_e+dr_plus(f1_in, Df0Dtheta_in, kappa2, source_in, theta, pert_params,_r_e);
      r_m = _r_e+dr_minus(f1_in, Df0Dtheta_in, kappa2, source_in, theta, pert_params,_r_e);
      if(file_in==NULL){
        //printf("%E %E\n",r_p*cos(theta),r_p*sin(theta));
        printf("%E %E\n",r_m*cos(theta),r_m*sin(theta));
      } else {
        //fprintf(file_in,"%E %E\n",r_p*cos(theta),r_p*sin(theta));
        fprintf(file_in,"%E %E\n",r_m*cos(theta),r_m*sin(theta));
      }
    }
    theta-= 6.283185308/npts;
  }





}


/** A functions to print the arcs separately
*
*  \param source_in a struct that define elliptical sources
*  \param Df0Dtheta_in function related to the perturbed potential
*  \param f1_in function related to the perturbed potential
*  \param pert_params[] a vector that contains all the perturbation parameters
*  \param kappa2 \f$  \kappa_2 = 1 - \left[\frac{d^2 \phi_0(r)}{dr^2}\right]_{r=r_e} \f$
*  \param _r_e Einstein Radius (if is not write, the code assumes \f$ R_{_{\mathrm{E}}}=1 \f$ )
*  \param npts number of theta divisions 
*  \param file_in a file to write the arc's contours
*  \return nothing
*
*  \sa f_type, dr_plus, dr_minus
*/
void plot_arcs_sep(elliptical_source source_in, f_type f1_in, f_type Df0Dtheta_in, double pert_params[], double kappa2, double _r_e, int npts=500, FILE *file_in=NULL){
  double theta = 0.0;
  double r_p = 0.0;
  double r_m = 0.0;

  double arcs[10][2];
  theta_find(Df0Dtheta_in, source_in,  pert_params, _r_e, arcs);

//  printf("\%E %E\n",arcs[0][0],arcs[0][1]);
//  printf("\%E %E\n",arcs[1][0],arcs[1][1]);
//  printf("\%E %E\n",arcs[2][0],arcs[2][1]);
//  printf("\%E %E\n",arcs[3][0],arcs[3][1]);

  if(arcs[ int(arcs[0][0]) ][1]< arcs[ int(arcs[0][0]) ][0]) arcs[ int(arcs[0][0]) ][1] += 6.283185308;
  if(arcs[ int(arcs[0][0]+1) ][1]< arcs[ int(arcs[0][0]+1) ][0]) arcs[ int(arcs[0][0]+1) ][1] += 6.283185308;

  printf("\%E %E\n",arcs[0][0],arcs[0][1]);
  printf("\%E %E\n",arcs[1][0],arcs[1][1]);
  printf("\%E %E\n",arcs[2][0],arcs[2][1]);
  printf("\%E %E\n",arcs[3][0],arcs[3][1]);

  double theta_min[int(arcs[0][0])];
  double theta_max[int(arcs[0][0])];
  double delta_theta[int(arcs[0][0])];

  for(int i=1;i<=arcs[0][0];i++){
    theta_min[i]  = arcs[i][0] - arcs[0][1];
    theta_max[i]  = arcs[i][1] + arcs[0][1];
    delta_theta[i]= (theta_max[i]-theta_min[i]);

  }

  for(int i=1;i<=arcs[0][0]+1;i++){


//printing arc 'i' in anti-clock direction
    theta = theta_min[i];

    for(int j=0;j<=npts/2;j++){
      if(arg_sqrt(Df0Dtheta_in, source_in, theta, pert_params, _r_e)>0.0){
        r_p = _r_e+dr_plus(f1_in, Df0Dtheta_in, kappa2, source_in, theta, pert_params,_r_e);
        //r_m = _r_e+dr_minus(f1_in, Df0Dtheta_in, kappa2, source_in, theta, pert_params,_r_e);
        if(file_in==NULL){
          printf("%E %E\n",r_p*cos(theta),r_p*sin(theta));
          //printf("%E %E\n",r_m*cos(theta),r_m*sin(theta));
        } else {
          fprintf(file_in,"%E %E\n",r_p*cos(theta),r_p*sin(theta));
          //fprintf(file_in,"%E %E\n",r_m*cos(theta),r_m*sin(theta));
        }
      }
    
      theta += delta_theta[i]/(npts/2);
    }


//printing arc 'i' in clock direction
    theta = theta_max[i];

    for(int j=0;j<=npts/2;j++){

      if(arg_sqrt(Df0Dtheta_in, source_in, theta, pert_params, _r_e)>0.0){
        //r_p = _r_e+dr_plus(f1_in, Df0Dtheta_in, kappa2, source_in, theta, pert_params,_r_e);
        r_m = _r_e+dr_minus(f1_in, Df0Dtheta_in, kappa2, source_in, theta, pert_params,_r_e);
        if(file_in==NULL){
          //printf("%E %E\n",r_p*cos(theta),r_p*sin(theta));
          printf("%E %E\n",r_m*cos(theta),r_m*sin(theta));
        } else {
          //fprintf(file_in,"%E %E\n",r_p*cos(theta),r_p*sin(theta));
          fprintf(file_in,"%E %E\n",r_m*cos(theta),r_m*sin(theta));
        }
      }
    
      theta -= delta_theta[i]/(npts/2);
    }

//separing the arcs
        if(file_in==NULL){
          printf("\n\n");
        } else {
          fprintf(file_in,"\n\n");
        }


  }


}

#endif
