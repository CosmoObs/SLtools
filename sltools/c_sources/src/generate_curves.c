/** @file
* Module to generate curves useful for the perturbative approach, as
* 
*  Tangential Critical Curve and  Tangential Caustic. This function need
* - main field of the perturbative approach, i.e, \f$ f_n(\theta)\f$ and their derivatives
* - pert_params[] :a vector that contains all the perturbation parameters
* - npts: number of theta divisions 
* - Einstein Radius: \f$ R_{_{\mathrm{E}}}\f$ 
* - two files to write the tangential curves 
*
*
*  Elliptical source, contour plot. It need
*  - structure corresponding to the source: size, ellipticity, orientation,etc
*  - npts: number of theta divisions
*  - a file to write the plot in a text file
*
*
*  Curves of Constant Distortion. This function need
* - main field of the perturbative approach, i.e, \f$ f_n(\theta)\f$ and their derivatives
* - pert_params[] :a vector that contains all the perturbation parameters
* - npts: number of theta divisions 
* - Einstein Radius: \f$ R_{_{\mathrm{E}}}\f$ 
* - A POSITIVE VALUE of the length-to-width ratio threshold, i.e \f$ R_{\rm th} \f$
* - two files to write theses curves in a text file (with the following format x1, x2, y1, y2), each for both value of \f$ R_{\rm th}\f$
*
* This function computes automatically both positive and negative constant distortion curves
*
*/

/** @package perturbative_method
* 
* To plot the critical and caustic curves (model independent)
* 
* To plot of an elliptical source, not aligned to the main axis, in its corresponding plane 
*/

#include "generate_curves.h"
#include "config.h"

/** A functions to print the curves model (tangential critical curve,tangential caustic curves)
*
*  \param Df0Dtheta_in function related to the perturbed potential
*  \param f1_in function related to the perturbed potential
*  \param D2f0Dtheta2_in function related to the perturbed potential
*  \param pert_params[] a vector that contains all the perturbation parameters
*  \param kappa2 \f$  \kappa_2 = 1 - \left[\frac{d^2 \phi_0(r)}{dr^2}\right]_{r=r_e} \f$
*  \param npts number of theta divisions 
*  \param _r_e Einstein Radius (\f$ R_{_{\mathrm{E}}}\f$ )
*  \param file_in1 a file to write tangential critical curve
*  \param file_in2 a file to write tangential caustic
*  \return curve datas in the screen or the curve data in a text file
*
*  \sa f_type, r_crit, caustic_y1, caustic_y2, y1_src,y2_src
*/

void plot_curves(f_type f1_in, f_type Df0Dtheta_in, f_type D2f0Dtheta2_in, double pert_params[], double kappa2, double _r_e, int npts, FILE *file_in1){
  
  double theta = 0.0;
  double r_tcrit=0.0;
  double y1_caust=0.0;
  double y2_caust=0.0;

  for(int i=0;i<=npts;i++){
  r_tcrit= r_crit(f1_in,D2f0Dtheta2_in,kappa2,theta,pert_params,_r_e); 
  y1_caust=caustic_y1(Df0Dtheta_in,D2f0Dtheta2_in,theta,pert_params,_r_e);
  y2_caust=caustic_y2(Df0Dtheta_in,D2f0Dtheta2_in,theta,pert_params,_r_e);

  if(file_in1==NULL){
        printf("%f %f %f %f \n",r_tcrit*cos(theta),r_tcrit*sin(theta),y1_caust,y2_caust);
      } else {
        fprintf(file_in1,"%f %f %f %f\n",r_tcrit*cos(theta),r_tcrit*sin(theta),y1_caust,y2_caust);
      }

   theta+= 6.283185308/npts;
}

}

/** A functions to print the source contour (the plot in the source plane...)
*
*  \param source_in a struct that define elliptical sources
*  \param npts number of theta divisions (If itsn't given, the code assumes npts=1000)
*  \param file_in a file to write the source contour
*  \return source curve data in the screen or the source curve data in a text file
*/


void plot_sources(elliptical_source source_in, int npts, FILE *file_in){
  double y1s_plot=0.0;
  double y2s_plot=0.0;
  double theta=0.0;
  
  for(int i=0;i<=npts;i++){
    
    y1s_plot=y1_src(source_in, theta);
    y2s_plot=y2_src(source_in, theta);
    
    if(file_in==NULL){
        printf("%E %E\n",y1s_plot,y2s_plot);
      } else {
        fprintf(file_in,"%f %f\n",y1s_plot,y2s_plot);
      }
  
    theta+= 6.283185308/npts;
}
}

/** A functions to print the constant distortion curves in both planes
*
*  \param f1_in function related to the perturbed potential
*  \param Df0Dtheta_in function related to the perturbed potential
*  \param D2f0Dtheta2_in function related to the perturbed potential
*  \param pert_params[] a vector that contains all the perturbation parameters
*  \param kappa2 \f$  \kappa_2 = 1 - \left[\frac{d^2 \phi_0(r)}{dr^2}\right]_{r=r_e} \f$
*  \param npts number of theta divisions 
*  \param _r_e Einstein Radius (\f$ R_{_{\mathrm{E}}}\f$ )
*  \param file_in1 a file to write the constant distortion curve for \f$ R_{\rm th}>0 \f$ 
*  \param file_in2 a file to write the constant distortion curve for \f$ R_{\rm th}<0 \f$ 
*  \return curve datas in the screen or the curve data in a text file
*
*  \sa f_type, x_th, y1_th,y2_th
*/

void plot_rth_curves(f_type f1_in, f_type Df0Dtheta_in, f_type D2f0Dtheta2_in, double pert_params[], double kappa2, double _r_e, double r_th, int npts, FILE *file_in1, FILE *file_in2){
  
  double theta = 0.0;
  double r_thp=0.0;
  double y1_rthp=0.0;
  double y2_rthp=0.0;
  double r_thn=0.0;
  double y1_rthn=0.0;
  double y2_rthn=0.0;
  
  for(int i=0;i<=npts;i++){
    
  r_thp= x_th(f1_in,D2f0Dtheta2_in, kappa2, theta,pert_params,r_th,_r_e)+_r_e; 
  y1_rthp=y1_th(f1_in,Df0Dtheta_in,D2f0Dtheta2_in,theta,kappa2,pert_params,_r_e,r_th);
  y2_rthp=y2_th(f1_in,Df0Dtheta_in,D2f0Dtheta2_in,theta,kappa2,pert_params,_r_e,r_th);

  r_thn= x_th(f1_in,D2f0Dtheta2_in, kappa2, theta,pert_params,-r_th,_r_e)+_r_e; 
  y1_rthn=y1_th(f1_in,Df0Dtheta_in,D2f0Dtheta2_in,theta,kappa2,pert_params,_r_e,-r_th);
  y2_rthn=y2_th(f1_in,Df0Dtheta_in,D2f0Dtheta2_in,theta,kappa2,pert_params,_r_e,-r_th);
  
  if(file_in1==NULL){
        printf("%f %f %f %f\n",r_thp*cos(theta),r_thp*sin(theta), y1_rthp,y2_rthp);
      } else {
        fprintf(file_in1,"%f %f %f %f \n",r_thp*cos(theta),r_thp*sin(theta), y1_rthp,y2_rthp);
      }

  if(file_in2==NULL){
        printf("%f %f %f %f\n",r_thn*cos(theta),r_thn*sin(theta), y1_rthn,y2_rthn);
      } else {
        fprintf(file_in2,"%f %f %f %f \n",r_thn*cos(theta),r_thn*sin(theta), y1_rthn,y2_rthn);
      }

   theta+= 6.283185308/npts;
}

}
