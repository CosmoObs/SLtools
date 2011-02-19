/** @file
* This module compute the cross section for deformation arcs 
* - \f$ \bar{f}_n(\theta)\f$  and their derivatives
* - \f$ r_{\rm crit}(\theta)\f$  Radial coordinate of the Tangential Critical Curve
* All of these expressions are described in .../sltools/PerturbativeMethod/writeups/Report_on_Perturbative_Method.pdf
*/

/** @package perturbative_method */

#ifndef SIGMA_DEFARCS_H
#include "perturbative_method.h"
#include "general_methods.h"

double sigma_defarcs(f_type f1_in, f_type D2f0Dtheta2_in, double pert_params[], double kappa2, double _r_e, double r_th, int npts=50){
  double rth_fact=(fabs(r_th)*fabs(r_th)+1.)/(pow(fabs(r_th)*fabs(r_th)-1,2));
  double x[npts], w[npts];
  double Pi = 3.1415926535897932384626433832795;
  gauleg(0, 2.0*Pi, x, w, npts);

  double out = 0.0;
  double r_crit_tmp = 0.0;
  
  for(int i=0;i<npts;i++){
    r_crit_tmp = r_crit(f1_in, D2f0Dtheta2_in, kappa2, x[i], pert_params, _r_e);
    out+=w[i]*r_crit_tmp*r_crit_tmp;
  }
  return kappa2*kappa2*rth_fact*out;
}

double sigma_defarcs_sis(double _r_e, double r_th){
double Pi = 3.1415926535897932384626433832795;
double factor=(fabs(r_th*r_th)+1)/pow(fabs(r_th*r_th)-1 ,2);
return 2.*Pi*_r_e*_r_e*factor;
}
#endif 
