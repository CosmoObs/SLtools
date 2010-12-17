/** @file
* Function to calculate the substructure fields  (valid for any radial profile)
*
//   ***********************  [Have to fix all Doxygen stuff still -mg]
**/

// By MSSG, based on pseudo_elliptical_model.h HSDM code
// Start: Nov 25, 2010

# ifndef SUBSTRUCTURE_MODEL_H
# define SUBSTRUCTURE_MODEL_H

#include <cmath> 

typedef double (*f_type2)   (double r, double subst_params[]); 



double f1_sub(f_type2 alpha_in, double R_E, double theta, double subst_params[]){
  
  double rp = subst_params[0];
  double thetap = subst_params[1];
  double mp = subst_params[3];
  double pot_params[]={mp,subst_params[4]};
  
  double zetaRE = sqrt(R_E*R_E - 2*R_E*rp*cos(theta-thetap) + rp*rp);
  
  double f1_sub_tmp = alpha_in(zetaRE,pot_params)*(R_E-rp*(cos(theta-thetap))) / zetaRE;

  double pot_params[]= thetap, rp ;

return f1_sub_tmp;

}



double df0dtheta_sub(f_type2 alpha_in, double R_E, double theta, double subst_params[]){

  double  rp = subst_params[0];
  double thetap = subst_params[1];
  double mp = subst_params[3];
  double pot_params[]={mp,subst_params[4]};
  
  double zetaRE = sqrt(R_E*R_E - 2*R_E*rp*cos(theta-thetap) + rp*rp);
  
  double OmegaE =  R_E*rp*sin(theta-thetap); 

  double  df0dt_sub_tmp = alpha_in(zetaRE,pot_params)*OmegaE/ zetaRE;

  return df0dt_sub_tmp;

}


double d2f0dtheta2_sub(f_type2 alpha_in, double shear_in, double R_E, double theta, double subst_params[]){
  
  double  rp = subst_params[0];
  double thetap = subst_params[1];
  double mp = subst_params[3];
  double pot_params[]={mp,subst_params[4]};
  
  double zetaRE = sqrt(R_E*R_E - 2*R_E*rp*cos(theta-thetap) + rp*rp);
  
  double kappaE = mp / (2*zetae);

  double OmegaE =  R_E*rp*sin(theta-thetap); 
  double gammaE =  shear_in(zetae,pot_params); 

  double  d2f0dt2_sub_tmp = (alpha_in(zetaRE,pot_params)*(R_E*rp*(cos(theta-thetap))) / zetaRE) - 2*gammaE*(OmegaE*OmegaE/(zetaRE*zetaRE));


  return d2f0dt2_sub_tmp;

}

}
