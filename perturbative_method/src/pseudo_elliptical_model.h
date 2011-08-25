/** @file
* Function to calculate the main perturbative fields related with Pseudo-Elliptical Models (valid for any radial profile)
*
// Note: (pert_params[0],pert_params[1], pert_params[2]....) must be (ellipticity, flag control, mass and so on..)  
//
// The three kind of parameterization accepted in this functions are:
//
// For pert_params[1]=1, Is used \f$ a_{1\eta}=1-\eta, a_{2\eta}=1+\eta \f$, i.e Angle Deflection Model
// 
// For pert_params[1]=2, Is used \f$ a_{1\eta}=1-\eta, a_{2\eta}=1/(1-\eta) \f$, i.e Standard Parameterization
//
// For pert_params[1]=3, Is used \f$ a_{1\eta}=1, a_{2\eta}=1/(1-\eta)^2 \f$, i.e Keeton's Parameterization
//
// example: for the PNFW model, useful for comparison with gravlens, we must define pert_params[]={\f$ \eta \f$ ,3, \f$ \kappa_s \f$, \f$ r_s \f$}
*/

/** @package pseudo_elliptical_models
* Package to compute quantities related to Pseudo-Elliptical Models in the framework of Perturbative Approach
*
* Detailed description we find in .../sltools/PerturbativeMethod/writeups/Report_on_Perturbative_Method.pdf
*/
# ifndef PE_MODEL_H
# define PE_MODEL_H
#include <math.h> 

typedef double (*f_type2)   (double r, double pert_params[]); 

double f1_pe(f_type2 alpha_in, double R_E, double theta, double pert_params[]);

double df0dtheta_pe(f_type2 alpha_in,double R_E, double theta, double pert_params[]);

double d2f0dtheta2_pe(f_type2 alpha_in,f_type2 shear_in, double R_E, double theta, double pert_params[]);

#endif
