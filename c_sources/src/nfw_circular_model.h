/** @file
* Package to calcule quantities related to the NFW model
*
* - Lensing Function (Angle Deflection, Convergence, Shear)
*
* - Einstein Radius (Use bracketing_lambda_t and gsl functions)
*
* - Function useful to the Perturbative Approach \f$ \kappa_2 \f$
*
* - Function related with the third derivative of the lensing potential \f$ c_3 \f$ 
* 
*  This module need 
* - pot_params[]={\f$\kappa_s, r_s \f$}
*/

/** @package circular_models
*  Package to compute quantities related to the circular NFW model.
*
*
*/
#ifndef NFW_CIRCULAR_MODEL_H
#define NFW_CIRCULAR_MODEL_H

double F(double X);


double Flinha(double X);


double conv_nfw_circ(double r, double pot_params[]);


double conv_nfw_circ_prime(double r, double pot_params[]);


double alpha_nfw_circ(double r, double pot_params[]);


double shear_nfw_circ(double r, double pot_params[]);


double kappa2_nfw(double pot_params[],double _r_e_nfw);


double c3_nfw(double pot_params[],double _r_e_nfw);


double lambda_t_nfw_circ(double r, double pot_params[]);


double re_find_func_nfw(double r, void *p);


double re_find_func_nfw(double r, double params[]);


double bracketing_lambda_t(double f(double r, double params[]), 
                double params[], double out[], double z, double dz);

double r_e_nfw_find(double pot_params[], double *est_err_out, double x_lo,
                 double x_hi, int max_iter, double relative_error, int v);



#endif
