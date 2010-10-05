/** @file
* Example of doxygen documentation for C functions FIXME. 
*/

/** @package nfw_model
*  Package to compute quantities related to NFW models.
*
*  Detailed descrition FIXME
*
*/
#ifndef EXTERNAL_SHEAR_H
#define EXTERNAL_SHEAR_H

double f1_exsh(double theta,double s_ext_params[]){
 double kex=s_ext_params[0], gex=s_ext_params[1],theta_g=s_ext_params[3];
 double bar_t=2*(theta-theta_g),f1ex;
 return f1ex=_r_e*(kex-gex*cos(bar_t));
}  

double df0dtheta_exsh(double theta,double s_ext_params[]){
   double kex=s_ext_params[0], gex=s_ext_params[1],theta_g=s_ext_params[3];
   double bar_t=2*(theta-theta_g), df0ex;
   return df0ex=pow(_r_e,2)*gex*sin(bar_t); 
}

double df0dtheta2_exsh(double theta,double s_ext_params[]){
   double kex=s_ext_params[0], gex=s_ext_params[1],theta_g=s_ext_params[3];
   double bar_t=2*(theta-theta_g),ddf0ex;
   return ddf0ex=2.*pow(_r_e,2)*gex*cos(bar_t); 
}

#endif