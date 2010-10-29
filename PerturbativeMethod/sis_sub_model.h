/** @file
* Example of doxygen documentation for C functions FIXME. 
*/

/** @package sis_sub_model
*  SIS as a substructure
*
*  Detailed descrition FIXME
*
*/



#ifndef SIS_SUB_MODEL_H
#define SIS_SUB_MODEL_H

double _r_e_model = 1.0;


//pot_params[0] = mp
//pot_params[1] = rp
//pot_params[2] = thetap
/**  \f$ f_1(\theta) \f$ for SIS as a substructure
*
*  \param theta angular coordinate in the lens plane
*  \param pot_params[0] = mp
*  \param pot_params[1] = rp
*  \param pot_params[2] = thetap
*  \return \f$ f_1(\theta) \f$

*  \sa Df0Dtheta_pert_SIS, f_type
*/
double f1_pert_SIS(double theta, double pot_params[]){	
  double mp=pot_params[0], rp=pot_params[1],theta_p=pot_params[2];
  double re=_r_e_model, bar_theta=theta-theta_p;
  double tmp1 = mp*(re-rp*cos(bar_theta));
  double tmp2 = sqrt(pow(re,2)-2.0*rp*re*cos(bar_theta)+pow(rp,2));

  return tmp1/tmp2 ;
}

/**  \f$ \frac{d f_0(\theta)}{d \theta} \f$ for SIS as a substructure
*
*  \param theta angular coordinate in the lens plane
*  \param pot_params[0] = mp
*  \param pot_params[1] = rp
*  \param pot_params[2] = thetap
*  \return \f$ \frac{d f_0(\theta)}{d \theta} \f$

*  \sa f1_pert_SIS, f_type
*/
double Df0Dtheta_pert_SIS(double theta, double pot_params[]){
    double mp=pot_params[0], rp=pot_params[1],theta_p=pot_params[2];
  double re=_r_e_model, bar_theta=theta-theta_p;
  double tmp1 = mp*re*rp*sin(bar_theta);
  double tmp2 = sqrt(pow(re,2)-2.0*rp*re*cos(bar_theta)+pow(rp,2));
  return tmp1/tmp2;
}
/**********************************************************************************************************************/

double D2f0Dtheta2_pert_SIS(double theta, double pot_params[]){
double mp=pot_params[0],rp=pot_params[1],thetap=pot_params[2];
double bar_theta=theta-thetap, re=_r_e_model;
double Rsq=pow(re,2)-2.0*rp*re*cos(bar_theta)+pow(rp,2);
return (mp*re*rp/sqrt(Rsq))*(cos(bar_theta)-re*rp*pow(sin(bar_theta),2)/Rsq);
}


double kappa_2_SIS(double pot_params[]){
return 1.0;
}
#endif
