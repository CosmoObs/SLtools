/** @file
* This module is useful to compute the main fields of the perturbative approach, considering external shear + continous matter
*
* Basically, it computes \f$ f_0(\theta)\f$, \f$ f_1(\theta)\f$ and their derivatives.
*
* Its main parameters are ext_params[], where
*
* - s_ext_params[0]= continuous matter value, i.e, \f$ \kappa_{\rm ext}\f$ 
*
* - s_ext_params[1]= external shear value, i.e., \f$ \gamma_{\rm ext}\f$ 
*
* - s_ext_params[2]= orientation of the external shear, i.e., \f$ \theta_{\gamma}\f$ 
*
* - r_e_ : Einstein Radius of the unperturbed model, i.e, \f$ R_{\rm E}\f$
*/

/** @package external_shear
*
* This package compute:
*
* - quantities related with perturbed potential due an external shear + continous matter
*
*/
#ifndef EXTERNAL_SHEAR_H
#define EXTERNAL_SHEAR_H
/**  First perturbative field for the external shear + continous matter model.
*
*
*  \f$ f_{1}(\theta) \equiv R_{\rm E}\left[\kappa_{\rm ext}-\gamma_{\rm ext}\cos{2(\theta-\theta_\gamma)} \right]\f$
*
* \param  theta : angular coordinate
*
* \param  s_ext_params[] : parameters of the external shear + continuous matter.
*  \param _r_e : Einstein Radius of the unperturbed model
* \return \f$ f_1(\theta)\f$
**/
double f1_exsh(double theta,double s_ext_params[],double _r_e){
 double kex=s_ext_params[0], gex=s_ext_params[1],theta_g=s_ext_params[3];
 double bar_t=2*(theta-theta_g),f1ex;
 return f1ex=_r_e*(kex-gex*cos(bar_t));
}  

/**  Second perturbative field for the external shear + continous matter model.
*
*
*  \f$ \frac{d f_0(\theta)}{d\theta} \equiv R^2_{\rm E}\gamma_{\rm ext}\sin{2(\theta-\theta_\gamma)}\f$
*
* \param  theta : angular coordinate
*
* \param  s_ext_params[] : parameters of the external shear + continuous matter.
*  \param _r_e : Einstein Radius of the unperturbed model
* \return \f$ \frac{d f_0(\theta)}{d\theta}\f$
**/
double df0dtheta_exsh(double theta,double s_ext_params[],double _r_e){
   double kex=s_ext_params[0], gex=s_ext_params[1],theta_g=s_ext_params[3];
   double bar_t=2*(theta-theta_g), df0ex;
   return df0ex=pow(_r_e,2)*gex*sin(bar_t); 
}
/**  Third perturbative field for the external shear + continous matter model.
*
*
*  \f$ \frac{d^2 f_0(\theta)}{d\theta^2} \equiv 2\,R^2_{\rm E}\gamma_{\rm ext}\cos{2(\theta-\theta_\gamma)}\f$
*
* \param  theta : angular coordinate
*
* \param  s_ext_params[] : parameters of the external shear + continuous matter.
*  \param _r_e : Einstein Radius of the unperturbed model
* \return \f$ \frac{d^2 f_0(\theta)}{d\theta^2}\f$
**/
double df0dtheta2_exsh(double theta,double s_ext_params[],double _r_e){
   double kex=s_ext_params[0], gex=s_ext_params[1],theta_g=s_ext_params[3];
   double bar_t=2*(theta-theta_g),ddf0ex;
   return ddf0ex=2.*pow(_r_e,2)*gex*cos(bar_t); 
}

#endif