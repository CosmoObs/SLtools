/** @file
*  Package to compute the goodness of fit, to compare the solution of the perturbative method with fully numerical solution of the tangential  critical and caustic curves FIXME
*/

/** @package goodness_of_fit
*  Package to compute the goodness of fit, to compare the solution of the perturbative method with fully numerical solution of the tangential  critical and caustic curves
*
*  Detailed descrition FIXME
*
*/


#ifndef GOF_PER_METHOD_H
#define GOF_PER_METHOD_H


#include <cmath> 
#include "perturbative_method.h"
#define TAM_MAX 10000

/** A function to compute the goodness of fit to the interesting curves (tangential critical or tangential caustic curve) of any model, aiming to compare the perturbative method solution to the fully numerical solution
*
* Knowing the data set defining the points of the curves (both tangential critical or tangential caustic curve), we have
*
* the radial coordinate of the point (\f$ x_{\rm cord}, y_{\rm cord}\f$) \f$ r_{\rm curv}=\sqrt{x_{\rm cord}^2+y_{\rm cord}^2} \f$
*
* the angular coordinate of the (\f$ x_{\rm cord}, y_{\rm cord}\f$) \f$ \phi_=\atan{y_{\rm cord}/y_{\rm cord}} \f$
*
* the weight due the non uniform distribution of points \f$ w(i)=\phi(i)-\phi(i-1)\f$
*
* the radial coordinate of the curve, obtaining with the perturbative method \f$ r_{\rm pert} \f$ calculated for each angular coordinate \f$\phi(i) \f$
*  \param f1_in function related to the perturbed potential
*  \param Df0Dtheta_in function related to the perturbed potential
*   \param D2f0Dtheta2_in function related to the perturbed potential (second derivative)
*  \param kappa2 \f$  \kappa_2 = 1 - \left[\frac{d^2 \phi_0(r)}{dr^2}\right]_{r=r_e} \f$
*  \param pert_params[] a vector that contains all the perturbation parameters
*  \param _r_e Einstein Radius (\f$ R_{_{\mathrm{E}}}\f$ )
*   \param x_cord[] vector containing the coordinate x of the data set
*   \param y_cord[] vector containing the coordinate y of the data set
*  \param cflag :flag to decide which of the comparison will be do (If cflag=1 comparison with the tangential critical curve, If cflag=2 comparison with the tangential caustic curve)
*  \param n_points (number of points of the data set)
*  \return \f$ gof:=\dfrac{\sum_{i=1}^N w(i)\left[r_{\rm curv}^2(i)-r_{\rm pert}^2(i)\right]}{\sum_{i=1}^N w(i)r_{\rm curv}^2(i)} \f$
*  \sa r_crit, r_caust
*
*/


double gof_per_method(f_type f1_in, f_type Df0Dtheta_in,f_type D2f0Dtheta2_in,double kappa2,double pert_params[],double _r_e, int cflag){

    float xct_cord[TAM_MAX],yct_cord[TAM_MAX];
    float xca_cord[TAM_MAX],yca_cord[TAM_MAX];
    
    FILE *incurv1;
    incurv1=fopen("cc.dat","r");
    int ict=0;
    if (incurv1== NULL) 
    {
      printf("Nao foi possivel abrir o arquivo.\n");
      exit(1);
     }else{
      while (! feof(incurv1))
      {
        ict++;
        fscanf(incurv1," %f %f",&xct_cord[ict], &yct_cord[ict]);	
      }}
      int n_crit = ict++-1;
      printf("The number of points is %i\n",n_crit);
    

    FILE *incurv2;
    incurv2=fopen("ca.dat","r");
    int ica=0;

    if (incurv2== NULL) 
    {
      printf("Nao foi possivel abrir o arquivo.\n");
      exit(1);
     }else{
      while (! feof(incurv2))
      {
        ica++;
        fscanf(incurv2," %f %f",&xca_cord[ica], &yca_cord[ica]);
	
      }}
      int n_ca = ica++-1;
      printf("The number of points is %i\n",n_ca);

     int n_points=n_crit;

  fclose(incurv1);
  fclose(incurv2);
  double phi[TAM_MAX],phi_ca[TAM_MAX];
  double r_ct[TAM_MAX],r_ca[TAM_MAX],phi_ca_pm[TAM_MAX];

  FILE *angles;
  angles=fopen("ca_angles.dat","w");

  for(int i=1;i<=n_points;i++){
    phi[i]=atan(yct_cord[i]/xct_cord[i]);
    phi_ca[i]=atan(yca_cord[i]/xca_cord[i]);
    phi_ca_pm[i]=atan(fabs(caustic_y2(Df0Dtheta_in,D2f0Dtheta2_in,phi[i],pert_params,_r_e))/fabs(caustic_y1(Df0Dtheta_in,D2f0Dtheta2_in,phi[i],pert_params,_r_e)));
    r_ct[i]=sqrt(pow(xct_cord[i],2)+pow(yct_cord[i],2));
    r_ca[i]=sqrt(pow(xca_cord[i],2)+pow(yca_cord[i],2));
    fprintf(angles,"%f %f\n",phi_ca[i], phi_ca_pm[i]);
  }
  fclose(angles);
//   wphi[1]=fabs(phi[1]-phi[2]);

  double wphi[TAM_MAX],wphi_ca[TAM_MAX];
  double dx_ct[TAM_MAX],dx_ca[TAM_MAX];
  double dy_ct[TAM_MAX],dy_ca[TAM_MAX];

  for (int j=1;j<n_points; j++){
    int j1=j+1;
    wphi[j]=fabs(phi[j1]-phi[j]);
    wphi_ca[j]=fabs(phi_ca[j1]-phi_ca[j]);
    dx_ct[j]=fabs(xct_cord[j1]-xct_cord[j]);
    dx_ca[j]=fabs(xca_cord[j1]-xca_cord[j]);
    dy_ct[j]=fabs(yct_cord[j1]-yct_cord[j]);
    dy_ca[j]=fabs(yca_cord[j1]-yca_cord[j]);
  }
  wphi[n_points]=fabs(phi[n_points]-phi[n_points-1]);
  wphi_ca[n_points]=fabs(phi_ca[n_points]-phi_ca[n_points-1]);

//   double area_ct=0.0, area_ca=0.0;
// 
//   for(int j2=1;j2<n_points;j2++){
//     area_ct+=dx_ct[j2]*dy_ct[j2];
//     area_ca+=dx_ca[j2]*dy_ca[j2];
//   }
// 
//   double cnor_ct=area_ct;
//   double cnor_ca=area_ca;
  
  double sum_ct, sum_nor_ct, sum_ca, sum_nor_ca;
  sum_ct=0.0, sum_nor_ct=0.0;
  sum_ca=0.0, sum_nor_ca=0.0;
  double theta[TAM_MAX],r_pm_ct[TAM_MAX],r_pm_ca[TAM_MAX];
  
    for(int k=1; k<=n_points; k++){
      theta[k]=phi[k]; 
      r_pm_ct[k]=r_crit(f1_in, D2f0Dtheta2_in, kappa2, theta[k], pert_params,_r_e);
      sum_ct+= wphi[k]*pow((r_ct[k]-r_pm_ct[k]),2);
      sum_nor_ct+= wphi[k]*pow(r_ct[k],2);
//
      r_pm_ca[k]=r_caust(Df0Dtheta_in,D2f0Dtheta2_in,theta[k], pert_params, _r_e);
      sum_ca+= wphi_ca[k]*pow((r_ca[k]-r_pm_ca[k]),2);
      sum_nor_ca+=wphi_ca[k]*pow(r_ca[k],2);
    }
    double fom_crit=(sum_ct/sum_nor_ct);
    double fom_ca=(sum_ca/sum_nor_ca);
    
//     double fom_crit=sum_ct/cnor_ct;
//     double fom_ca=sum_ca/cnor_ca;

    if(pert_params[0]<=1E-7){
        fom_crit=0.0;
        fom_ca=0.0;
    }

    double fom;
    if(cflag==1){
        fom=fom_crit;
    }else{
        fom=fom_ca;
    }
    
    return fom;
}

#endif