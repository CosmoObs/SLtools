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
#define TAM_MAX 1000

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

    float x_cord[TAM_MAX],y_cord[TAM_MAX];

    FILE *encurv;
    if(cflag==1){encurv=fopen("cc.dat","r");}
    else{encurv=fopen("ca.dat","r");}
    int i=0;
    if (encurv== NULL) 
    {
      printf("Nao foi possivel abrir o arquivo.\n");
//        break;
//        exit(1);
     }else{
      while (! feof(encurv))
      {
        i++;
        fscanf(encurv," %f %f",&x_cord[i], &y_cord[i]);
	
      }}
      int n_points = i++-1;
      printf("The number of points is %i\n",n_points);


  double d_nor=sqrt(pow((x_cord[1]-x_cord[n_points]),2)+pow((y_cord[1]-y_cord[n_points]),2) );
  double phi[TAM_MAX],r_cord[TAM_MAX],wphi[TAM_MAX];  
  double r_pertm[TAM_MAX],theta[TAM_MAX];

  for(int i=1;i<=n_points;i++){
    phi[i]=atan(y_cord[i]/x_cord[i]);
    r_cord[i]=sqrt(pow(x_cord[i],2)+pow(y_cord[i],2));
  }

  wphi[1]=fabs(phi[1]-phi[2]);

  for (int j=2;j<=n_points; j++){
    int j1=j-1;
    wphi[j]=fabs(phi[j]-phi[j1]);
  }

  double sum1=0.0, sum_nor=0.0;
    for(int k=1; k<=n_points; k++){
      theta[k]=phi[k]; 
      if(cflag==1){ r_pertm[k]=r_crit(f1_in, D2f0Dtheta2_in, kappa2, theta[k], pert_params,_r_e);}
      else{ r_pertm[k]=r_caust(Df0Dtheta_in,D2f0Dtheta2_in,theta[k], pert_params, _r_e);}
      sum1+= wphi[k]*pow((r_cord[k]-r_pertm[k]),2);
     sum_nor+=wphi[k]*pow(r_cord[k],2);
    }
//      double fom;
//     if(cflag==1){fom=sum1/sum_nor;}else{fom=sum1;}
     double fom=sum1/sum_nor;
//    double fom=sum1/d_nor;  
    return fom;
}

#endif
