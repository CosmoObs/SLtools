/** @file
*  Functions that are useful to compare the perturbative method with the exact solution. 
*
*  These functions are:
*
*  \f$ \mathcal{D}^2\f$ for critical curves \sa d2_crit
*
*  \f$ \mathcal{D}^2\f$ for caustics \sa d2_caust
*
*  \f$ D_{\rm max}:=3|c_3|dr_{\rm max}^2\f$ \sa lin_grad
*/

/** @package comparison_with_exact_solution
*
*  This Package compare the solution of the perturbative method with exact solution, through the calculation of the mean weighted squared fractional difference \f$ \mathcal{D}^2 \f$ and by the condition of the linearity of the gradient of the lensing potential \f$ D_{\rm max} \f$.
*
* It needs to know:
*
* The main fields of the Perturbative Approach, i.e., \f$ f_0(\theta)\f$, \f$ f_1(\theta) \f$ and their derivatives.
*
* The Einstein radius for a model \f$R_{\rm E}\f$, the second derivative of the potential \f$ \kappa_2\f$ and the third derivative of the potential \f$c_{3}\f$
*
* The other quantities as critical radius \f$ r_{\rm crit}\f$ and the radial coordinates of the caustic \f$r_{\rm caust}\f$
*
* A data file containing the exact solution for tangential caustic and critical curve of a model.
*
*
*/

#ifndef COMPARSION_FNS_H
#define COMPARSION_FNS_H

#include <cmath> 
#include "perturbative_method.h"
#define TAM_MAX 10000

/** A function to compute the mwsrdf for tangential critical curve, aiming to compare  the perturbative method solution with the fully numerical solution
*
* Knowing the data set defining the points of the tangential critical curve , we have
*
* the radial coordinate of the point (\f$ x_{1,\rm crit}, x_{2,\rm crit}\f$) \f$ R_{\rm fns}=\sqrt{x_{1,\rm crit}^2+x_{2,\rm crit}^2} \f$
*
* the weight due the non uniform distribution of points \f$ w(i)=\theta_L(i+1)-\theta_L(i)\f$
*
* the radial coordinate of the tangential critical curve, obtaining with the perturbative method \f$ R_{\rm crit} \f$ calculated for each angular coordinate \f$ \theta_L(i) \f$
*
*  \param f1_in function related to the perturbed potential
*   \param D2f0Dtheta2_in function related to the perturbed potential (second derivative)
*  \param kappa2 \f$  \kappa_2 = 1 - \left[\frac{d^2 \phi_0(r)}{dr^2}\right]_{r=r_e} \f$
*  \param pert_params[] a vector that contains all the perturbation parameters
*  \param _r_e Einstein Radius (\f$ R_{_{\mathrm{E}}}\f$ )
*   \param xct_cord[] vector containing the coordinate \f$ x_1 \f$ of the data set
*   \param yct_cord[] vector containing the coordinate \f$ x_2 \f$ of the data set
*  \param n_points (number of points of the data set)
*  \return \f$ \mathcal{D}^2:=\sum_{i=1}^N w(i)\left[R_{\rm fns}^2(i)-R_{\rm crit}^2(i)\right]/\sum_{i=1}^N w(i)R_{\rm fns}^2(i) \f$
*  \sa r_crit
*
*/

double d2_crit(f_type f1_in, f_type D2f0Dtheta2_in,double kappa2, double pert_params[],double _r_e){

    float xct_cord[TAM_MAX],yct_cord[TAM_MAX];
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
      int n_points=n_crit;

      fclose(incurv1);

      double phi[TAM_MAX];
      double r_ct[TAM_MAX];

      for(int i=1;i<=n_points;i++){
        phi[i]=atan(yct_cord[i]/xct_cord[i]);
        r_ct[i]=sqrt(pow(xct_cord[i],2)+pow(yct_cord[i],2));
      }
      double wphi[TAM_MAX];
      int j1=0;
      for (int j=1;j<n_points; j++){
        j1=j+1;
        wphi[j]=fabs(phi[j1]-phi[j]);
      }
      wphi[n_points]=fabs(phi[n_points]-phi[n_points-1]);
      double sum_ct, sum_nor_ct;
      sum_ct=0.0, sum_nor_ct=0.0;
      double theta[TAM_MAX],r_pm_ct[TAM_MAX];

      for(int k=1; k<=n_points; k++){
          theta[k]=phi[k]; 
          r_pm_ct[k]=r_crit(f1_in, D2f0Dtheta2_in, kappa2, theta[k], pert_params,_r_e);
          sum_ct+= wphi[k]*pow((r_ct[k]-r_pm_ct[k]),2);
          sum_nor_ct+= wphi[k]*pow(r_ct[k],2);
      }
      
      double fom_crit=(sum_ct/sum_nor_ct);

      return fom_crit;
}

/** A function to compute the mwsrdf for tangential caustic curve, aiming to compare  the perturbative method solution with the fully numerical solution
*
* Knowing the data set defining the points of the tangential critical and caustic curves , we have
*
* the radial coordinate of the point (\f$ y_{1,\rm caust}, y_{2,\rm caust}\f$) \f$ r_{caust,\rm fns}=\sqrt{y_{1,\rm crit}^2+y_{2,\rm crit}^2} \f$
*
* the angular coordinate of the points (\f$ x_{1,\rm crit}, x_{2,\rm crit}\f$) \f$ \theta_L=\arctan{x_{2,\rm crit}/x_{1,\rm crit}} \f$
*
* the angular coordinate (\f$ y_{1,\rm caust}, y_{2,\rm caust}\f$) \f$ \theta_s=\arctan{y_{2,\rm caust}/y_{1,\rm caustt}} \f$
*
* the weight due the non uniform distribution of points \f$ w_s(i)=\theta_s(i+1)-\theta_s(i)\f$
*
* the radial coordinate of the tangential caustic curve, obtaining with the perturbative method \f$ R_{\rm caust} \f$ calculated for each angular coordinate \f$ \theta_L(i) \f$, such that \f$ \theta_s(\theta_L(i)) \f$ obtained from the perturbative approach match with the angle obtained with the fully numerical solution.
*
*  \param f1_in function related to the perturbed potential
*  \param Df0Dtheta_in function related to the perturbed potential
*   \param D2f0Dtheta2_in function related to the perturbed potential (second derivative)
*  \param kappa2 \f$  \kappa_2 = 1 - \left[\frac{d^2 \phi_0(r)}{dr^2}\right]_{r=r_e} \f$
*  \param pert_params[] a vector that contains all the perturbation parameters
*  \param _r_e Einstein Radius (\f$ R_{_{\mathrm{E}}}\f$ )
*   \param xct_cord[] vector containing the coordinate \f$ x_1 \f$ of the data set
*   \param yct_cord[] vector containing the coordinate \f$ x_2 \f$ of the data set
*  \param n_points (number of points of the data set)
*  \return \f$ \mathcal{D}^2:=\sum_{i=1}^N w(i)\left[R_{\rm fns}^2(i)-R_{\rm caust}^2(i)\right]/\sum_{i=1}^N w(i)R_{\rm fns}^2(i) \f$
*  \sa r_crit, r_caust
*
*/

double d2_caust(f_type f1_in, f_type Df0Dtheta_in, f_type D2f0Dtheta2_in, double kappa2, double pert_params[],double _r_e){

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
    double r_ca[TAM_MAX];  
    FILE *angles;
    angles=fopen("ca_angles.dat","w");
    for(int i=1;i<=n_points;i++){
    phi[i]=atan(yct_cord[i]/xct_cord[i]);
    phi_ca[i]=atan(yca_cord[i]/xca_cord[i]);
    r_ca[i]=sqrt(pow(xca_cord[i],2)+pow(yca_cord[i],2));
    }

//  computing the statistical weight
    double wphi_ca[TAM_MAX];
    double dphi_u[TAM_MAX],dphi_l[TAM_MAX];
    int i2=0;

    for(int i1=1;i1<n_points; i1++){
      i2=i1+1;
      wphi_ca[i1]=fabs(phi_ca[i2]-phi_ca[i1]);
     dphi_l[i1]=phi_ca[i1]-phi_ca[i1]/100.;
      dphi_u[i1]=phi_ca[i1]+phi_ca[i1]/100.;
    }
    wphi_ca[n_points]=fabs(phi_ca[n_points]-phi_ca[n_points-1]);
    dphi_u[n_points]=phi_ca[n_points]+phi_ca[n_points]/200.;
    
    // computng the matching angle in the source plane, to propose n_per=5*n_crit  
    int n_per=5*n_points;
    double y1_ca,y2_ca;
//     int icount=0;
    double theta_ca,theta_l=0.0;
    double theta_pm[TAM_MAX],phi_lens[TAM_MAX];
    double half_pi= 6.283185308/4.0;
    for(int i_per=1; i_per<=n_per;i_per++){
      y1_ca=caustic_y1(Df0Dtheta_in,D2f0Dtheta2_in,theta_l, pert_params, _r_e);
      y2_ca=caustic_y2(Df0Dtheta_in,D2f0Dtheta2_in,theta_l, pert_params, _r_e);
      theta_ca=atan(fabs(y2_ca)/fabs(y1_ca));
      for(int i_p2=1; i_p2<=n_points;i_p2++){
        if( (theta_ca>dphi_l[i_p2]) && (theta_ca<=dphi_u[i_p2])){
            theta_pm[i_p2]=theta_ca;
            phi_lens[i_p2]=theta_l;
        }
      }
      theta_l+=half_pi/n_per;
      }

    for(int k1=1; k1<=n_points;k1++){
    fprintf(angles,"%f %f %f %f \n",phi_ca[k1],theta_pm[k1],phi[k1],phi_lens[k1]);
    }
    
    fclose(angles);


    double sum_ca=0.0, sum_nor_ca=0.0;
    double r_pm_ca[TAM_MAX];

    for(int k=1; k<=n_points; k++){
      r_pm_ca[k]=r_caust(Df0Dtheta_in,D2f0Dtheta2_in,phi_lens[k], pert_params, _r_e);
      sum_ca+= wphi_ca[k]*pow((r_ca[k]-r_pm_ca[k]),2);
      sum_nor_ca+=wphi_ca[k]*pow(r_ca[k],2);
    }
    double fom_ca=(sum_ca/sum_nor_ca);

    if(pert_params[0]<=1E-7){fom_ca=0.0;}

    return fom_ca;
}
/** A function to compute the second term of the Taylor's Series along the tangential critical curve, considering its maximum value.
*
* First, it need the third derivative of the unperturbed potential, i.e.
*
* \f$ C_3(r)= \frac{4}{r_s(X^2_{\mathrm{E}}-1)}\left[\frac{1}{X_{\mathrm{E}}}-\frac{3}{2}\frac{X_{\mathrm{E}}\kappa(R_{\mathrm{E}})}{2\kappa_s}  \right] +2 \frac{\gamma(R_{\mathrm{E}})}{R_{\mathrm{E}}} \f$ 
*
*
* the radial coordinate of the tangential caustic curve, obtaining with the perturbative method \f$ R_{\rm caust} \f$ calculated for each angular coordinate \f$ \theta_L(i) \f$, such that \f$ \theta_s(\theta_L(i)) \f$ obtained from the perturbative approach match with the angle obtained with the fully numerical solution.
*
*  \param f1_in : function related to the perturbed potential
*   \param D2f0Dtheta2_in : function related to the perturbed potential (second derivative)
*   \param c3_model : The third derivative of the unperturbed potential
*  \param kappa2 \f$ :  \kappa_2 = 1 - \left[\frac{d^2 \phi_0(r)}{dr^2}\right]_{r=r_e} \f$
*  \param pert_params[] a vector that contains all the perturbation parameters
*  \param _r_e : Einstein Radius (\f$ R_{_{\mathrm{E}}}\f$ )
*  \param n_points (number of points of the data set)
*  \return \f$ D_{\rm max}:=3|c_3|dr_{\rm max}^2\f$
*  \sa r_crit
*
*/

double lin_grad(f_type f1_in,f_type D2f0Dtheta2_in, double c3_model, double kappa2_model, double pert_params[],double _r_e,int npts=100){

    double half_pi= 6.283185308/4.0;
    double theta=0.0;
    double r_temp= r_crit(f1_in,D2f0Dtheta2_in,kappa2_model,0.0,pert_params,_r_e);
    double r_temp2=r_crit(f1_in,D2f0Dtheta2_in,kappa2_model,half_pi,pert_params,_r_e);
    
    for(int i=npts; i>=1; i--){
      if(r_temp2<r_temp){ r_temp2=r_crit(f1_in,D2f0Dtheta2_in,kappa2_model,theta,pert_params,_r_e);}
    theta-=half_pi/npts;
    }
    double r_crit_max=r_temp2;
    double dr_max=r_crit_max-_r_e;
//     printf("the maximum deviation is %f \n",dr_max);
    return 3.0*fabs(c3_model)*pow((dr_max),2);

}

#endif
