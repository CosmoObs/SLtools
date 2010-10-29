/** @file
* Package useful to compute image properties, such number, center coordinates, curvature radius, length, width, length-to-width ratio.
*
* For more details, see Chap 2.5 of the Report on Perturbative Approach (in writeups folder)
*/

/** @package arc_properties
* Package to calculate several image properties.
*
*
* Detailde description FIXME
*
*/

#include <cstdlib>
#include <cstdio>
#include <fstream>        //funcoes de entrada e saida para arquivos
#include <iostream>
#include <math.h>

#ifndef PROPERTIES_ARCS_H
#define PROPERTIES_ARCS_H

#include "perturbative_method.h"
#include "theta_find.h"
# include "../numerical_methods/solving_system.h"

/** A function to calculate numerical integration for the following functions (It uses the Simpson´s Method )
*
* \f$ I_{1}= \int_{\theta_{\rm inf}}^{\theta_{\rm sup}}\bar{r}(\theta)d\theta \f$
*
* \f$ I_{2}= \int_{\theta_{\rm inf}}^{\theta_{\rm sup}}W(\theta)\bar{r}(\theta)d\theta \f$
*
* where \f$ W(\theta) \f$ and \f$ \bar{r}(\theta)\f$ are functions of \f$ f_1 \f$ and \f$ df_0/d\theta \f$. 
*
*  \param f1_in function related to the perturbed potential
*  \param Df0Dtheta_in function related to the perturbed potential
*  \param kappa2 \f$  \kappa_2 = 1 - \left[\frac{d^2 \phi_0(r)}{dr^2}\right]_{r=r_e} \f$
*  \param source_in a struct that define elliptical sources
*  \param pert_params[] a vector that contains all the perturbation parameters
*  \param _r_e Einstein Radius (\f$ R_{_{\mathrm{E}}}\f$ )
*  \param x_lim[2]: Vector containing the limits of the integral (x_lim[0]: lower limit, x_lim[1]: upper limit)
*  \param iflag :flag to decide which of the integrals will be performed (If iflag=1 performs \f$ I_1\f$, else \f$ I_2 \f$)
*  \param npts (number to divide the range of limits of the integral, if it isn´t given, the code assumes 2001. Note: must be a non number)
*  \return \f$ I_1 \f$ or  \f$ I_2 \f$
*  \sa wr_theta, r_mean
*
*/

double integrator(f_type f1_in, f_type Df0Dtheta_in, double kappa2, elliptical_source source_in, double pert_params[], double _r_e, double x_lim[], int iflag, int npts=2001){
  double x[npts], ftemp[npts], result;
  double dx=(x_lim[1]-x_lim[0])/(npts-1);
//   printf("the angle limits are %f %f\n",x_lim[0],x_lim[1]);
  for (int i=1; i<=npts;i++){
      x[i]=x_lim[0]+(i-1)*dx;
      if(iflag==1){
        ftemp[i]=r_mean(f1_in, Df0Dtheta_in, kappa2, source_in, x[i], pert_params, _r_e);
        }else{
        ftemp[i]=wr_theta(f1_in,Df0Dtheta_in, kappa2, source_in, x[i], pert_params, _r_e);}
      }
  int mpts=(npts-1)/2;
  double sum_1=0.0, sum_2=0.0, sum_3=0.0;

  for (int k=1; k<=mpts; k++){
      int k1=2*k-1, k2=2*k,k3=2*k+1;
      sum_1+=ftemp[k1];
      sum_2+=ftemp[k2];
      sum_3+=ftemp[k3];
    }
  result=(dx/3.0)*(sum_1+4.0*sum_2+sum_3);
  return result;  

}
/** Function to calculate the length of the image for the method that take into account the distance between 3 point: 
* 
* Knowing \f$ \theta_{\rm inf}\f$ and \f$ \theta_{\rm sup}\f$, we compute the mean value \f$ \theta_g= 0.5(\theta_{\rm inf}+\theta_{\rm sup})\f$
*
* For the angles above, we compute \f$ \bar{r}(\theta)\f$ and
* \f$ L(\theta_1,\theta_2)=\sqrt{\bar{r}_1^2+\bar{r}_2^2+\bar{r}_1\bar{r}_2\cos{(\theta_2-\theta_1)}}\f$.
*
*  \param f1_in function related to the perturbed potential
*  \param Df0Dtheta_in function related to the perturbed potential
*  \param source_in a struct that define elliptical sources
*  \param kappa2 \f$  \kappa_2 = 1 - \left[\frac{d^2 \phi_0(r)}{dr^2}\right]_{r=r_e} \f$
*  \param pert_params[] a vector that contains all the perturbation parameters
*  \param x_lim[2]: Vector containing the angle of the extrema of the image (x_lim[0]: lower angle, x_lim[1]: upper angle)
*  \param _r_e Einstein Radius (\f$ R_{_{\mathrm{E}}}\f$ )
*  \return \f$ \mathcal{L}_2= L(\theta_{\rm inf},\theta_g)+L(\theta_g,\theta_{\rm sup})\f$
*
*/
double method_L2(f_type f1_in, f_type Df0Dtheta_in, double kappa2, elliptical_source source_in, double x_lim[], double pert_params[],double _r_e){
  double theta_1=x_lim[0],theta_2=x_lim[1];
  double theta_g=0.5*(theta_1+theta_2);
  double r_mean1=r_mean(f1_in, Df0Dtheta_in, kappa2, source_in, theta_1, pert_params, _r_e);
  double r_meang=r_mean(f1_in, Df0Dtheta_in, kappa2, source_in, theta_g, pert_params, _r_e);
  double r_mean2=r_mean(f1_in, Df0Dtheta_in, kappa2, source_in, theta_2, pert_params, _r_e);
  double l1g=sqrt(pow(r_mean1,2)+pow(r_meang,2)-2*r_mean1*r_meang*cos(theta_g-theta_1));
  double l2g=sqrt(pow(r_mean2,2)+pow(r_meang,2)-2*r_mean2*r_meang*cos(theta_2-theta_g));
  return l1g+l2g;
}

/** Function to calculate the length of the image, taking into account the length of an arc of circunference passing for three points
* 
* Knowing \f$ \theta_{\rm inf}\f$ and \f$ \theta_{\rm sup}\f$, we compute the mean value \f$ \theta_g= 0.5(\theta_{\rm inf}+\theta_{\rm sup})\f$
*
* For the angles above, we compute \f$ \bar{r}(\theta)\f$, \f$(x_{1i},x_{2i})=(\bar{r}(\theta_i)\cos(\theta_i),\bar{r}(\theta_i)\sin(\theta_i))\f$ and  \f$ L_1=\sqrt{\bar{r}_{\rm inf}^2+\bar{r}_{\rm sup}^2+\bar{r}_{\rm inf}\bar{r}_{\rm sup}\cos{(\theta_{\rm sup}-\theta_{\rm inf})}}\f$.
*
* We solve the following system equation
* 
* \f$ x^2_{1i}+x^2_{2i}+Ax_{1i}+Bx_{2i}+C=0 \f$ (\f$i=1,2,3 \f$) 
*
* (i.e, we obtain \f$ A,B \f$ and \f$ C \f$)and theferore we obtain the center of the circle and its radius, passing for these three points, i.e.,
*
* \f$ x_c=-0.5A \f$, \f$ y_c = -0.5B \f$, \f$ R_c=\sqrt{x_c^2+y_c^2-C} \f$
*
*  \param f1_in function related to the perturbed potential
*  \param Df0Dtheta_in function related to the perturbed potential
*  \param kappa2 \f$  \kappa_2 = 1 - \left[\frac{d^2 \phi_0(r)}{dr^2}\right]_{r=R_{\rm E}} \f$
*  \param source_in a struct that define elliptical sources
*  \param theta_lim1:  Lower angle, i.e., \f$ \theta_{\rm inf}\f$ (It is obtained from theta_find.h function)
*  \param theta_lim2:  Upper angle, i.e., \f$ \theta_{\rm sup}\f$ (It is obtained from theta_find.h function)
*  \param pert_params[] a vector that contains all the perturbation parameters
*  \param circ_params[] : output vector containing \f$ circ\_par[0]= x_c \f$, \f$ circ\_par[1]= y_c \f$ and \f$ circ\_par[2]=R_c \f$
*  \param _r_e Einstein Radius (\f$ R_{_{\mathrm{E}}}\f$ )
*  \return \f$ \mathcal{L}_3= \arccos{\left(1-\dfrac{L_1^2}{2 R_c^2} \right)}R_c\f$
*
*/
double method_L3(f_type f1_in, f_type Df0Dtheta_in, double kappa2, elliptical_source source_in, double x_lim[], double pert_params[],double circ_par[],double _r_e){
  double theta_1=x_lim[0],theta_2=x_lim[1];
  double theta_g=0.5*(theta_1+theta_2);
  double r_mean1=r_mean(f1_in, Df0Dtheta_in, kappa2, source_in, theta_1, pert_params, _r_e);
  double r_meang=r_mean(f1_in, Df0Dtheta_in, kappa2, source_in, theta_g, pert_params, _r_e);
  double r_mean2=r_mean(f1_in, Df0Dtheta_in, kappa2, source_in, theta_2, pert_params, _r_e);
  double A[3][3];
  double B[3],X[3];
  A[0][0]=r_mean1*cos(theta_1), A[0][1]=r_mean1*sin(theta_1), A[0][2]=1.0, B[0]=-pow(r_mean1,2);
  A[1][0]=r_meang*cos(theta_g), A[1][1]=r_meang*sin(theta_g), A[1][2]=1.0, B[1]=-pow(r_meang,2);
  A[2][0]=r_mean2*cos(theta_2), A[2][1]=r_mean2*sin(theta_2), A[2][2]=1.0, B[2]=-pow(r_mean2,2);
  solving_system(A,B,X);
  double x_c=-0.5*X[0], y_c=-0.5*X[1], r_c=sqrt(pow(x_c,2)+pow(y_c,2)-X[2]);
  circ_par[0]=x_c, circ_par[1]=y_c, circ_par[2]=r_c;
  double d12=sqrt(pow(A[0][0]-A[2][0],2)+pow(A[0][1]-A[2][1],2));
  double arg_acos=1.0-(pow(d12,2)/(2.0*pow(r_c,2)));
  double  l3= r_c*acos(arg_acos) ;
  return l3;
}


/** A function to calculate the length of the image usig a  numerical integration of the function 
* \f$ I_{1}= \int_{\theta_{\rm inf}}^{\theta_{\rm sup}}\bar{r}(\theta)d\theta \f$
*
*  \param f1_in function related to the perturbed potential
*  \param Df0Dtheta_in function related to the perturbed potential
*  \param kappa2 \f$  \kappa_2 = 1 - \left[\frac{d^2 \phi_0(r)}{dr^2}\right]_{r=r_e} \f$
*  \param source_in a struct that define elliptical sources
*  \param pert_params[] a vector that contains all the perturbation parameters
*  \param _r_e Einstein Radius (\f$ R_{_{\mathrm{E}}}\f$ )
*  \param x_lim[2]: Vector containing the limits of the integral (x_lim[0]: lower limit, x_lim[1]: upper limit)
*  \param npts (number to divide the range of limits of the integral, if it isn´t given, the code assumes 2001. Note: must be a non number)
*  \return \f$ \mathcal{L}_4= I_1 \f$
*  \sa integrator
*
*/

double method_L4(f_type f1_in, f_type Df0Dtheta_in, double kappa2, elliptical_source source_in, double pert_params[], double _r_e, double x_lim[2], int npts=2001){
  int iflag1=1;
  double l4= integrator(f1_in,Df0Dtheta_in,kappa2, source_in, pert_params, _r_e, x_lim, iflag1, npts);
  return l4; 
}


/** A function to calculate the area of the image usig a  numerical integration of the function 
* \f$ I_{2}= \int_{\theta_{\rm inf}}^{\theta_{\rm sup}}W(\theta)\bar{r}(\theta)d\theta \f$
*
*  \param f1_in function related to the perturbed potential
*  \param Df0Dtheta_in function related to the perturbed potential
*  \param kappa2 \f$  \kappa_2 = 1 - \left[\frac{d^2 \phi_0(r)}{dr^2}\right]_{r=r_e} \f$
*  \param source_in a struct that define elliptical sources
*  \param pert_params[] a vector that contains all the perturbation parameters
*  \param _r_e Einstein Radius (\f$ R_{_{\mathrm{E}}}\f$ )
*  \param x_lim[2]: Vector containing the limits of the integral (x_lim[0]: lower limit, x_lim[1]: upper limit)
*  \param npts (number to divide the range of limits of the integral, if it isn´t given, the code assumes 2001. Note: must be a non number)
*  \return \f$ \mathcal{A}_k= I_2 \f$
*  \sa integrator
*
*/

double area_img(f_type f1_in, f_type Df0Dtheta_in, double kappa2, elliptical_source source_in, double pert_params[], double _r_e, double x_lim[], int npts=2001){
  int iflag2=2;
  double a_arc = integrator(f1_in, Df0Dtheta_in, kappa2, source_in, pert_params, _r_e, x_lim, iflag2, npts);
  return a_arc;
}

/** A function to calculate the width of the image usig a definition of the mean with
*
* Knowing \f$ \theta_{\rm inf}\f$ and \f$ \theta_{\rm sup}\f$, we compute the mean value \f$ \theta_g= 0.5(\theta_{\rm inf}+\theta_{\rm sup})\f$
* 
* For this mean angle, we compute the mean width of the image \f$ W(\theta_g)\f$.
*
*  \param f1_in function related to the perturbed potential
*  \param Df0Dtheta_in function related to the perturbed potential
*  \param kappa2 \f$  \kappa_2 = 1 - \left[\frac{d^2 \phi_0(r)}{dr^2}\right]_{r=r_e} \f$
*  \param source_in a struct that define elliptical sources
*  \param theta_lim1:  Lower angle, i.e., \f$ \theta_{\rm inf}\f$ (It is obtained from theta_find.h function)
*  \param theta_lim2:  Upper angle, i.e., \f$ \theta_{\rm sup}\f$ (It is obtained from theta_find.h function)
*  \param pert_params[] a vector that contains all the perturbation parameters
*  \param _r_e Einstein Radius (\f$ R_{_{\mathrm{E}}}\f$ )
*  \return \f$ W_q = W(\theta_g) \f$
*  \sa w_mean
*
*/

double width_q(f_type f1_in, f_type Df0Dtheta_in, double kappa2, elliptical_source source_in, double x_lim[], double pert_params[],double _r_e){
   double theta_1=x_lim[0],theta_2=x_lim[1];
   double theta_g=0.5*(theta_1+theta_2);
   return w_mean(f1_in, Df0Dtheta_in, kappa2, source_in, theta_g, pert_params, _r_e);
}
/** A function to calculate the  lenght-to-width ratio of an image using 
* \f$ \dfrac{L_k}{W_k};=\dfrac{\pi L_k^ 2}{4\mathcal{A}}\f$, where \f$ \mathcal{A} \f$ is the area of the image
*
*  \param f1_in function related to the perturbed potential
*  \param Df0Dtheta_in function related to the perturbed potential
*  \param kappa2 \f$  \kappa_2 = 1 - \left[\frac{d^2 \phi_0(r)}{dr^2}\right]_{r=r_e} \f$
*  \param source_in a struct that define elliptical sources
*  \param pert_params[] a vector that contains all the perturbation parameters
*  \param _r_e Einstein Radius (\f$ R_{_{\mathrm{E}}}\f$ )
*  \param x_lim[2]: Vector containing the limits of the integral (x_lim[0]: lower limit, x_lim[1]: upper limit)
*  \param iflag_k: flag control to calculate the \f$ L/W \f$ for the definition of \f$ L_k \f$
*   If(iflag_k=2,3,4, it uses \f$ L_2 \f$,\f$ L_3 \f$,\f$ L_4 \f$ )
*  \param npts (number to divide the range of limits of the integral, if it isn´t given, the code assumes 2001. Note: must be a non number)
*  \return \f$ L/W \f$
*  \sa area_img, method_L2, method_L3, method_L4
*
*/

double lw_k(f_type f1_in, f_type Df0Dtheta_in, double kappa2, elliptical_source source_in, double pert_params[], double _r_e, double x_lim[2], int iflag_k, int npts=2001){
  double length_temp;
  double twpi= 6.283185308, cte=twpi/8.; 
  double circ_par[3];
  if(iflag_k == 2){length_temp=method_L2(f1_in, Df0Dtheta_in, kappa2, source_in, x_lim, pert_params, _r_e);}
  if(iflag_k == 3){length_temp=method_L3(f1_in, Df0Dtheta_in, kappa2, source_in, x_lim, pert_params, circ_par, _r_e);}
  if(iflag_k == 4){length_temp=method_L4(f1_in, Df0Dtheta_in, kappa2, source_in, pert_params, _r_e, x_lim,npts);}
  return  (cte*pow(length_temp,2))/area_img(f1_in, Df0Dtheta_in, kappa2, source_in, pert_params, _r_e, x_lim, npts);
}
/** A function to calculate some image properties such:
* Number of images and for each image (identified by Image ID) is calculated
*
* coordinate center of the image \f$ x_c, y_c \f$ and its curvature radius \f$ R_c \f$
*
* Length (\f$ L_2, L_3, L_4 \f$)
*
* Width_q (\f$ W_q \f$)
*
* Length-to-width ratio using \f$W_q (we shall call as \f$ R_q \f$) \f$, i.,e \f$ R_{q,i}=L_i/W_q \f$ (\f$ i=2,3,4 \f$)  
*
* Length-to-width ratio using the area of the image (we shall call as \f$ R_k \f$) where \f$ R_k=\left(\dfrac{L}{W} \right)_k \f$ for (\f$ i=2,3,4 \f$) )
* 
*  \param f1_in function related to the perturbed potential
*  \param Df0Dtheta_in function related to the perturbed potential
*  \param kappa2 \f$  \kappa_2 = 1 - \left[\frac{d^2 \phi_0(r)}{dr^2}\right]_{r=r_e} \f$
*  \param source_in a struct that define elliptical sources
*  \param pert_params[] a vector that contains all the perturbation parameters
*  \param _r_e Einstein Radius (\f$ R_{_{\mathrm{E}}}\f$ )
*  \param npts (number to divide the range of limits of the integral, if it isn´t given, the code assumes 2001. Note: must be a non number)
*  \param ifile_in File to write all the results. If it isn´t given, the results will be written in the screen.
*  \return ifile_in 
*  \sa lw_k, width_q, method_L2, method_L3, method_L4
*/

void arc_measures(f_type f1_in, f_type Df0Dtheta_in, double kappa2, elliptical_source source_in, double pert_params[], double _r_e, int npts=2001, FILE *file_in=NULL){

  FILE *outcir = fopen ("circle_img.dat","w"); 
  double twpi= 6.283185308; 
  double arcs[10][2];
  double circ_par[3];
  double xarg_min,xarg_sup;
  theta_find(Df0Dtheta_in, source_in, pert_params,_r_e, arcs);
  int N_img=arcs[0][0];   
  double theta_inf[N_img];
  double theta_sup[N_img]; 
//
  double len_1,len_2, len_2b, len_3;  
  double w_q, lw_1,lw_2, lw_2b, lw_3; 
  double x_lim[2];

  if(file_in==NULL){
    printf("\t\tImage properties such number,angles, length, width and length-to-width ratios\n\n");
    printf("Number of images is %i\n\n",N_img);
  }else{
    
    fprintf(file_in,"\t\t Image properties such number, angles, lengths, widths and length-to-width ratios\n");
    fprintf(file_in,"Number of images is %i\n\n",N_img);
    }

  for(int i=1;i<=arcs[0][0];i++){
    xarg_min= arg_sqrt(Df0Dtheta_in, source_in, arcs[i][0], pert_params,_r_e);
    xarg_sup= arg_sqrt(Df0Dtheta_in, source_in, arcs[i][1], pert_params,_r_e);
    theta_inf[i]=arcs[i][0], theta_sup[i]=arcs[i][1]; 
    if( (theta_inf[i]> ((3./4.)*twpi)) && (theta_inf[i]< twpi) &&  (theta_sup[i]< (twpi/4.))){theta_sup[i]=theta_sup[i]+twpi;}
    x_lim[0]=theta_inf[i], x_lim[1]=theta_sup[i];  
// Computing the first measurement for the length, using the method L2 of the Report
    len_1= method_L2(f1_in,Df0Dtheta_in, kappa2, source_in, x_lim, pert_params,_r_e);
// Computing the first measurement for the length, using the method L3 of the Report
    len_2= method_L3(f1_in,Df0Dtheta_in, kappa2, source_in, x_lim, pert_params, circ_par, _r_e);
// Computing the first measurement for the length, using the method L4 of the Report
    len_3=method_L4(f1_in,Df0Dtheta_in, kappa2, source_in, pert_params, _r_e, x_lim,npts);
//  Measuring the widht of the arc using the Eq. 2.69 of the Report
    w_q=width_q(f1_in,Df0Dtheta_in, kappa2, source_in, x_lim, pert_params,_r_e);
//  Measuring the length-to-width ratio of the arcs using the formula of the area
    lw_1=lw_k(f1_in, Df0Dtheta_in, kappa2, source_in,pert_params, _r_e, x_lim,2,npts);
    lw_2=lw_k(f1_in, Df0Dtheta_in, kappa2, source_in,pert_params, _r_e, x_lim,3,npts);
    lw_3=lw_k(f1_in, Df0Dtheta_in, kappa2, source_in,pert_params, _r_e, x_lim,4,npts);
    double w1=len_1/lw_1,w2=len_2/lw_2,w3=len_3/lw_3;

    double theta2 = 0.0;
    int ncirc=200;
    for(int i1=0;i1<=ncirc;i1++){
            fprintf(outcir," %f %f\n",circ_par[0]+circ_par[2]*cos(theta2),circ_par[1]+circ_par[2]*sin(theta2));
            theta2+= 6.283185308/ncirc;
    }
    fprintf(outcir,"\n");
    
    
    if(file_in==NULL){
        printf("image ID: %i\n",i);
        printf("angle_inf %f arg_sqrt_inf:  %f angle_sup:  %f arg_sqrt_sup: %f\n",theta_inf[i],xarg_min,theta_sup[i],xarg_sup);
	printf("image center coordinates (x0,y0): (%f, %f) curvature radius: %f\n",circ_par[0],circ_par[1],circ_par[2]);
	printf("length_1: %f, length_2: %f, length_3: %f\n",len_1,len_2,len_3 );
        printf("width_q: %f, width_1: %f, width_2: %f, width_3 %f\n",w_q,w1,w2,w3);
        printf("LW_q1: %f, LW_q2: %f, LW_q3 %f\n",len_1/w_q,len_2/w_q,len_3/w_q);
        printf("LW_k1: %f, LW_k2: %f,  LW_k3 %f\n\n",lw_1,lw_2,lw_3);
    } 
    else {
        fprintf(file_in,"image ID: %i\n",i);
        fprintf(file_in,"angle_inf %f arg_sqrt_inf:  %f angle_sup:  %f arg_sqrt_sup: %f\n",theta_inf[i],xarg_min,theta_sup[i],xarg_sup);  
	fprintf(file_in,"image center coordinates (x0,y0): (%f, %f), curvature radius: %f\n",circ_par[0],circ_par[1],circ_par[2]);
        fprintf(file_in,"length_1: %f, length_2: %f, length_3: %f\n",len_1,len_2,len_3 );
        fprintf(file_in,"width_q: %f, width_1: %f, width_2: %f, width_3 %f\n",w_q,w1,w2,w3);
        fprintf(file_in,"LW_q1: %f, LW_q2: %f, LW_q3 %f\n",len_1/w_q,len_2/w_q,len_3/w_q);
        fprintf(file_in,"LW_k1: %f, LW_k2: %f,  LW_k3 %f\n\n",lw_1,lw_2,lw_3);
	
    }
      

  }

  
}

# endif