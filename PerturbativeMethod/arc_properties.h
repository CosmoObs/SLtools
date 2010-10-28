#include <cstdlib>
#include <cstdio>
#include <fstream>        //funcoes de entrada e saida para arquivos
//#include <cstdio>
#include <iostream>
#include <math.h>

#ifndef PROPERTIES_ARCS_H
#define PROPERTIES_ARCS_H

#include "perturbative_method.h"
# include "../numerical_methods/solving_system.h"

double integrator(f_type f1_in, f_type Df0Dtheta_in, double kappa2, elliptical_source source_in, double pert_params[], double _r_e, double x_lim[2], int iflag, int npts=2001){
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
//   printf("the sums are %f %f and %f \n", sum_1, sum_2, sum_3);
  result=(dx/3.0)*(sum_1+4.0*sum_2+sum_3);
//   printf("The result of the integral is %f\n",result);
  return result;  

}


double method_L2(f_type f1_in, f_type Df0Dtheta_in, double kappa2, elliptical_source source_in, double theta_lim1, double theta_lim2, double pert_params[],double _r_e){
  double theta_1=theta_lim1,theta_2=theta_lim2;
  double theta_g=0.5*(theta_1+theta_2);
  double r_mean1=r_mean(f1_in, Df0Dtheta_in, kappa2, source_in, theta_1, pert_params, _r_e);
  double r_meang=r_mean(f1_in, Df0Dtheta_in, kappa2, source_in, theta_g, pert_params, _r_e);
  double r_mean2=r_mean(f1_in, Df0Dtheta_in, kappa2, source_in, theta_2, pert_params, _r_e);
//   printf("some results %f %f %f %f %f %f\n",theta_1,r_mean1,theta_g,r_meang,theta_2,r_mean2);
  double l1g=sqrt(pow(r_mean1,2)+pow(r_meang,2)-2*r_mean1*r_meang*cos(theta_g-theta_1));
  double l2g=sqrt(pow(r_mean2,2)+pow(r_meang,2)-2*r_mean2*r_meang*cos(theta_2-theta_g));
  return l1g+l2g;
}

double method_L3(f_type f1_in, f_type Df0Dtheta_in, double kappa2, elliptical_source source_in, double theta_lim1, double theta_lim2, double pert_params[],double _r_e){
  double theta_1=theta_lim1,theta_2=theta_lim2;
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
//   double l3=r_c*fabs(theta_2-theta_1); 
  double d12=sqrt(pow(A[0][0]-A[2][0],2)+pow(A[0][1]-A[2][1],2));
  double arg_acos=1.0-(pow(d12,2)/(2.0*pow(r_c,2)));
  double  l3= r_c*acos(arg_acos) ;
  return l3;
}


double method_L3b(f_type f1_in, f_type Df0Dtheta_in, double kappa2, elliptical_source source_in, double theta_lim1, double theta_lim2, double pert_params[],double _r_e){
  double theta_1=theta_lim1,theta_2=theta_lim2;
  double theta_g=0.5*(theta_1+theta_2);
  double r_mean1=r_mean(f1_in, Df0Dtheta_in, kappa2, source_in, theta_1, pert_params, _r_e);
  double r_meang=r_mean(f1_in, Df0Dtheta_in, kappa2, source_in, theta_g, pert_params, _r_e);
  double r_mean2=r_mean(f1_in, Df0Dtheta_in, kappa2, source_in, theta_2, pert_params, _r_e);
  double x1=r_mean1*cos(theta_1), y1=r_mean1*sin(theta_1);
  double x2=r_mean2*cos(theta_2), y2=r_mean2*sin(theta_2);
  double x3=r_meang*cos(theta_g), y3=r_meang*sin(theta_g);
  double d12=sqrt(pow(x1-x2,2)+pow(y1-y2,2));
  double d23=sqrt(pow(x2-x3,2)+pow(y2-y3,2));  
  double d13=sqrt(pow(x3-x1,2)+pow(y3-y1,2));
  double vec_temp=fabs( (x1-x2)*(y2-y3) - (x2-x3)*(y1-y2) );
  double radius=(d12*d23*d13)/(2.*vec_temp);
//   double  l3_temp=radius*fabs(theta_2-theta_1);
  double arg_acos=1.0-(pow(d12,2)/(2.0*pow(radius,2)));
  double  l3_temp= radius*acos(arg_acos) ;
// L3 = arccos(1 - (L1**2)/(2*(R**2)) )*R 
  return l3_temp;
}

double method_L4(f_type f1_in, f_type Df0Dtheta_in, double kappa2, elliptical_source source_in, double pert_params[], double _r_e, double x_lim[2], int npts=2001){
  int iflag1=1;
  double l4= integrator(f1_in,Df0Dtheta_in,kappa2, source_in, pert_params, _r_e, x_lim, iflag1, npts);
  return l4; 
}

double area_img(f_type f1_in, f_type Df0Dtheta_in, double kappa2, elliptical_source source_in, double pert_params[], double _r_e, double x_lim[2], int npts=2001){
  int iflag2=2;
  double a_arc = integrator(f1_in, Df0Dtheta_in, kappa2, source_in, pert_params, _r_e, x_lim, iflag2, npts);
  return a_arc;
}

double width_q(f_type f1_in, f_type Df0Dtheta_in, double kappa2, elliptical_source source_in, double theta_lim1, double theta_lim2, double pert_params[],double _r_e){
   double theta_1=theta_lim1,theta_2=theta_lim2;
   double theta_g=0.5*(theta_1+theta_2);
   return w_mean(f1_in, Df0Dtheta_in, kappa2, source_in, theta_g, pert_params, _r_e);
}

double lw_k(f_type f1_in, f_type Df0Dtheta_in, double kappa2, elliptical_source source_in, double pert_params[], double _r_e, double x_lim[2], int iflag_k, int npts=2001){
  double length_temp;
  double twpi= 6.283185308, cte=twpi/8.; 
  if(iflag_k == 1){length_temp=method_L2(f1_in, Df0Dtheta_in, kappa2, source_in, x_lim[0], x_lim[1], pert_params, _r_e);}
  if(iflag_k == 2){length_temp=method_L3(f1_in, Df0Dtheta_in, kappa2, source_in, x_lim[0], x_lim[1], pert_params, _r_e);}
  if(iflag_k == 3){length_temp=method_L3b(f1_in, Df0Dtheta_in, kappa2, source_in, x_lim[0], x_lim[1], pert_params, _r_e);}
  if(iflag_k == 4){length_temp=method_L4(f1_in, Df0Dtheta_in, kappa2, source_in, pert_params, _r_e, x_lim,npts);}

  return  (cte*pow(length_temp,2))/area_img(f1_in, Df0Dtheta_in, kappa2, source_in, pert_params, _r_e, x_lim, npts);
}

void arc_measures(f_type f1_in, f_type Df0Dtheta_in, double kappa2, elliptical_source source_in, double pert_params[], double _r_e, int npts=2001, FILE *file_in=NULL){

  double twpi= 6.283185308; 
  double arcs[10][2];
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
    len_1= method_L2(f1_in,Df0Dtheta_in, kappa2, source_in, theta_inf[i],theta_sup[i], pert_params,_r_e);
// Computing the first measurement for the length, using the method L3 of the Report
    len_2= method_L3(f1_in,Df0Dtheta_in, kappa2, source_in, theta_inf[i],theta_sup[i], pert_params,_r_e);
// Computing the first measurement for the length, using the method L3 of the Report, using Pedro's Python Function
    len_2b= method_L3b(f1_in,Df0Dtheta_in, kappa2, source_in, theta_inf[i],theta_sup[i], pert_params,_r_e);
// Computing the first measurement for the length, using the method L4 of the Report
    len_3=method_L4(f1_in,Df0Dtheta_in, kappa2, source_in, pert_params, _r_e, x_lim,npts);
//  Measuring the widht of the arc using the Eq. 2.69 of the Report
    w_q=width_q(f1_in,Df0Dtheta_in, kappa2, source_in, theta_inf[i],theta_sup[i], pert_params,_r_e);
//  Measuring the length-to-width ratio of the arcs using the formula of the area
    lw_1=lw_k(f1_in, Df0Dtheta_in, kappa2, source_in,pert_params, _r_e, x_lim,1,npts);
    lw_2=lw_k(f1_in, Df0Dtheta_in, kappa2, source_in,pert_params, _r_e, x_lim,2,npts);
    lw_2b=lw_k(f1_in, Df0Dtheta_in, kappa2, source_in,pert_params, _r_e, x_lim,3,npts);
    lw_3=lw_k(f1_in, Df0Dtheta_in, kappa2, source_in,pert_params, _r_e, x_lim,4,npts);
    double w1=len_1/lw_1,w2=len_2/lw_2,w2b=len_2b/lw_2b,w3=len_3/lw_3;

    if(file_in==NULL){
        printf("image ID: %i\n",i);
        printf("angle_inf %f arg_sqrt_inf:  %f angle_sup:  %f arg_sqrt_sup: %f\n",theta_inf[i],xarg_min,theta_sup[i],xarg_sup);
        printf("length_1: %f, length_2: %f, length_2(b): %f, length_3: %f\n",len_1,len_2,len_2b,len_3 );
        printf("width_q: %f, width_1: %f, width_2: %f, width_2(b): %f, width_3 %f\n",w_q,w1,w2,w2b,w3);
        printf("LW_q1: %f, LW_q2: %f, LW_q2(b): %f, LW_q3 %f\n",len_1/w_q,len_2/w_q,len_2b/w_q,len_3/w_q);
        printf("LW_k1: %f, LW_k2: %f, LW_k2(b): %f, LW_k3 %f\n\n",lw_1,lw_2,lw_2b,lw_3);
    } 
    else {
        fprintf(file_in,"image ID: %i\n",i);
        fprintf(file_in,"angle_inf %f arg_sqrt_inf:  %f angle_sup:  %f arg_sqrt_sup: %f\n",theta_inf[i],xarg_min,theta_sup[i],xarg_sup);  
        fprintf(file_in,"length_1: %f, length_2: %f, length_2(b): %f, length_3: %f\n",len_1,len_2,len_2b,len_3);
        fprintf(file_in,"width_q: %f, width_1: %f, width_2: %f, width_2(b): %f, width_3 %f\n",w_q,w1,w2,w2b,w3);
        fprintf(file_in,"LW_q1: %f, LW_q2: %f, LW_q2(b): %f, LW_q3 %f\n",len_1/w_q,len_2/w_q,len_2b/w_q,len_3/w_q);
        fprintf(file_in,"LW_k1: %f, LW_k2: %f, LW_k2(b): %f, LW_k3 %f\n\n",lw_1,lw_2,lw_2b,lw_3);
    }
      

  }

  
}

# endif