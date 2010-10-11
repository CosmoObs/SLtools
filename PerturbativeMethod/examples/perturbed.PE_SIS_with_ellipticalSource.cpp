/** @file
 * Code to compute arc points for SIS central potential and perturber, with an elliptical source
*/

/** @package perturbedSIE_with_ellipticalSource 
*  Code to compute quantities related to the Perturbative Method discussed in  http://adsabs.harvard.edu/abs/2007MNRAS.382L..58A and further works 
*
*  Detailed description FIXME
*
*/

// Orig by G. Caminha
// Sept 2010
// Modifications and comments by MSSG

#include <cstdio>
#include <cmath> 
#include <cstdlib>
#include <iostream>

#include "../perturbative_method.h"
#include "../sis_sub_model.h" 

//To compile: g++ -Wall perturbedSIE_with_ellipticalSource.cpp


//Quantities related to the potential
/*********************************************************************************************************************/
//!  Function defined in http://arxiv.org/abs/0804.4277, for SIE lens model
//!
//!  \f$ f_{n}(\theta) \equiv \frac{1}{n!} \left[\frac{\partial\psi(r,\theta)}{\partial r}\right]_{r=1} \f$
/*!
  \param theta angular coordinate in the lens plane
  \param pot_params[] potential parameters
  \return 0.0
*/

//double f1_SIE(double theta, double pot_params[]){
//  return 0.0;
// }

//double Df0Dtheta_SIE(double theta, double pot_params[]){
//  return 0.0;
//}

/*********************************************************************************************************************/

// ********************** Eq. 8, Alard 2008
//pot_params[0] = mp
//pot_params[1] = rp
//pot_params[2] = thetap
double f1_SIE(double theta, double pot_params[]){
  double mp=pot_params[0],rp=pot_params[1],thetap=pot_params[2], eta=pot_params[3];
  double bar_theta=theta-thetap;
  double Rsq=1.0-2.0*rp*cos(bar_theta)+pow(rp,2);
  double denr = sqrt(Rsq) ; 
  double numr = mp*(1 - rp*cos(bar_theta)) ;
  return -(eta/2)*cos(2.0*theta)+numr/denr; 
}

double Df0Dtheta_SIE(double theta, double pot_params[]){
  double mp=pot_params[0],rp=pot_params[1],thetap=pot_params[2], eta=pot_params[3];
  double bar_theta=theta-thetap;
  double Rsq=1.0-2.0*rp*cos(bar_theta)+pow(rp,2);
  double numr = mp*rp*sin(bar_theta);
  double denr = sqrt(Rsq) ; 
  return eta*sin(2.0*theta)+ numr/denr; 
}


double D2f0Dtheta2_SIE(double theta, double pot_params[]){
double mp=pot_params[0],rp=pot_params[1],thetap=pot_params[2], eta=pot_params[3];
double bar_theta=theta-thetap;
double Rsq=1.0-2.0*rp*cos(bar_theta)+pow(rp,2);
return 2.0*eta*cos(2.0*theta)+(mp*rp/sqrt(Rsq))*(cos(bar_theta)-rp*(pow(sin(bar_theta),2))/Rsq);
}



/*********************************************************************************************************************/


void print_arcs_SIE(elliptical_source source_in, int npts,  double mp, double rp, double thetap, double central_eta){
  
  double theta = 0.0;
  double r_p = 0.0;
  double r_m = 0.0;
  double pot_params[4];
  //  pot_params[0] = 0.005;
  //  pot_params[1] = 1.3;
  //  pot_params[2] = 0.0;	
//  pot_params[2] = 6.283185308/20.;
//  pot_params[3]=0.1; //this correspond to the ellipticity of the central lensing potential

  pot_params[0] = mp;
  pot_params[1] = rp;
  pot_params[2] = thetap;	
  pot_params[3]=central_eta; //this correspond to the ellipticity of the central lensing potential

// Added by Habib Dumet-Montoya
  double r_tcrit=0.0;
  double y1_caust=0.0;
  double y2_caust=0.0;
  double src_y1=0.0;
  double src_y2=0.0;

FILE *outarcfile = fopen ("arcs_sis.dat" , "w");
FILE *outcritcurvefile = fopen ("tang_crit.dat" , "w");
FILE *outcausticfile = fopen ("tang_caust.dat" , "w");
FILE *outsrcfile = fopen ("src_plot.dat" , "w");


  for(int i=0;i<=npts;i++){

    if(arg_sqrt(Df0Dtheta_SIE, source_in, theta, pot_params)>0.0){
      r_p = 1.0+dr_plus(f1_SIE, Df0Dtheta_SIE, 1.0, source_in, theta, pot_params);
      r_m = 1.0+dr_minus(f1_SIE, Df0Dtheta_SIE, 1.0, source_in, theta, pot_params);
//      printf("%E %E\n",r_p*cos(theta),r_p*sin(theta));
//      printf("%E %E\n",r_m*cos(theta),r_m*sin(theta));

	fprintf(outarcfile,"%E %E\n",r_p*cos(theta),r_p*sin(theta));
        fprintf(outarcfile,"%E %E\n",r_m*cos(theta),r_m*sin(theta));	
    }

// Delete later..
//  Defining the tangential critical radius
     r_tcrit= r_crit(f1_SIE,D2f0Dtheta2_SIE,1.0,theta,pot_params); 
// Writing in a file the tangential critical curve
     fprintf(outcritcurvefile,"%E %E\n",r_tcrit*cos(theta),r_tcrit*sin(theta));   
//  Defining the polar equation of the tangential caustic
     y1_caust=caustic_y1(Df0Dtheta_SIE,D2f0Dtheta2_SIE,theta,pot_params);
     y2_caust=caustic_y2(Df0Dtheta_SIE,D2f0Dtheta2_SIE,theta,pot_params);
//   Writing in a file the tangential caustic line
     fprintf(outcausticfile,"%E %E\n",y1_caust,y2_caust); 
//   Writing the source in a file 
     fprintf(outsrcfile,"%E %E\n",y1_src(source_in, theta),y2_src(source_in, theta));  

    theta+= 6.283185308/npts;
  }


}

int main(int argc, char *argv[]){

  int npts;

  elliptical_source source;

  double     mp;
  double     rp ;
  double     thetap;
  
  double     central_eta;


  if (argc < 2) {
    // printf("No input args, using defaults; To input, use this format: \n ");
    // printf("Num pts   Source parameters: x0   y0  R0 (size)  eta (ellip)   theta (ang. position)\n");

    npts = 4000;
    
    source.x0 = 0.00;
    source.y0 = 0.10;
    source.R0 = 0.03;
    source.eta0 = 0.0;
    source.theta0 =  0.0;
    
  }
  else{
    npts = atoi(argv[1]);

    source.x0 = atof(argv[2]);
    source.y0 = atof(argv[3]);
    source.R0 = atof(argv[4]);
    source.eta0 = atof(argv[5]);
    source.theta0 =  atof(argv[6]);

    source.x0 = atof(argv[2]);
    source.y0 = atof(argv[3]);
    source.R0 = atof(argv[4]);
    source.eta0 = atof(argv[5]);
    source.theta0 =  atof(argv[6]);

     mp = atof(argv[7]) ;
     rp = atof(argv[8]);
     thetap = atof(argv[9]);

     central_eta = atof(argv[10]);
  }

  /*   printf(" Using: npts %i  \n", npts);
   printf(" x0,y0 = (%f,%f) \n", source.x0,source.y0);
   printf(" R0, eta0, theta0 = (%f,%f,%f) \n",source.R0, source.eta0, source.theta0 );
  */

  print_arcs_SIE(source, npts, mp,rp,thetap, central_eta);

  //  double pot_params[] = {0.03, 1.3, 0.0};
  //  printf("%f\n",f1_SIE(3.14159/2, pot_params));
  return 0;

}
