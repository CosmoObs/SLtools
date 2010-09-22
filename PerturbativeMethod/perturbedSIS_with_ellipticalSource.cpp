/** @file
 * Code to compute arc points for SIS central potential and perturber, with an elliptical source
*/

/** @package perturbedSIS_with_ellipticalSource 
*  Code to compute quantities related to the Perturbative Method discussed in  http://adsabs.harvard.edu/abs/2007MNRAS.382L..58A and further works 
*
*  Detailed description FIXME
*
*/

// Orig by G. Caminha
// Sept 2010
// Modifications and comments by MSS Gill

#include <cstdio>
#include <cmath> 
#include <cstdlib>
#include <iostream>

#include "perturbative_method_functions.h"  // This is needed to define the elliptical_source structure, and get the arg of the sqrt for the arc point solns

//To compile: g++ -Wall perturbedSIS_with_ellipticalSource.cpp


//Quantities related to the potential
/*********************************************************************************************************************/
//!  Function defined in http://arxiv.org/abs/0804.4277 (Peirani et al.) or Alard 2008, for SIS lens model
//!
//!  \f$ f_{n}(\theta) \equiv \frac{1}{n!} \left[\frac{\partial\psi(r,\theta)}{\partial r}\right]_{r=1} \f$
/*!
  \param theta angular coordinate in the lens plane
  \param pot_params[] potential parameters
  \return 0.0
*/

/*********************************************************************************************************************/

double twopi = 6.283185308;

/////// Now define the two f(perturber) functions we'll need to calculate the arc properties -- eq. 9 in Alard08

//pot_params[0] = mp in units of central potential mass
//pot_params[1] = rp distance to perturber in units of r_E
//pot_params[2] = thetap angle of perturber
double f1_SIS(double theta, double pot_params[]){
  double numr = pot_params[0] * (1.0 - pot_params[1]*cos(theta-pot_params[2]) );
  double denr = sqrt ( 1.0 - 2.0*pot_params[1]*cos(theta-pot_params[2]) + pot_params[1]*pot_params[1]);
  return numr/denr; 
}


double Df0Dtheta_SIS(double theta, double pot_params[]){
  double numr = pot_params[0] * (pot_params[1]*sin(theta-pot_params[2]) );
  double denr = sqrt ( 1.0 - 2.0*pot_params[1]*cos(theta-pot_params[2]) + pot_params[1]*pot_params[1]);
  return numr/denr;
}




/*********************************************************************************************************************/

void print_arc_points(elliptical_source source_in, int npts,  double mp, double rp, double thetap){

  
  // Assign properties of the perturber into the pot_params vector
 double pot_params[3];
 pot_params[0] = mp;
 pot_params[1] = rp;
 pot_params[2] = thetap;


 double theta = 0.0; // This is the general angle for where we are evaluating things in the lens plane
  double r_p = 0.0;  // Outer arc edge value 
  double r_m = 0.0;  // Inner arc edge value 

  double kappa2 = 1; // Always unity for SIS potentials

  // Run over all theta, cranking out the points
  for(int i=0;i<=npts;i++){

    if(arg_sqrt(Df0Dtheta_SIS, source_in, theta, pot_params)>0.0){ // Df0Dtheta_SIS is defined in  perturbative_method_functions.h (pmf.h)
      r_p = 1.0+x_plus(f1_SIS, Df0Dtheta_SIS, kappa2, source_in, theta, pot_params); // Outer arc edge value (x_plus is defined in pmf.h too)
      r_m = 1.0+x_minus(f1_SIS, Df0Dtheta_SIS, kappa2, source_in, theta, pot_params);// Inner arc edge value (x_minus is defined in pmf.h too)
      printf("%E %E\n",r_p*cos(theta),r_p*sin(theta)); // Print outer arc edge value in Cartesian coords
      printf("%E %E\n",r_m*cos(theta),r_m*sin(theta)); // Print inner arc edge value in Cartesian coords
    }
    theta+= twopi/npts;
  }

}

int main(int argc, char *argv[]){

  int npts;
  
  elliptical_source source;  // This structure comes from perturbative_method.h
  
  // Perturber properties
  double mp = 0.1 ; //  mass of perturber in units of central potential mass
  double rp = 0.1 ; //  distance to perturber in units of r_E
  double thetap = 0;//  thetap angle of perturber    
  
  if (argc < 2) {
    // printf("No input args, using defaults; To input, use this format: \n ");
    // printf(" Num pts around 2pi   Source parameters: x0   y0  R0 (size)  eta (ellip)   theta (ang. position)\n");
    // printf(" Perturber params: mp   rp   thetap  \n");

    npts = 100;
    
    // Source defaults
    source.x0 = 1;
    source.y0 = 0.0;
    source.R0 = 0.9;
    source.eta0 = 0.0;
    source.theta0 =  0.0;
    // Perturber property defaults
    mp = 0.1 ;
    rp = 0.1 ;
    thetap = 0;

  }
  else{ // Get all the properties from the command line
    npts = atoi(argv[1]);

    source.x0 = atof(argv[2]);
    source.y0 = atof(argv[3]);
    source.R0 = atof(argv[4]);
    source.eta0 = atof(argv[5]);
    source.theta0 =  atof(argv[6]);

     mp = atof(argv[7]) ;
     rp = atof(argv[8]);
     thetap = atof(argv[9]);
  }

  /*   printf(" Using: npts %i  \n", npts);
   printf(" x0,y0 = (%f,%f) \n", source.x0,source.y0);
   printf(" R0, eta0, theta0 = (%f,%f,%f) \n",source.R0, source.eta0, source.theta0 );
  */



// Run over entire all theta values in image plane, in nsteps steps
// Outputting just the values of the Einstein radius 
  printf("\n\n");
  double theta = 0.0;
  double R_E  = 1.0;
  for(int i=0;i<=npts;i++){
    printf("%E %E\n",R_E*cos(theta), R_E*sin(theta));
    theta += twopi/npts;
  }

  // Call the function that will print out the arc points
  print_arc_points(source, npts, mp,rp,thetap);

  return 0;
}
