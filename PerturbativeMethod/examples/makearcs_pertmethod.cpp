/* Arc making code to replicate MM (MMakler) Math'ca code
Started by G. Caminha on  Aug 31, 2010
Commented and small changes by MSSG
 */

// Needed headers, in the g++ vs. gcc style
#include <cstdio>
#include <cmath> 
#include <ctime>
#include <iostream>

/********************************** Define needed functions *************************************/

// R_Einstein is implicitly set to unity throughout this code so far (prob missing in a few places)

double R_E = 1.0;

// All the below uses Alard 2007 notation

// Define the f0 and f1 functions, 0th and 1st order derivatives of psi, the perturbed part of the potential
// For SIS, psi = 0, so f0 = f1 = 0
double f0(double theta){
  return 0.0;
}

double f1(double theta){
  return 0.0;
}

// Kappa2 = 1-C_2, where C_2 = 1/2 * d^2 phi_0 / dr^2 (phi_0 =
// axisymmetric part of the unperturbed potential).
// Since phi_0 = constant * r for SIS, C_2 = 0, and Kappa2 = 1
double k2(){
  return 1.0;
}

// Now take the source center off-center in the source plane, to the
// point (x0, y0), this modifies the f0,1 functions generally into
// f0,1_bar functions
double f0_bar(double theta, double x_0, double y_0){
  return f0(theta) + x_0*cos(theta) + y_0*sin(theta);
}

double f1_bar(double theta, double x_0, double y_0){
  return f1(theta) + x_0*cos(theta) + y_0*sin(theta);
}


// Define the 1st deriv in theta of the f0bar funct
// SIS: f0 = 0, f0bar = x_0*cos(theta) + y_0*sin(theta)
double Df0_bar_Dtheta(double theta, double x_0, double y_0){
  return y_0*cos(theta) - x_0*sin(theta);
}

// Now get the arc radii in the image plane, which we decided we would
// call in the better notation "x" vs. "dr" as Alard -- corresp. to
// Alard eq. 12.   R0 is the size of the circular source

// This is the outer arc radius 
double x_plus(double theta, double x_0, double y_0, double R0){
  double DfbDt = Df0_bar_Dtheta(theta, x_0, y_0);
  double sqrtarg = R0*R0 - DfbDt*DfbDt;
  if (sqrtarg < 0.0) return 0.0;

  return 1.0/k2() * ( f1_bar(theta,x_0,y_0) +  sqrt(sqrtarg) );
}

// This is the inner arc radius 
double x_minus(double theta, double x_0, double y_0, double R0){
  double DfbDt = Df0_bar_Dtheta(theta, x_0, y_0);
  double sqrtarg = R0*R0 - DfbDt*DfbDt;
  if (sqrtarg < 0.0) return 0.0;

  return 1.0/k2() * ( f1_bar(theta,x_0,y_0) -  sqrt(sqrtarg) );
}

// How many steps to break the circle up into
int nsteps = 500;

/**********************************  Begin main code ******************************/
int main(){

  // Define some numbers, same as what MM has in his code
  double x_0 = 1;
  double y_0 = 0 ; 
  double R_0 = 1.0/15.0;

  // We'll need this
  double theta = 0.0;

//  printf("%f\n",x_plus(theta , x_0, y_0, R_0));  // Keep in -- to check selected numerical values vs. MM code

// Run over entire all theta values in image plane, in nsteps steps
// Outputting inner radius values in Cartesian coords
  for(int i=0;i<=nsteps;i++){
    //if( x_minus(theta , x_0, y_0, R_0) > 0.0 ){
      printf("%E %E\n",(x_minus(theta , x_0, y_0, R_0)+1)*cos(theta), (x_minus(theta , x_0, y_0, R_0)+1)*sin(theta));
    //}
      theta += 2*3.1415/nsteps ;
  }

// Run over entire all theta values in image plane, in nsteps steps
// Outputting outer radius values in Cartesian coords
  printf("\n\n");
  theta = 0.0;
  for(int i=0;i<=nsteps;i++){
    //if( x_plus(theta , x_0, y_0, R_0) > 0.0 ){
      printf("%E %E\n",(x_plus(theta , x_0, y_0, R_0)+1)*cos(theta), (x_plus(theta , x_0, y_0, R_0)+1)*sin(theta));
    //}
      theta += 2*3.1415/nsteps ;
  }

// Run over entire all theta values in image plane, in nsteps steps
// Outputting just the values of the Einstein radius 
  printf("\n\n");
  theta = 0.0;
  for(int i=0;i<=nsteps;i++){
    printf("%E %E\n",R_E*cos(theta), R_E*sin(theta));
    theta += 2*3.1415/100.0;
  }

  // End code
  return 0;
}
