#include <cstdlib>
#include <cstdio>

#include "../perturbative_method.h"
#include "../pnfw_model.h"
#include "../theta_find.h"

void print_arcs_nfw(elliptical_source source_in, double pot_params[],  int npts=200){

  double theta = 0.0;
  double r_p = 0.0;
  double r_m = 0.0;
// Added by Habib Dumet-Montoya
  double r_tcrit=0.0;
  double y1_caust=0.0;
  double y2_caust=0.0;
  double src_y1=0.0;
  double src_y2=0.0;

//pot_paramas[0]=r_s, pot_params[1]=ks,pot_params[2]=elipticity,
  //double pot_params[] = {0.9,1.0,0.2};
// Added by Habib Dumet-Montoya
FILE *outarc = fopen ("arcs_pnfw.dat" , "w");
FILE *oute = fopen ("einstein_circle.dat" , "w");
FILE *outtc = fopen ("tang_crit.dat" , "w");
FILE *outcau = fopen ("tang_caust.dat" , "w");
FILE *outsrc = fopen ("src_plot.dat" , "w");

    for(int i=0;i<=npts;i++){

      if(arg_sqrt(Df0Dtheta_nfw, source_in, theta, pot_params)>0.0){
      r_p = _r_e+pot_params[0]*dr_plus(f1_nfw, Df0Dtheta_nfw, kappa2_nfw(pot_params), source_in, theta, pot_params);
      r_m = _r_e+pot_params[0]*dr_minus(f1_nfw, Df0Dtheta_nfw,kappa2_nfw(pot_params), source_in, theta, pot_params);
//       printf("%E %E\n",r_p*cos(theta),r_p*sin(theta));
//       printf("%E %E\n",r_m*cos(theta),r_m*sin(theta));
      fprintf(outarc,"%E %E\n",r_p*cos(theta),r_p*sin(theta)); // printing to a file
      fprintf(outarc,"%E %E\n",r_m*cos(theta),r_m*sin(theta)); // printing to a file
    }
// Added by Habib: The tangential Critical Curve.
// Printing the Einstein Ring to a File
     fprintf(oute,"%E %E\n",_r_e*cos(theta),_r_e*sin(theta)); 
//  Defining the tangential critical radius
     r_tcrit= r_crit(f1_nfw,D2f0Dtheta2_nfw,kappa2_nfw(pot_params),theta,pot_params); 
// Writing in a file the tangential critical curve
     fprintf(outtc,"%E %E\n",r_tcrit*cos(theta),r_tcrit*sin(theta));   
//  Defining the polar equation of the tangential caustic
     y1_caust=caustic_y1(Df0Dtheta_nfw,D2f0Dtheta2_nfw,theta,pot_params);
     y2_caust=caustic_y2(Df0Dtheta_nfw,D2f0Dtheta2_nfw,theta,pot_params);
//   Writing in a file the tangential caustic line
     fprintf(outcau,"%E %E\n",y1_caust,y2_caust); 
//   Writing the source in a file 
     fprintf(outsrc,"%E %E\n",y1_src(source_in, theta),y2_src(source_in, theta));  
//
    theta+= 6.283185308/npts;
  }
}

int main(){

// Added by Habib Dumet-Montoya
  FILE *outls = fopen ("lensing_data.dat" , "w");  

  int npts = 1000;

  double twpi= 6.283185308;  
  double pot_params[] = {1.0,0.8,0.1};

  elliptical_source source;
/*
  source.x0 = 0.04;
  source.y0 = 0.04;
  source.R0 = 0.01;
  source.eta0 = 0.0;
  source.theta0 = 0.0;
*/
  source.x0 = (r_e_nfw(pot_params)/10);
  source.y0 = 0.0;
  source.R0 = (sqrt(2.0)/50.0)*r_e_nfw(pot_params);
  source.eta0 = 0.0;
  source.theta0 = 0.0;

  //_r_e = 1.1172980994467465;

  
// Added by Habib Dumet-Montoya
  fprintf(outls," PNFW Model Parameters: \n");  
  fprintf(outls," Scale Radius = %f\n",pot_params[0]);
  fprintf(outls," Characteristic Convergence = %f\n",pot_params[1]);
  fprintf(outls," Ellipticity of the Lensing Potential =  %f\n\n",pot_params[2]);
  fprintf(outls," Source Data: \n");  
  fprintf(outls," Source Position %f %f\n",source.x0,source.y0);
  fprintf(outls," Source Radius %f\n",source.R0);
  fprintf(outls," Source Ellipticity and Inclination %f %f\n\n",source.eta0,source.theta0);

  //print_arcs_nfw(source,npts);

  //for(int i=0;i<=1000;i++) printf("%E %E\n",0.9+0.5*double(i)/1000.0, lambda_t(0.9+0.5*double(i)/1000.0,pot_params));

  double *out = (double*)malloc(2.0*sizeof(double));
  bracketing_lambda_t( lambda_t_nfw, pot_params, out);
  printf("# Range of the root = %E  %E\n",out[0], out[1] );
  double z_out=  root_find(alpha_nfw, shear_nfw, pot_params,out[0] ,out[1]);
  fprintf(outls," Einstein Radius = %E\n",z_out);
  _r_e = z_out;
//   _r_e=r_e_nfw(pot_params);
  //   fprintf(outls, "Einstein Radius= %f\n",r_e_nfw(pot_params));
  theta_find(Df0Dtheta_nfw, source,pot_params );
  print_arcs_nfw(source,pot_params,npts);
}
