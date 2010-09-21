#include <cstdlib>
#include <cstdio>

#include "../perturbative_method.h"
#include "../pnfw_model.h"
#include "../theta_find.h"

void print_arcs_nfw(elliptical_source source_in, int npts=200){

  double theta = 0.0;
  double r_p = 0.0;
  double r_m = 0.0;

//pot_paramas[0]=r_s, pot_params[1]=ks,pot_params[2]=elipticity,
  double pot_params[] = {0.9,1.0,0.2};

    for(int i=0;i<=npts;i++){

      if(arg_sqrt(Df0Dtheta_nfw, source_in, theta, pot_params)>0.0){
      r_p = _r_e+dr_plus(f1_nfw, Df0Dtheta_nfw, kappa2_nfw(pot_params), source_in, theta, pot_params);
      r_m = _r_e+dr_minus(f1_nfw, Df0Dtheta_nfw,kappa2_nfw(pot_params), source_in, theta, pot_params);
      printf("%E %E\n",r_p*cos(theta),r_p*sin(theta));
      printf("%E %E\n",r_m*cos(theta),r_m*sin(theta));
    }
//     printf("%E %E\n",r_e*cos(theta),r_e*sin(theta));
    theta+= 6.283185308/npts;
  }
}

int main(){

  //int npts = 1000;

  elliptical_source source;

  source.x0 = 0.1;
  source.y0 = 0.0;
  source.R0 = 0.05;
  source.eta0 = 0.0;
  source.theta0 =  0.0;
  _r_e = 1.1172980994467465;

  double pot_params[] = {0.,0.5,0.2};

  //print_arcs_nfw(source,npts);
  theta_find(Df0Dtheta_nfw, source,pot_params );
  //for(int i=0;i<=1000;i++) printf("%E %E\n",0.9+0.5*double(i)/1000.0, lambda_t(0.9+0.5*double(i)/1000.0,pot_params));


  intervalo_teste( lambda_t, pot_params, 2.0);

  root_find(alpha_nfw, conv_nfw, pot_params, 0.552000, 0.552000);

}
