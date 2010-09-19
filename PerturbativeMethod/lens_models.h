#ifndef LENSMODELS_H
#define LENSMODELS_H


double f1_SIS(double theta, double pot_params[]){
  return 0.0;
}

double Df0Dtheta_SIS(double theta, double pot_params[]){
  return 0.0;
}



//pot_params[0] = mp
//pot_params[1] = rp
//pot_params[2] = thetap
double f1_pert_SIS(double theta, double pot_params[]){
  double tmp1 = pot_params[0] * (1.0 - pot_params[1]*cos(theta-pot_params[2]) );
  double tmp2 = sqrt ( 1.0 - 2.0*pot_params[1]*cos(theta-pot_params[2]) + pot_params[1]*pot_params[1]);

  return tmp1/tmp2 ;
}

double Df0Dtheta_pert_SIS(double theta, double pot_params[]){
  double tmp1 = pot_params[0] * (pot_params[1]*sin(theta-pot_params[2]) );
  double tmp2 = sqrt ( 1.0 - 2.0*pot_params[1]*cos(theta-pot_params[2]) + pot_params[1]*pot_params[1]);
  return tmp1/tmp2;
}




#endif
