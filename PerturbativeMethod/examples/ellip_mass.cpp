#include "elliptical_mass.h"
#include "enfw_model.h"


// to compile make: g++ -Wall -o ellip ellip_mass.cpp `pkg-config gsl --cflags --libs`

int main(){
  double theta =0.0;
  double pert_params[] = {0.1,1, 1.0, 1.0};
  double r_e = 1.0;

  printf("%E\n",nfw_ellip_f1(theta, pert_params, r_e));

  return 0;
}
