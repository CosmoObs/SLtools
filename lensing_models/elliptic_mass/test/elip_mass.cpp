#include <cstdio>
#include <cmath>


#include "nfw_circular_model.h"
#include "integrals.h"

// g++ -Wall elip_mass.cpp `pkg-config gsl --cflags --libs`
int main(){
  double X1 = 0.01;
  double X2 = 0.02;
  double u = 1E-10;
  double conv_params[] = {0.7,1.0};
  double q = 0.9;
  double a = 1.0/sqrt(q);
  double b = sqrt(q);

  double out;
  for (int i=0;i<=100000;i++){
    out = kCSI2U(conv_nfw_circ,X1,X2, conv_params, u, a, b);
    printf("%E %E\n",u,out);
    u += 0.01/100000.0;
  }
  return 0;
}
