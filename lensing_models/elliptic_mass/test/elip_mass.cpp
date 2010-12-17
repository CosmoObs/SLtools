#include <cstdio>
#include <cmath>


#include "nfw_circular_model.h"
#include "integrals.h"
#include "substructure_model.h"

// g++ -Wall elip_mass.cpp `pkg-config gsl --cflags --libs`
int main(){
  double X1 = 0.1;
  double X2 = 0.0;
  double u = 12345;
  double conv_params[] = {0.7,0.9};
  double q = 0.9;
  double a = 1.0/sqrt(q);
  double b = sqrt(q);

  double out;
  for (int i=0;i<100;i++){
    out = J(2,conv_nfw_circ,X1,X2, conv_params, a, b);
    printf("%E %E\n",X1,out);
    out = K(1,conv_nfw_circ,X1,X2, conv_params, a, b);
    printf("%E %E\n",X1,out);
    X1 += 1./100.0;
  }
  return 0;
}
