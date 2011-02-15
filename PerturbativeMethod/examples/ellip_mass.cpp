#include "elliptical_mass.h"
#include "enfw_model.h"


// to compile make: g++ -Wall -o ellip ellip_mass.cpp `pkg-config gsl --cflags --libs`

int main(){

  double eps = 0.5;
  double ks = 1.1;
  double rs = 1.0;
  double pert_params[] = {ks, rs};
  double r_e = 1.3;

  double a = 1.0/sqrt(1.0-eps);
  double b = sqrt(1.0-eps);
  double theta = 0.0;

  double X1 = 0.0;
  double X2 = 0.0;

  FILE *file1  = fopen("out1.txt","w+");
  FILE *file2  = fopen("out2.txt","w+");
  FILE *file3  = fopen("out3.txt","w+");
  FILE *file4  = fopen("out4.txt","w+");
  FILE *file5  = fopen("out5.txt","w+");
  for(int i=0;i<=100;i++){
    X1 = r_e*cos(theta);
    X2 = r_e*sin(theta);
    fprintf(file1,"%E %E\n",theta, J(0, conv_nfw_circ, X1, X2, pert_params, a, b));
    fprintf(file2,"%E %E\n",theta, J(1, conv_nfw_circ, X1, X2, pert_params, a, b));

    fprintf(file3,"%E %E\n",theta, K(0, conv_nfw_circ_prime, X1, X2, pert_params, a, b));
    fprintf(file4,"%E %E\n",theta, K(1, conv_nfw_circ_prime, X1, X2, pert_params, a, b));
    fprintf(file5,"%E %E\n",theta, K(2, conv_nfw_circ_prime, X1, X2, pert_params, a, b));
    theta += 2.0*3.14159/100.0;
  }

  fclose(file1);
  fclose(file2);
  fclose(file3);
  fclose(file4);
  fclose(file5);
//verifica funcoes "f"s
/*  printf("%E\n",nfw_ellip_f1(theta, pert_params, r_e));

  printf("%E\n",f_1_ellip_mass(conv_nfw_circ, theta, pert_params, r_e, a, b));
  printf("%E\n",Df0Dtheta_ellip_mass(conv_nfw_circ, theta, pert_params, r_e, a, b));
  printf("%E\n",D2f0Dtheta2_ellip_mass(conv_nfw_circ, theta, pert_params, r_e, a, b));
*/

  return 0;
}
