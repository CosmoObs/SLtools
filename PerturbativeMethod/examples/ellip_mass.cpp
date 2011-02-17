#include "elliptical_mass.h"
#include "enfw_model.h"

#include "pnfw_model2.h"
#include "siep_model.h"
#include "generate_curves.h"
#include "generate_arcs.h"
//#include "curvas.h"
// to compile make: g++ -Wall -o ellip ellip_mass.cpp `pkg-config gsl --cflags --libs`

int main(){

//verifica as funções Jn e Kn
/*  double eps = 0.;
  double ks = 1.0;
  double rs = 1.0;
  double pert_params[] = {eps,1.0,ks, rs};
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
  fclose(file5);*/



//verifica funcoes "f"s

/*  double eps = 0.0;
  double ks = 1.;
  double rs = 1.;
  double pert_params[] = {eps,2,ks, rs};
  double pert_params2[] = {ks, rs};
  double r_e = r_e_nfw_find(pert_params2);

  double a = 1.0/sqrt(1.0-eps);
  double b = sqrt(1.0-eps);
  double theta = 0.1;

  printf("Modelos Elipticos:\n");
  printf("  usando funções finais:\n");
  printf("             f1 = %E\n",nfw_ellip_f1(theta, pert_params, r_e));
  printf("      Df0Dtheta = %E\n",nfw_ellip_Df0Dtheta(theta, pert_params, r_e));
  printf("    D2f0Dtheta2 = %E\n",nfw_ellip_D2f0Dtheta2(theta, pert_params, r_e));

  printf("  usando funções gerais:\n");
  printf("             f1 = %E\n",f_1_ellip_mass(conv_nfw_circ, theta, pert_params2, r_e, a, b));
  printf("      Df0Dtheta = %E\n",Df0Dtheta_ellip_mass(conv_nfw_circ, theta, pert_params2, r_e, a, b));
  printf("    D2f0Dtheta2 = %E\n",D2f0Dtheta2_ellip_mass(conv_nfw_circ, theta, pert_params2, r_e, a, b));

  printf("\n\nModelos Pseudo-Elipticos: PNFW\n");
  printf("  usando funções finais:\n");
  printf("             f1 = %E\n",f1_pnfw(theta, pert_params, r_e));
  printf("      Df0Dtheta = %E\n",Df0Dtheta_pnfw(theta, pert_params, r_e));
  printf("    D2f0Dtheta2 = %E\n",D2f0Dtheta2_pnfw(theta, pert_params, r_e));

  printf("\n\nModelos Pseudo-Elipticos: PSIS\n");
  printf("  usando funções finais:\n");
  printf("             f1 = %E\n",f1_psis(theta, pert_params, r_e));
  printf("      Df0Dtheta = %E\n",Df0Dtheta_psis(theta, pert_params, r_e));
  printf("    D2f0Dtheta2 = %E\n",D2f0Dtheta2_psis(theta, pert_params, r_e));
  printf("\n\n\n");


  printf("Plotanto CC e CA\n");
  r_e= r_e_nfw_find(pert_params2);
  double kappa_2 = kappa2_nfw(pert_params2,r_e);
  printf("Einstein Radius = %f\n",r_e);
  printf("kappa_2 = %f\n",kappa_2);

  elliptical_source source;
  source.x0 = (r_e/10.);
  source.y0 = source.x0;
  printf(" Source Position %f %f\n",source.x0,source.y0);
  source.R0 = (sqrt(2.0)/25.0)*r_e;
  source.eta0 = 0.;
  source.theta0 = 2*3.14159/10.0;
  
  FILE *cc = fopen ("tang_crit.dat" , "w");
  FILE *ca = fopen ("tang_caust.dat" , "w");
  FILE *arcs=fopen ("arcs.dat", "w");
  FILE *src =fopen ("src.dat","w");
  plot_curves(nfw_ellip_f1, nfw_ellip_Df0Dtheta, nfw_ellip_D2f0Dtheta2, pert_params, kappa_2, r_e, 100, cc, ca);
  plot_arcs_sep(source,nfw_ellip_f1,nfw_ellip_Df0Dtheta, pert_params, kappa_2, r_e,100, arcs);
  plot_sources(source, 100,src);

  FILE *ccp = fopen ("tang_critp.dat" , "w");
  FILE *cap = fopen ("tang_caustp.dat" , "w");
  FILE *arcsp=fopen ("arcsp.dat", "w");
  FILE *srcp =fopen ("srcp.dat","w");
  plot_curves(f1_pnfw, Df0Dtheta_pnfw, D2f0Dtheta2_pnfw, pert_params, kappa_2, r_e, 100, ccp, cap);
  plot_arcs_sep(source,f1_pnfw,Df0Dtheta_pnfw, pert_params, kappa_2, r_e,100, arcsp);
  plot_sources(source, 100,srcp);*/



//graficos de CC e CA exatos (numericos)

  double q=1.0;
  double ks = 1.0;
  int I=1;
  int Npt = 100;


 // CAF(q, ks, I, Npt, "CA1.txt", "CA2.txt");
//  CCL(q, ks, I, Npt, "CC1.txt", "CC2.txt");





  return 0;
}
