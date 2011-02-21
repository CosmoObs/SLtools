#include "elliptical_mass.h"
#include "enfw_model.h"
#include "sie_model.h"


#include "sis_model.h"

#include "pnfw_model2.h"
#include "siep_model.h"
#include "generate_curves.h"
#include "generate_arcs.h"


#include "lente.h"
#include "integrais.h"
#include "propriedades.h"
#include "curvas.h"
#include "area.h"


// to compile make: g++ -Wall -o ellip ellip_mass.cpp `pkg-config gsl --cflags --libs`

int main(){

//verifica as funções Jn e Kn para NFW
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


//verifica as funções Jn e Kn para SIS
  double eps = 0.5;
  double b_sis = 1.0;
  //double pert_params[] = {eps,2.0,b_sis};
  double pert_params[] = {b_sis};
  double r_e = b_sis;

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
    fprintf(file1,"%E %E\n",theta, J(0, conv_sis_circ, X1, X2, pert_params, a, b) - r_e/b_sis*J(0, conv_sis_circ, X1, X2, pert_params, 1, 1)+1);
    fprintf(file2,"%E %E\n",theta, J(1, conv_sis_circ, X1, X2, pert_params, a, b) - r_e/b_sis*J(0, conv_sis_circ, X1, X2, pert_params, 1, 1)+1);

    fprintf(file3,"%E %E\n",theta, r_e*r_e*K(0, conv_sis_circ_prime, X1, X2, pert_params, a, b));
    fprintf(file4,"%E %E\n",theta, r_e*r_e*K(1, conv_sis_circ_prime, X1, X2, pert_params, a, b));
    fprintf(file5,"%E %E\n",theta, r_e*r_e*K(2, conv_sis_circ_prime, X1, X2, pert_params, a, b));
    theta += 2.0*3.14159/100.0;
  }
    fprintf(file1,"\n\n");
    fprintf(file2,"\n\n");
    fprintf(file3,"\n\n");
    fprintf(file4,"\n\n");
    fprintf(file5,"\n\n");
  theta = 0.01;
  for(int i=0;i<=100;i++){
    fprintf(file1,"%E %E\n",theta, J0_sie_analytic(r_e, theta, b_sis, b, a));
    fprintf(file2,"%E %E\n",theta, J1_sie_analytic(r_e, theta, b_sis, b, a));

    fprintf(file3,"%E %E\n",theta, K0_sie_analytic(r_e, theta, b_sis, a, b));
    fprintf(file4,"%E %E\n",theta, K1_sie_analytic(r_e, theta, b_sis, a, b));
    fprintf(file5,"%E %E\n",theta, K2_sie_analytic(r_e, theta, b_sis, a, b));
    theta += 2.0*3.14159/100.0;
  }




  fclose(file1);
  fclose(file2);
  fclose(file3);
  fclose(file4);
  fclose(file5);

/********************************************************************************************************************/
/********************************************************************************************************************/
//verifica funcoes "f"s para NFW

/*  double eps = 0.5;
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
  printf("             f1 = %E\n",f1_enfw(theta, pert_params, r_e));
  printf("      Df0Dtheta = %E\n",Df0Dtheta_enfw(theta, pert_params, r_e));
  printf("    D2f0Dtheta2 = %E\n",D2f0Dtheta2_enfw(theta, pert_params, r_e));

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
  plot_curves(f1_enfw, Df0Dtheta_enfw, D2f0Dtheta2_enfw, pert_params, kappa_2, r_e, 500, cc, ca);
  plot_arcs_sep(source,f1_enfw,Df0Dtheta_enfw, pert_params, kappa_2, r_e,500, arcs);
  plot_sources(source, 100,src);

  FILE *ccp = fopen ("tang_critp.dat" , "w");
  FILE *cap = fopen ("tang_caustp.dat" , "w");
  FILE *arcsp=fopen ("arcsp.dat", "w");
  FILE *srcp =fopen ("srcp.dat","w");
  plot_curves(f1_pnfw, Df0Dtheta_pnfw, D2f0Dtheta2_pnfw, pert_params, kappa_2, r_e, 500, ccp, cap);
  plot_arcs_sep(source,f1_pnfw,Df0Dtheta_pnfw, pert_params, kappa_2, r_e,500, arcsp);
  plot_sources(source, 100,srcp);



  //graficos de CC e CA exatos (numericos)

  double q=1.0 - eps;
 // double ks = 1.0;
  int I=1;
  int Npt = 500;


  CAF(q, ks, I, Npt, "CA1.txt", "CA2.txt");
  CCL(q, ks, I, Npt, "CC1.txt", "CC2.txt");
*/



/********************************************************************************************************************/
/********************************************************************************************************************/
//verifica funcoes "f"s para SIS

/*  double eps = 0.5;
  double b_sis = 1.9;
  double pert_params[] = {eps,2,b_sis};
  //double pert_params2[] = {b_sis};
  double r_e = b_sis;

  //double a = 1.0/sqrt(1.0-eps);
  //double b = sqrt(1.0-eps);
  double theta = 0.0001;

  FILE *file1  = fopen("out1.txt","w+");
  FILE *file2  = fopen("out2.txt","w+");
  FILE *file3  = fopen("out3.txt","w+");
  FILE *file4  = fopen("out4.txt","w+");
  for(int i=0;i<200;i++){
    fprintf(file1,"%E %E\n",theta,f1_sie_analytic(theta, pert_params, r_e));
    fprintf(file2,"%E %E\n",theta,Df0Dtheta_sie_analytic(theta, pert_params, r_e));
    fprintf(file3,"%E %E\n",theta,D2f0Dtheta2_sie_analytic(theta, pert_params, r_e));

    fprintf(file4,"%E %E\n",theta,f1_sie_analytic(theta, pert_params, r_e) - f1_sie(theta, pert_params, r_e));
    theta+=2*3.14159/200.0;
  }
  fprintf(file1,"\n\n");
  fprintf(file2,"\n\n");
  fprintf(file3,"\n\n");
  fprintf(file4,"\n\n");
  theta = 0.0001;
  for(int i=0;i<200;i++){
    fprintf(file1,"%E %E\n",theta,f1_sie(theta, pert_params, r_e));
    fprintf(file2,"%E %E\n",theta,Df0Dtheta_sie(theta, pert_params, r_e));
    fprintf(file3,"%E %E\n",theta,D2f0Dtheta2_sie(theta, pert_params, r_e));

    fprintf(file4,"%E %E\n",theta,Df0Dtheta_sie_analytic(theta, pert_params, r_e) - Df0Dtheta_sie(theta, pert_params, r_e));
    theta+=2*3.14159/200.0;
  }
  fprintf(file4,"\n\n");
  theta = 0.0001;
  for(int i=0;i<200;i++){
    fprintf(file4,"%E %E\n",theta,D2f0Dtheta2_sie_analytic(theta, pert_params, r_e) - D2f0Dtheta2_sie(theta, pert_params, r_e));
    theta+=2*3.14159/200.0;
  }

  fclose(file1);
  fclose(file2);
  fclose(file3);
  fclose(file4);*/

/*  printf("Plotanto CC e CA\n");
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
  plot_curves(f1_enfw, Df0Dtheta_enfw, D2f0Dtheta2_enfw, pert_params, kappa_2, r_e, 500, cc, ca);
  plot_arcs_sep(source,f1_enfw,Df0Dtheta_enfw, pert_params, kappa_2, r_e,500, arcs);
  plot_sources(source, 100,src);

  FILE *ccp = fopen ("tang_critp.dat" , "w");
  FILE *cap = fopen ("tang_caustp.dat" , "w");
  FILE *arcsp=fopen ("arcsp.dat", "w");
  FILE *srcp =fopen ("srcp.dat","w");
  plot_curves(f1_pnfw, Df0Dtheta_pnfw, D2f0Dtheta2_pnfw, pert_params, kappa_2, r_e, 500, ccp, cap);
  plot_arcs_sep(source,f1_pnfw,Df0Dtheta_pnfw, pert_params, kappa_2, r_e,500, arcsp);
  plot_sources(source, 100,srcp);



  //graficos de CC e CA exatos (numericos)

  double q=1.0 - eps;
 // double ks = 1.0;
  int I=1;
  int Npt = 500;


  CAF(q, ks, I, Npt, "CA1.txt", "CA2.txt");
  CCL(q, ks, I, Npt, "CC1.txt", "CC2.txt");
*/
  return 0;
}
