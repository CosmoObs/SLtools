#include <cmath>
#include <iostream>
#include <cstdlib>


#include "classes.h"

double W=-1.0; //parametro da equacao de estado p=w\rho

double invE(double z, double Om0 = 0.3, double Ol0 = 0.7){//inverse of function E;E=H/H_0
  double zm1 =z+1.0;
  double zm12=zm1*zm1;
  double zm13=zm1*zm12;
  double Ok0 = 1.0 - Om0 - Ol0;
  return 1.0/sqrt( Om0*zm13 + Ok0*zm12 + Ol0*pow(zm1,3.0*(1.0+W)) );
}

double F(double z, double Om0 = 0.3, double Ol0 = 0.7){
  //double zm1 =z+1.0;
  //double zm12=zm1*zm1;
  //double zm13=zm1*zm12;
  //double Ok0 = 1.0 - Om0 - Ol0;

  //double Oz = Om0*zm13*pow(invE(z,Om0,Ol0),2.0);


  return pow(invE(z,Om0,Ol0),-2.0);
  //para rho_ref = rho_crit , delta=200 e delta=deltaC
}



double g(double C){
  return C*C/( log(1.0+C) - C/(1.0+C) );
}

double DELTA(double z, double Om0 = 0.3, double Ol0 = 0.7){//DELTA_200!!!!!
  return 200.0;
}

double Kappa_s_incomplete(double zl, double M, double C, double Om0 = 0.3, double Ol0 = 0.7){

  return 0.0007362440639914055*
         (g(C)*pow(M,0.3333333333333333)*pow(F(zl,Om0,Ol0)*DELTA(zl,Om0,Ol0),0.6666666666666666))/(1.0 + zl);
  //em funcao de F
}

double S_k(double zl, double zs, double Om0 = 0.3, double Ol0 = 0.7){
  Integral out(invE,zl,zs,20);
  double Ok = 1.0-Om0-Ol0;
  double sqrtOk = sqrt(fabs(Ok));

  if(Ok > 10E-6) return sinh(sqrtOk*out.calcIntegral(Om0,Ol0,1))/sqrtOk;

  if(Ok < -10E-6) return sin(sqrtOk*out.calcIntegral(Om0,Ol0,1))/sqrtOk;

  if(fabs(Ok)<=10E-6) return out.calcIntegral(Om0,Ol0,1);

  printf("ERROR in X\n");
  return 0.0;
}


double I(double zl, double zs, double Om0 = 0.3, double Ol0 = 0.7){return S_k(zl,zs,Om0,Ol0);}

double kappa_s(double zl, double zs, double M, double C, double Om0, double Ol0){
  return Kappa_s_incomplete(zl,M,C,Om0,Ol0)*I(zl,zs,Om0,Ol0)*I(0.0,zl,Om0,Ol0)/I(0.0,zs,Om0,Ol0);
}

double x_s(double zl, double M, double C, double Om0, double Ol0){
  return 5.062041306*
      pow(M,0.3333333333333333)*pow(F(zl,Om0,Ol0)*DELTA(zl,Om0,Ol0),-0.3333333333333333)*(1.0+zl)/(I(0.0,zl,Om0,Ol0)*C);
   //em funcao de F, em MINUTOS DE ARCO
}




int main(int argc, char **argv){

  double Zl, Zs, M, C, Om0, Ol0;

  if(argc < 7){
    printf("ERRO em %s-> menos de 6 argumentos de entrada\n",argv[0]);
    FILE *out = fopen("out.txt","w+");
    fprintf(out,"ERRO em %s-> menos de 6 argumentos de entrada\n",argv[0]);
    fclose(out);
    return 0;
   }

  Zl = atof(argv[1]);
  Zs = atof(argv[2]);
  M  = atof(argv[3]);
  C  = atof(argv[4]);
  Om0= atof(argv[5]);
  Ol0= atof(argv[6]);

//  printf("Zl=%f,Zs=%f,M=%f,C=%f,Om0=%f,Ol0=%f\n",Zl,Zs,M,C,Om0,Ol0);

//  FILE *out = fopen("out.txt","w+");

//  fprintf(out,"Zl=%f,Zs=%f,M=%f,C=%f,Om0=%f,Ol0=%f\n",Zl,Zs,M,C,Om0,Ol0);
//  fprintf(out,"%E\t%E\n",kappa_s(Zl,Zs,M,C,Om0,Ol0),x_s(Zl,M,C,Om0,Ol0));
  printf("kappa_s= %f\tx_s= %f\n",kappa_s(Zl,Zs,M,C,Om0,Ol0),x_s(Zl,M,C,Om0,Ol0));

//  fclose(out);

/*  printf("Entre com Zl, Zs, M (em unidades de 10^14M_solh^(-1)), C, Om0, Ol0\n");
  printf("Zl =");scanf("%lf",&Zl );
  printf("Zs =");scanf("%lf",&Zs );
  printf("M  =");scanf("%lf",&M  );
  printf("C  =");scanf("%lf",&C  );
  printf("Om0=");scanf("%lf",&Om0);
  printf("Ol0=");scanf("%lf",&Ol0);
*/
  //for(int i=0;i<10000;i++) kappa_s(1,5,3,3,0.3,0.7);

//  printf("kappa_s=%E\n",kappa_s(Zl,Zs,M,C,Om0,Ol0));
//  printf("x_s    =%E Minutos de Arco\n",x_s(Zl,M,C,Om0,Ol0));


  return 0;
}
