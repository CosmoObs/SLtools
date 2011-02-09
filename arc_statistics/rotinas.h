#ifndef ROTINAS_H
#define ROTINAS_H

#include <cstring>
#include <cstdlib>

#include "integration.h"
#include <gsl/gsl_spline.h>



//return the inverse of the adimensional Hublle paremeter
//for LCDM: cosmological_params[0] -> Mass density
//          cosmological_params[1] -> 'Curvature density'
//          cosmological_params[2] -> 'Cosmological constant density'
double E_inv(double z, double cosmological_params[]){//(parametro de Hubble adimensional)^(-1)

  double zm1 = z + 1.0;
  double zm1_2 = zm1*zm1;
  double zm1_3 = zm1_2*zm1;
  //printf("%E\t%E\t%E\n", cosmological_params[0], cosmological_params[1], cosmological_params[2]);

  double out = 1.0/sqrt(cosmological_params[0]*zm1_3 + cosmological_params[1]*zm1_2 + cosmological_params[2]);

  return out;
}

inline double Int_E_inv(double z_init, double z_final, double cosmological_params[]){ return IntGaulegSub50(E_inv, cosmological_params, z_init, z_final);}

double I(double z_init, double z_final, double cosmological_params[]){//integral usada para calculos de distancias cosmologicas

  double out = 0.0;

  if(cosmological_params[1] > 1E-6){
    double SQRTA1 = sqrt(fabs(cosmological_params[1]));

    out =  sinh( SQRTA1 * Int_E_inv(z_init, z_final, cosmological_params) )/SQRTA1;
  }
  else if(cosmological_params[1] < -1E-6){
    double SQRTA1 = sqrt(fabs(cosmological_params[1]));

    out =  sin( SQRTA1 * Int_E_inv(z_init, z_final, cosmological_params) )/SQRTA1;
  }
  else if(fabs(cosmological_params[1]) <= 1E-6){
    out = Int_E_inv(z_init, z_final, cosmological_params) ;
  }
  else printf("ERROR\n");

  return out;
}

//return the concentration parameter
double conc(double z_l, double M, double cosmological_params[], long int seed = 0){//compute de concentrarion parameter
  //double c_0 = 5.97;
  //double alc = -0.102;
  //double M_st = 8.0;
  //return c_0/(1.0+z)*pow(M/M_st,alc);
  //printf("%E\n",4.67*pow(M,-0.11)*pow(invE(z,Om0,Ol0),0.666666667));
  //printf("%E\t%E\n",4.67*pow(M,-0.11)*pow(invE(z,Om0,Ol0),0.666666667),4.67*pow(M,-0.11)/(1.0+z));
  //return 4.67*pow(M,-0.11)*pow(invE(z,Om0,Ol0),0.666666667);
  //neto, E

  //return 4.67*pow(M,-0.11)/(1.0+z);
  //neto, (1+z)

  //return 4.3*pow(M,-0.098)*pow(E_inv(z_l,cosmological_params),0.666666667);
  //maccio, WMAP3 com E;***

  //return 4.3*pow(M,-0.098)/(1.0+z);
  //maccio, WMAP3 com 1/(1+z);


  //return 9.59/(1.0 + z)*pow(M,-0.102);
  //dolag

  //return ( 7.03289888/(1.0 + z_l) )*pow(M,-0.13);
  //Bullock 2001

  //return (7.03289888)*pow(M,-0.13)*pow(invE(z,Om0,Ol0),0.666666667);
  //Bullock 2001 com E(z)!!!!!

  //return 3.0;


  //Utiliza os valores dos best-fits de Gao-2008
  M = M*1e14;
  double x[5], y[5];
  double A[5];
  double B[5];

  A[0] = -0.138;
  A[1] = -0.125;
  A[2] = -0.092;
  A[3] = -0.031;
  A[4] = -0.004;

  B[0] = 2.646;
  B[1] = 2.372;
  B[2] = 1.891;
  B[3] = 0.985;
  B[4] = 0.577;

  for(int i=0;i<5;i++) y[i] = pow(10.0,B[i])*pow(M,A[i]);

  x[0] = 0.0;
  x[1] = 0.5;
  x[2] = 1.0;
  x[3] = 2.0;
  x[4] = 3.0;


/*  printf("%E %E\n",x[0], y[0] = pow(10.0,2.646)*pow(M,-0.138));
  printf("%E %E\n",x[1], y[1] = pow(10.0,2.372)*pow(M,-0.125));
  printf("%E %E\n",x[2], y[2] = pow(10.0,1.891)*pow(M,-0.092));
  printf("%E %E\n",x[3], y[3] = pow(10.0,0.985)*pow(M,-0.031));
  printf("%E %E\n",x[4], y[4] = pow(10.0,0.577)*pow(M,-0.004));*/

  gsl_interp_accel *acc = gsl_interp_accel_alloc ();
  gsl_spline *spline = gsl_spline_alloc (gsl_interp_polynomial, 5);
  gsl_spline_init (spline, x, y, 5);

  return gsl_spline_eval (spline, z_l, acc);

  //chama o programa do angelo para calcular conc
/*  M = M*1e14;
  char str [120];
  char str2 [120];

  FILE * out1  = fopen ("myfile.txt" , "w+");
  fprintf(out1,"./cmz -m %e -z %e",M,z_l);
  rewind(out1);
  fgets (str, 120 , out1);
  fclose (out1);

  FILE * out11  = fopen ("myfile.txt" , "w+");
  fprintf(out11,"%li",seed);
  rewind(out11);
  fgets (str2, 120 , out1);

  strcat(str," -s ");
  strcat(str,str2);
  strcat(str," > cmz-out.txt");
  system(str);
  FILE * out2  = fopen ("cmz-out.txt" , "r");
  fgets(str2, 120, out2);
  //puts(str);
  //puts(str2);
  fclose(out11);
  fclose(out2);
  return atof(str2);*/
}

double g(double conc){
  double UMC = 1.0 + conc;
  return conc*conc/( log(UMC) - conc/UMC );
}

double DELTA(double z, double cosmological_params[]){
  return 200;
}

double F(double z, double cosmological_params[]){
  //double zm1 =z+1.0;
  //double zm12=zm1*zm1;
 // double zm13=zm1*zm12;
  //double Ok0 = 1.0 - Om0 - Ol0;

 // double Oz = Om0*pow(1.0+z,3.0)*pow(invE(z,Om0,Ol0),2.0);


  return pow(E_inv(z,cosmological_params),-2.0);
  //para rho_ref = rho_crit , delta=200 e delta=deltaC

  //return Om0*pow(1.0+z,3.0);
  //para rho_ref = rho_m , delta=200 e delta=deltaC

  //return Om0*pow(1.0+z,3.0) + Ol0;
  //para rho_ref = rho_u , delta=200 e delta=deltaC


  //return pow(invE(z,Om0,Ol0),-2.0)/Oz;
  //para rho_ref = rho_crit , elta=deltaVir

  //return Om0*pow(1.0+z,3.0)/Oz;
  //para rho_ref = rho_m , delta=deltaVir

  //return (Om0*pow(1.0+z,3.0) + Ol0)/Oz;
  //para rho_ref = rho_u , delta=deltaVir

}



inline double kappa_s(double zl, double zs, double M, double conc, double cosmological_params[]){return 0.0007362440639914055*(g(conc)*pow(M,0.3333333333333333)*pow(F(zl,cosmological_params)*DELTA(zl,cosmological_params),0.6666666666666666))/(1.0 + zl)*I(zl,zs,cosmological_params)*I(0.0,zl,cosmological_params)/I(0.0,zs,cosmological_params);}

double x_s(double zl, double M, double conc, double cosmological_params[]){ return 5.062041306*pow(M,0.3333333333333333)*pow(F(zl,cosmological_params)*DELTA(zl,cosmological_params),-0.3333333333333333)*(1.0+zl)/(I(0.0,zl,cosmological_params)*conc); }
//x_s em MINUTOS DE ARCO

double n(double z_source, double m_lim){//em MINUTOS DE ARCO
  double zm = 0.23*(m_lim-20.6);
  double N = 45.0*pow(zm,3.4);
  double zstar = zm/1.412;
  double zstar3 = zstar*zstar*zstar;

  double out = N * z_source*z_source/zstar3*exp(-pow(z_source/zstar,1.5));
  return out;
}

int DataFile(FILE *arquivo,int coluna1, int coluna2, int NumeroLinhas, double *saida1, double *saida2)
{

  char IN[500];
  double out[20];

  printf("Lendo dados do arquivo.\n");

  for(int i=1;i<=NumeroLinhas;i++)
  {
    fgets(IN,500,arquivo);
    sscanf(IN,"%lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg",
               &out[0] ,&out[1] ,&out[2] ,&out[3] ,&out[4], &out[5], &out[6], &out[7], &out[8], &out[9],
               &out[10],&out[11],&out[12],&out[13],&out[14],&out[15],&out[16],&out[17],&out[18],&out[19],&out[20]);

    saida1[i-1] = out[coluna1-1];
    saida2[i-1] = out[coluna2-1];
  }

  printf("->Leitura finalizada.\n");
  return 0;
}

int DataFile2Coluns(FILE *arquivo, int NumeroLinhas, double *saida1, double *saida2)
{

  char IN[500];
  double out[20];

  printf("Lendo dados do arquivo.\n");

  for(int i=1;i<=NumeroLinhas;i++)
  {
    fgets(IN,500,arquivo);
    sscanf(IN,"%lg %lg", &out[0] ,&out[1]);

    saida1[i-1] = out[0];
    saida2[i-1] = out[1];
  }

  printf("->Leitura finalizada.\n");
  return 0;
}

int DataFile3Coluns(FILE *arquivo, int NumeroLinhas, double *saida1, double *saida2, double *saida3)
{

  char IN[500];
  double out[20];

  printf("Lendo dados do arquivo.\n");

  for(int i=1;i<=NumeroLinhas;i++)
  {
    fgets(IN,500,arquivo);
    sscanf(IN,"%lg %lg %lg", &out[0] ,&out[1], &out[2]);

    saida1[i-1] = out[0];
    saida2[i-1] = out[1];
    saida3[i-1] = out[2];
  }

  printf("->Leitura finalizada.\n");
  return 0;
}

int DataFile2(FILE *arquivo, int LinhaLida, int NumeroLinhas, double **saida)
{

  char *IN = (char*) malloc(1000*sizeof(char));

  double *out = (double*) malloc(41*sizeof(double));


  printf("Lendo dados do arquivo.\n");

  for(int i=1;i<=NumeroLinhas;i++)
  {
    fgets(IN,1000,arquivo);
    sscanf(IN,"%lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg",
               &out[0] ,&out[1] ,&out[2] ,&out[3] ,&out[4], &out[5], &out[6], &out[7], &out[8], &out[9],
               &out[10],&out[11],&out[12],&out[13],&out[14],&out[15],&out[16],&out[17],&out[18],&out[19],
               &out[20],&out[21],&out[22],&out[23],&out[24],&out[25],&out[26],&out[27],&out[28],&out[29],
               &out[30],&out[31],&out[32],&out[33],&out[34],&out[35],&out[36],&out[37],&out[38],&out[39],
               &out[40],&out[41]);

    for(int j=0;j<42;j++) saida[i-1][j] = out[j];
   // printf("%i\n",i);
  }

  printf("->Leitura finalizada.\n");
  return 0;
}


int DataFile21(FILE *arquivo, int LinhaLida, int NumeroLinhas, double **saida)
{

  char *IN = (char*) malloc(1000*sizeof(char));

  double *out = (double*) malloc(3*sizeof(double));


  printf("Lendo dados do arquivo.\n");

  for(int i=1;i<=NumeroLinhas;i++)
  {
    fgets(IN,1000,arquivo);
    sscanf(IN,"%lg %lg %lg",&out[0] ,&out[1] ,&out[2]);

    for(int j=0;j<3;j++) saida[i-1][j] = out[j];
   // printf("%i\n",i);
  }

  printf("->Leitura finalizada.\n");
  return 0;
}

void fim(){
    std::cout << "tempo decorrido: " << double(clock())/(double(CLOCKS_PER_SEC)*60.0) << " minutos" << std::endl;
    std::cout << "FIM!" << std::endl;
}

#endif

