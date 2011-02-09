#include <cstdio>
#include <cmath>
#include <ctime>

#include <iostream>
#include <algorithm>
#include <gsl/gsl_statistics.h>


#define PI 3.1415926535897931159979634685441851615905761718750



#include "integration.h"

#include "rotinas.h"

#include "cross-section.h"

#include "integration-class.h"

#include "farcs.h"


//function to perfor tests
void test(){
  double zl, zs, Om, Ok, Ol, M, Conc, m_lim; 
  
  zl = 0.4;
  zs = 2.0;
  Om = 0.3;
  Ok = 0.0;
  Ol = 0.7;
  double cosmological_params[] = {Om, Ok, Ol};
  M = 5.0;
  Conc = conc(zl, M, cosmological_params, 0);
  m_lim = 24.0;

  printf("\n\n------------------------------Test program started------------------------------\n");
  printf("The imput parameters are:\n");
  printf("zl=%E, zs=%E, M=%E\n",zl,zs,M);
  printf("Om=%E, Ok=%E, Ol=%E\n",cosmological_params[0],cosmological_params[1],cosmological_params[2]);

  printf("m_lim=%E, ",m_lim);
  CrossSection A;
  A.cross_section_init(0.4);


  printf("\nSome Results:\n");
  double ks = kappa_s(zl, zs, M, Conc, cosmological_params);
  double xs = x_s(zl, M, Conc, cosmological_params);

  printf("NFW parameters-> kappa_s=%E, x_s=%E\n",ks, xs);
  printf("Interpolations-> area=%E, ",A.area(ks));
  printf("mag_med=%E, ",A.mag_med(ks));
  printf("mag_min=%E\n",A.mag_min(ks));

  double params[] = {zl, M, Conc, m_lim, cosmological_params[0], cosmological_params[1], cosmological_params[2]};
  printf("f_arc->  f_arcs_I_arg=%E,  f_arcs_I=%E\n", f_arcs_I_arg(zs, params , A), f_arcs_I(zl, M, Conc, m_lim, cosmological_params, A));
  printf("f_arc-> f_arcs_II_arg=%E, f_arcs_II=%E\n",f_arcs_II_arg(zs, params , A), f_arcs_II(zl, M, Conc, m_lim, cosmological_params, A));
  //printf("f_arc-> Int_farcs_III_arg_mu=%E\n", Int_f_arcs_III_arg_mu(zl, zs, M, Conc, m_lim, cosmological_params, A));
  //printf("f_arc-> farcs_III_arg=%E\n",f_arcs_III_arg(zs, params, A));
  printf("f_arc-> farcs_III=%E\n",f_arcs_III(zl, M, Conc, m_lim, cosmological_params, A));


 //um for para testar a distribuicao de f_arcs usando o programa cmz
 /* FILE * out  = fopen ("out.txt" , "w+");
  for(int i=1;i<=5000;i++){
    Conc = conc(zl, M, cosmological_params, 0);
    fprintf(out,"%E\n",f_arcs_II(zl, M, Conc, m_lim, cosmological_params, A));
  }
  fclose(out);*/

  //printf("f_arc-> f_arcs_III=%E\n",f_arcs_III(zl, M, Conc, m_lim, cosmological_params, A));

  printf("-----------------------------Test program finished!-----------------------------\n\n");
}

//compute f_arcs for the cluster in the file catalogoZ.txt an print the output in farcs-cmz-out.txt
void Conpute_Farcs_Cathalog(){
//Leitura do catalogo
  int nlinhas = 247697;

  double *Z = (double*) malloc(nlinhas*sizeof(double));
  double *M = (double*) malloc(nlinhas*sizeof(double));

  //int Zdiv=20, Mdiv=80; //numero de bins em Z e M

 // double  Zmin = 0.07000750303268433, Zmax = 1.4535800218582153;
 // double  Mmin = 49500001402880.0,    Mmax = 2146499966271488.0;

  FILE *FileIn  = fopen("catalogoZ.dat","r"); //abre o arquivo de entrada
  DataFile(FileIn, 2, 17 , nlinhas, M, Z); //enche os vetores M e Z
  fclose(FileIn);  //fecha o arquivo de entrada

  for (int i=0;i<=(nlinhas-1);i++) M[i] = M[i]/1E14; //Comverta a massa em unidades de 10^14
  //Mmin = Mmin/1E14;
  //Mmax = Mmax/1E14;

  //Testes para verificar os vetors M e Z
  printf("%.14f\t%f\n",Z[0],M[0]);
  printf("%.14f\t%f\n",Z[nlinhas-1],M[nlinhas-1]);
//Finaliza leitura do catalogo
fim();

  double zl = 0.4;
  double Om = 0.3;
  double Ok = 0.0;
  double Ol = 0.7;
  double cosmological_params[] = {Om, Ok, Ol};
  double Massa = 5.0;
  double Conc = conc(zl, Massa, cosmological_params, 0);
  double m_lim = 24.0;
  CrossSection A;
  A.cross_section_init(0.4);

  //for (int i=0;i<10000;i++) printf("farcs=%E\n",f_arcs_II(Z[i], M[i], Conc, m_lim, cosmological_params, A));

  FILE * out  = fopen ("farcs-cmz-out.txt" , "w+");

  for (int i=0;i<nlinhas;i++){
    fprintf(out,"%E\t%E\t",M[i],Z[i]);
    for (int j=1;j<=1;j++){
      Conc = conc(Z[i], M[i], cosmological_params);
      fprintf(out,"%E\t",f_arcs_II(Z[i], M[i], Conc, m_lim, cosmological_params, A));

    }
  fprintf(out,"\n");
  //fim();
  }


}

//compute f_arcs for the cluster in the file catalogoZ-DC5.txt an print the output in farcs-cmz-out-DC5.txt
void Conpute_Farcs_Cathalog_DC5(){

  int nlinhas = 3402;

  double *Z = (double*) malloc(nlinhas*sizeof(double));
  double *M = (double*) malloc(nlinhas*sizeof(double));

  FILE *FileIn  = fopen("catalogoZ-DC5.dat","r");
  DataFile2Coluns(FileIn, nlinhas, Z, M);
  //for (int i=0;i<nlinhas;i++)printf("%E %E\n",Z[i],M[i]);
  fclose(FileIn);  //fecha o arquivo de entrada

  for (int i=0;i<=(nlinhas-1);i++) M[i] = M[i]/1E14; //Comverta a massa em unidades de 10^14

  double zl = 0.4;
  double Om = 0.3;
  double Ok = 0.0;
  double Ol = 0.7;
  double cosmological_params[] = {Om, Ok, Ol};
  double Massa = 5.0;
  double Conc = conc(zl, Massa, cosmological_params, 0);
  double m_lim = 24.0;
  CrossSection A;
  A.cross_section_init(0.4);

  FILE * out  = fopen ("farcs-cmz-out-DC5.txt" , "w+");

  for (int i=0;i<nlinhas;i++){
    fprintf(out,"%E\t%E\t",M[i],Z[i]);
    for (int j=1;j<=40;j++){
      Conc = conc(Z[i], M[i], cosmological_params);
      fprintf(out,"%E\t",f_arcs_II(Z[i], M[i], Conc, m_lim, cosmological_params, A));

    }
  fprintf(out,"\n");
  //fim();
  }

  //fim();
}

//use the data in farcs-cmz-out-fnal-40.txt to compute N_arcs, Delta_N_arcs etc in bins of M
void statM(){

  int nlinhas = 247697;

  //matriz que recebera os dados do arquivo farcs-cmz-out-fnal-40.txt, saida[linha][coluna]
  double **saida;
  saida = (double**) malloc(nlinhas*sizeof(double**));
  for(int i=0; i<nlinhas; i++){
    saida[i] = (double*) malloc(41*sizeof(double*));
  }


  FILE *FileIn  = fopen("farcs-cmz-out-fnal-40.txt","r"); //abre o arquivo de entrada
  DataFile2(FileIn,20, nlinhas, saida);

  //for(int i=0;i<42;i++) printf("%E ",saida[99999][i]);

  double *farcs = (double*) malloc(39*sizeof(double));
  for(int i=0;i<40;i++) farcs[i]=0.0;

  double m_min = 0+0.4950000;
  double m_max = 1.0+0.4950000;
  double m_medd = m_min/2+m_max/2;

  double mean, var, std;

  double n_halos_bin = 0.0;

  for(int k=0;k<21;k++){

    for(int i=0;i<nlinhas;i++){
      if(saida[i][0]>m_min && saida[i][0]<=m_max){
        for(int j=0;j<40;j++) farcs[j]+=saida[i][j+2];
        n_halos_bin++;
      }
    }
    //printf("m_min=%.1f, m_med=%f, m_max=%.1f\n",m_min,m_medd,m_max);
    //for(int j=0;j<40;j++) printf("%E, ",farcs[j]);
    //printf("\n");
    mean     = gsl_stats_mean(farcs, 1, 40);
    var = gsl_stats_variance(farcs, 1, 40);
    std  =  gsl_stats_sd(farcs, 1, 40);
    printf("mas_med= %.1f ,media= %.2f ,var= %.2f ,s1= %.2f ,s2= %.2f ,s3= %.2f ,sT= %.2f\n",
           m_medd, mean, var, std, mean*sqrt(n_halos_bin)/n_halos_bin, sqrt(mean), sqrt(mean) + std + mean*sqrt(n_halos_bin)/n_halos_bin);




    m_min = m_min + 1.0;
    m_max = m_max + 1.0;
    m_medd = (m_min + m_max)/2.0;
    for(int i=0;i<40;i++) farcs[i]=0.0;
    n_halos_bin = 0;
  }
 // for(int j=0;j<40;j++) printf("%E\n",farcs[j]);

}


//use the data in farcs-cmz-out-fnal-40.txt to compute N_arcs, Delta_N_arcs etc in bins of Z
void statZ(){

  int nlinhas = 247697;

  //matriz que recebera os dados do arquivo farcs-cmz-out-fnal-40.txt, saida[linha][coluna]
  double **saida;
  saida = (double**) malloc(nlinhas*sizeof(double**));
  for(int i=0; i<nlinhas; i++){
    saida[i] = (double*) malloc(41*sizeof(double*));
  }


  FILE *FileIn  = fopen("farcs-cmz-out-fnal-40.txt","r"); //abre o arquivo de entrada
  DataFile2(FileIn,20, nlinhas, saida);

  //for(int i=0;i<42;i++) printf("%E ",saida[99999][i]);

  double *farcs = (double*) malloc(39*sizeof(double));
  for(int i=0;i<40;i++) farcs[i]=0.0;

  double zbin = 0.1;
  double z_min_cat = 9.589870E-03;
  double z_max_cat = 1.453580E+00;

  double z_min = z_min_cat;
  double z_max = z_min_cat+zbin;
  double z_medd = z_min/2+z_max/2;

  double mean, var, std;

  double n_halos_bin = 0.0;

  int n_bins_z = (z_max_cat-z_min_cat)/zbin;
  for(int k=0;k<n_bins_z;k++){

    for(int i=0;i<nlinhas;i++){
      if(saida[i][1]>z_min && saida[i][1]<=z_max){
        for(int j=0;j<40;j++) farcs[j]+=saida[i][j+2];
        n_halos_bin++;
      }
    }
    //printf("m_min=%.1f, m_med=%f, m_max=%.1f\n",m_min,m_medd,m_max);
    //for(int j=0;j<40;j++) printf("%E, ",farcs[j]);
    //printf("\n");
    mean     = gsl_stats_mean(farcs, 1, 40);
    var = gsl_stats_variance(farcs, 1, 40);
    std  =  gsl_stats_sd(farcs, 1, 40);
    printf("z_med= %.2f ,media= %.2f ,var= %.2f ,s1= %.2f ,s2= %.2f ,s3= %.2f ,s2+3= %.2f\n",
           z_medd, mean, var, std, mean*sqrt(n_halos_bin)/n_halos_bin, sqrt(mean), sqrt(mean) + std + mean*sqrt(n_halos_bin)/n_halos_bin);




    z_min = z_min + zbin;
    z_max = z_max + zbin;
    z_medd = (z_min + z_max)/2.0;
    for(int i=0;i<40;i++) farcs[i]=0.0;
    n_halos_bin = 0;
  }
 // for(int j=0;j<40;j++) printf("%E\n",farcs[j]);

}

//use the function c(M,Z) to compute N_arcs, Delta_N_arcs etc in bins of M
//statistica para apenas uma realizacao do catalogo
void stat_bf_M(){

  int nlinhas = 3402;//247697;

  //matriz que recebera os dados do arquivo farcs-cmz-out-fnal-40.txt, saida[linha][coluna]
  double **saida;
  saida = (double**) malloc(nlinhas*sizeof(double**));
  for(int i=0; i<nlinhas; i++){
    saida[i] = (double*) malloc(1*sizeof(double*));
  }


  FILE *FileIn  = fopen("arcs_in_halos-comp-AddArcs-Z-out.txt","r");//; //abre o arquivo de entrada

  DataFile21(FileIn,20, nlinhas, saida);

  //for(int i=0;i<42;i++) printf("%E ",saida[99999][i]);

  double *farcs = (double*) malloc(1*sizeof(double));
  for(int i=0;i<1;i++) farcs[i]=0.0;
  double mbin = 2.0;
  double m_min_cat = 0.50022;//9.589870E-03;
  double m_max_cat = 9.1595;//1.453580E+00;

  double m_min = m_min_cat;
  double m_max = m_min_cat+mbin;
  double m_medd = m_min/2+m_max/2;

  //double mean, var, std;

  double n_halos_bin = 0.0;

  int n_bins_m = (m_max_cat-m_min_cat)/mbin+1;
  for(int k=0;k<n_bins_m;k++){

    for(int i=0;i<nlinhas;i++){
      if(saida[i][0]>m_min && saida[i][0]<=m_max){
        for(int j=0;j<1;j++) farcs[j]+=saida[i][j+2];
        n_halos_bin++;
      }
    }
    //printf("m_min=%.1f, m_med=%f, m_max=%.1f\n",m_min,m_medd,m_max);
    //for(int j=0;j<40;j++) printf("%E, ",farcs[j]);
    //printf("\n");
    //mean     = gsl_stats_mean(farcs, 1, 40);
    //var = gsl_stats_variance(farcs, 1, 40);
    //std  =  gsl_stats_sd(farcs, 1, 40);
    //printf("mas_med= %.1f ,media= %.2f ,var= %.2f ,s1= %.2f ,s2= %.2f ,s3= %.2f ,sT= %.2f\n",
    //       m_medd, mean, var, std, mean*sqrt(n_halos_bin)/n_halos_bin, sqrt(mean), sqrt(mean) + std + mean*sqrt(n_halos_bin)/n_halos_bin);
  printf("%.2f %E\n",m_medd,*farcs);



    m_min = m_min + 1.0;
    m_max = m_max + 1.0;
    m_medd = (m_min + m_max)/2.0;
    for(int i=0;i<1;i++) farcs[i]=0.0;
    n_halos_bin = 0;
  }
 // for(int j=0;j<40;j++) printf("%E\n",farcs[j]);

}


//use the function c(M,Z) to compute N_arcs, Delta_N_arcs etc in bins of M
//statistica para apenas uma realizacao do catalogo
void stat_bf_Z(){

  int nlinhas = 3402;//247697;

  //matriz que recebera os dados do arquivo farcs-cmz-out-fnal-40.txt, saida[linha][coluna]
  double **saida;
  saida = (double**) malloc(nlinhas*sizeof(double**));
  for(int i=0; i<nlinhas; i++){
    saida[i] = (double*) malloc(1*sizeof(double*));
  }


  FILE *FileIn  = fopen("arcs_in_halos-comp-AddArcs-Z-out.txt","r");//("farcs-cmz-out-best-fit.txt","r"); //abre o arquivo de entrada

  DataFile21(FileIn,20, nlinhas, saida);

  //for(int i=0;i<42;i++) printf("%E ",saida[99999][i]);

  double *farcs = (double*) malloc(1*sizeof(double));
  for(int i=0;i<1;i++) farcs[i]=0.0;
  double zbin = 0.3;
  double z_min_cat = 0.03900000;//9.589870E-03;
  double z_max_cat = 1.37020000;//1.453580E+00;

  double z_min = z_min_cat;
  double z_max = z_min_cat+zbin;
  double z_medd = z_min/2+z_max/2;

  //double mean, var, std;

  double n_halos_bin = 0.0;

  int n_bins_z = (z_max_cat-z_min_cat)/zbin+1;
  for(int k=0;k<n_bins_z;k++){

    for(int i=0;i<nlinhas;i++){
      if(saida[i][1]>z_min && saida[i][1]<=z_max){
        for(int j=0;j<1;j++){ farcs[j]+=saida[i][j+2]; }
        n_halos_bin++;
      }
    }
    //printf("z_min=%f, z_med=%f, z_max=%.1f\n",z_min,z_medd,z_max);
    //for(int j=0;j<40;j++) printf("%E, ",farcs[j]);
    //printf("\n");
    //mean     = gsl_stats_mean(farcs, 1, 40);
    //var = gsl_stats_variance(farcs, 1, 40);
    //std  =  gsl_stats_sd(farcs, 1, 40);
    //printf("mas_med= %.1f ,media= %.2f ,var= %.2f ,s1= %.2f ,s2= %.2f ,s3= %.2f ,sT= %.2f\n",
    //       m_medd, mean, var, std, mean*sqrt(n_halos_bin)/n_halos_bin, sqrt(mean), sqrt(mean) + std + mean*sqrt(n_halos_bin)/n_halos_bin);
  printf("%.2f %E\n",z_medd,*farcs);



    z_min = z_min + zbin;
    z_max = z_max + zbin;
    z_medd = (z_min + z_max)/2.0;
    for(int i=0;i<1;i++) farcs[i]=0.0;
    n_halos_bin = 0;
  }
 // for(int j=0;j<40;j++) printf("%E\n",farcs[j]);

}

void stat_farcs_cosmology(){

  int nlinhas = 3402;//247697;//
  double *Z = (double*) malloc(nlinhas*sizeof(double));
  double *M = (double*) malloc(nlinhas*sizeof(double));
  double *f_arcs = (double*) malloc(nlinhas*sizeof(double));

  FILE *FileIn  = fopen("farcs-cmz-out-DC5-Om03-Ol07.txt","r"); //abre o arquivo de entrada
  DataFile3Coluns(FileIn, nlinhas, M, Z, f_arcs);

  double Z_max = *std::max_element(Z, Z+nlinhas);
  //double M_max = *std::max_element(M, M+nlinhas);

  double Z_min = *std::min_element(Z, Z+nlinhas);
  //double M_min = *std::min_element(M, M+nlinhas);

  //printf("%E %E %E %E\n",Z_max, Z_min, M_max, M_min);
  //printf("%E %E %E\n", M[412], Z[412], f_arcs[412]);

  double M_med = 1;
  double M_bin = 0.025;
  double M_min_bin = M_med-M_bin;
  double M_max_bin = M_med+M_bin;

  double Z_bin = 0.2;
  int n_bins = ((Z_max-Z_min)/Z_bin);
  double Z_med = Z_min + Z_bin/2.0;
  double Z_min_bin = Z_med-Z_bin/2.0;
  double Z_max_bin = Z_med+Z_bin/2.0;

  FILE *out  = fopen("farcs-cosmology-err1.txt","w+");

  for(int j=1;j<=n_bins+1;j++){

    int n_halos_bin = 0;
    for (int i=0;i<nlinhas;i++){
      if( M[i]>=M_min_bin && M[i]<M_max_bin && Z[i]>=Z_min_bin && Z[i]<Z_max_bin) n_halos_bin++;//printf("%E %E %E\n", M[i], Z[i], farcs[i]);
    }
    //printf("%i\n",n_halos_bin);
    double *vec_f_arcs_bin = (double*) malloc(n_halos_bin*sizeof(double));
    int integer=0;
    for (int i=0;i<nlinhas;i++){
      if( M[i]>=M_min_bin && M[i]<M_max_bin && Z[i]>=Z_min_bin && Z[i]<Z_max_bin){
        vec_f_arcs_bin[integer] = f_arcs[i];
        integer++;
      }
    }

    double n_arcs_bin = 0.0;
    for(int i=0;i<n_halos_bin;i++) n_arcs_bin+=vec_f_arcs_bin[i];
    //n_arcs_bin*=5000/220;
    //printf("%E\n",n_arcs_bin);

    double f_arcs_mean_bin = n_arcs_bin/n_halos_bin;
    double delta_f_arcs_bin = (n_arcs_bin<1.0)? 1.0/n_halos_bin : sqrt(n_arcs_bin)/n_halos_bin;
    printf("Z_med=%f, Z_min=%f, Z_max=%f\n",Z_med,Z_min_bin,Z_max_bin);
    printf("n_halos=%i, n_arcs=%E, f_arcs_bin=%E+-%E\n\n",n_halos_bin,n_arcs_bin,f_arcs_mean_bin,delta_f_arcs_bin);

    fprintf(out,"%E %E %E\n",Z_med,f_arcs_mean_bin,delta_f_arcs_bin);

    Z_med += Z_bin;
    Z_min_bin = Z_med-Z_bin/2.0;
    Z_max_bin = Z_med+Z_bin/2.0;

  }




  CrossSection A;
  A.cross_section_init(0.4);
  double Om = 0.4;
  double Ok = 0.0;
  double Ol = 0.6;
  double cosmological_params[] = {Om, Ok, Ol};

  double Conc = 0.0;
  double m_lim = 24.0;

  double zl = Z_min;
  double hzl = (Z_max-Z_min)/100.0;


  FILE *out2  = fopen("farcs-cosmology-Om04.txt","w+");
  for (int j=0;j<100;j++){
    Conc = conc(zl, M_max_bin, cosmological_params);
    fprintf(out2,"%E %E\n",zl,f_arcs_II( zl, M_max_bin, Conc, m_lim, cosmological_params, A));
    zl+=hzl;
  }
  zl = Z_min;
  fprintf(out2,"\n\n");
  for (int j=0;j<100;j++){
    Conc = conc(zl, M_min_bin, cosmological_params);
    fprintf(out2,"%E %E\n",zl,f_arcs_II( zl, M_min_bin, Conc, m_lim, cosmological_params, A));
    zl+=hzl;
  }
  zl = Z_min;


  Om = 0.3;
  Ok = 0.0;
  Ol = 0.7;
  cosmological_params[0] = Om;
  cosmological_params[1] = Ok;
  cosmological_params[2] = Ol;
  FILE *out3  = fopen("farcs-cosmology-Om03.txt","w+");
  for (int j=0;j<100;j++){
    Conc = conc(zl, M_max_bin, cosmological_params);
    fprintf(out3,"%E %E\n",zl,f_arcs_II( zl, M_max_bin, Conc, m_lim, cosmological_params, A));
    zl+=hzl;
  }
  zl = Z_min;
  fprintf(out3,"\n\n");
  for (int j=0;j<100;j++){
    Conc = conc(zl, M_min_bin, cosmological_params);
    fprintf(out3,"%E %E\n",zl,f_arcs_II( zl, M_min_bin, Conc, m_lim, cosmological_params, A));
    zl+=hzl;
  }
  zl = Z_min;

  Om = 0.2;
  Ok = 0.0;
  Ol = 0.8;
  cosmological_params[0] = Om;
  cosmological_params[1] = Ok;
  cosmological_params[2] = Ol;
  FILE *out4  = fopen("farcs-cosmology-Om02.txt","w+");
  for (int j=0;j<100;j++){
    Conc = conc(zl, M_max_bin, cosmological_params);
    fprintf(out4,"%E %E\n",zl,f_arcs_II( zl, M_max_bin, Conc, m_lim, cosmological_params, A));
    zl+=hzl;
  }
  zl = Z_min;
  fprintf(out4,"\n\n");
  for (int j=0;j<100;j++){
    Conc = conc(zl, M_min_bin, cosmological_params);
    fprintf(out4,"%E %E\n",zl,f_arcs_II( zl, M_min_bin, Conc, m_lim, cosmological_params, A));
    zl+=hzl;
  }

}


double Angular_Distance(double z_init, double z_final, double cosmological_params[]){

  double Ii = I(z_init, z_final, cosmological_params);
  //double c = 299792.458; //  km/s;
  //double H_0 = 2.33336E-18; //1/seconds

  return (z_final + 1.0) * Ii;
}

void Graficos(){
  CrossSection A;
  A.cross_section_init(0.4);
  double Om = 0.3;
  double Ok = 0.0;
  double Ol = 0.7;
  double cosmological_params[] = {Om, Ok, Ol};
  
  double zl = 0.0;
  double m_lim = 24.0;
  double M = 0.0;
  double Conc = 0.0;

  M = 1.0;
  Conc = conc(zl, M, cosmological_params);

  zl=0.0;
  FILE *File  = fopen("grafico1.txt","w+");
  printf("M=%f, m_lim=%f\n",M,m_lim);
  for(int i=0;i<500;i++){
    zl+=1.5/500.0;
    fprintf(File,"%E\t%E\n",zl,f_arcs_II( zl, M, Conc, m_lim, cosmological_params, A));
  }


  M+=1.0;
  zl=0.0;
  FILE *File2  = fopen("grafico2.txt","w+");
  printf("M=%f, m_lim=%f\n",M,m_lim);
  for(int i=0;i<500;i++){
    zl+=1.5/500.0;
    fprintf(File2,"%E\t%E\n",zl,f_arcs_II( zl, M, Conc, m_lim, cosmological_params, A));
  }

  M+=1.0;
  zl=0.0;
  FILE *File3  = fopen("grafico3.txt","w+");
  printf("M=%f, m_lim=%f\n",M,m_lim);
  for(int i=0;i<500;i++){
    zl+=1.5/500.0;
    fprintf(File3,"%E\t%E\n",zl,f_arcs_II( zl, M, Conc, m_lim, cosmological_params, A));
  }

  M+=1.0;
  zl=0.0;
  FILE *File4  = fopen("grafico4.txt","w+");
  printf("M=%f, m_lim=%f\n",M,m_lim);
  for(int i=0;i<500;i++){
    zl+=1.5/500.0;
    fprintf(File4,"%E\t%E\n",zl,f_arcs_II( zl, M, Conc, m_lim, cosmological_params, A));
  }

  M+=1.0;
  zl=0.0;
  FILE *File5  = fopen("grafico5.txt","w+");
  printf("M=%f, m_lim=%f\n",M,m_lim);
  for(int i=0;i<500;i++){
    zl+=1.5/500.0;
    fprintf(File5,"%E\t%E\n",zl,f_arcs_II( zl, M, Conc, m_lim, cosmological_params, A));
  }

  M+=1.0;
  zl=0.0;
  FILE *File6  = fopen("grafico6.txt","w+");
  printf("M=%f, m_lim=%f\n",M,m_lim);
  for(int i=0;i<500;i++){
    zl+=1.5/500.0;
    fprintf(File6,"%E\t%E\n",zl,f_arcs_II( zl, M, Conc, m_lim, cosmological_params, A));
  }

  M+=1.0;
  zl=0.0;
  FILE *File7  = fopen("grafico7.txt","w+");
  printf("M=%f, m_lim=%f\n",M,m_lim);
  for(int i=0;i<500;i++){
    zl+=1.5/500.0;
    fprintf(File7,"%E\t%E\n",zl,f_arcs_II( zl, M, Conc, m_lim, cosmological_params, A));
  }

  M+=1.0;
  zl=0.0;
  FILE *File8  = fopen("grafico8.txt","w+");
  printf("M=%f, m_lim=%f\n",M,m_lim);
  for(int i=0;i<500;i++){
    zl+=1.5/500.0;
    fprintf(File8,"%E\t%E\n",zl,f_arcs_II( zl, M, Conc, m_lim, cosmological_params, A));
  }

  M+=1.0;
  zl=0.0;
  FILE *File9  = fopen("grafico9.txt","w+");
  printf("M=%f, m_lim=%f\n",M,m_lim);
  for(int i=0;i<500;i++){
    zl+=1.5/500.0;
    fprintf(File9,"%E\t%E\n",zl,f_arcs_II( zl, M, Conc, m_lim, cosmological_params, A));
  }

  M+=1.0;
  zl=0.0;
  FILE *File10  = fopen("grafico10.txt","w+");
  printf("M=%f, m_lim=%f\n",M,m_lim);
  for(int i=0;i<500;i++){
    zl+=1.5/500.0;
    fprintf(File10,"%E\t%E\n",zl,f_arcs_II( zl, M, Conc, m_lim, cosmological_params, A));
  }
}

void stat_comp_addarcs(){
  int nlinhas = 3402;//247697;//
  double *zl = (double*) malloc(nlinhas*sizeof(double));
  double *M = (double*) malloc(nlinhas*sizeof(double));
  double *conc = (double*) malloc(nlinhas*sizeof(double));

  FILE *FileIn  = fopen("arcs_in_halos-comp-AddArcs-Z.txt","r"); //abre o arquivo de entrada
  DataFile3Coluns(FileIn, nlinhas, M, zl, conc);

  double farcs = 0.0;

  CrossSection A;
  A.cross_section_init(0.5);
  double Om = 0.3;
  double Ok = 0.0;
  double Ol = 0.7;
  double cosmological_params[] = {Om, Ok, Ol};
  
  double m_lim = 24.0;


  FILE *FileOut  = fopen("arcs_in_halos-comp-AddArcs-Z-out.txt","w+");
  for(int i=0;i<nlinhas;i++){
    farcs = f_arcs_II( zl[i], M[i], conc[i], m_lim, cosmological_params, A);
    printf("M=%E, zL=%E, conc=%E, f_arcs=%E\n",M[i], zl[i], conc[i], farcs);
    fprintf(FileOut,"%E %E %E\n",M[i], zl[i], farcs);

  }

}


int main(){
  double z_l = 0.0;
  double M = 5.0;

  double Om = 0.25;
  double Ok = 0.0;
  double Ol = 0.75;
  double cosmological_params[] = {Om, Ok, Ol};

  FILE * out1  = fopen ("C-M5-polynomial.txt" , "w+");
  for (int i=0;i<200;i++){
    fprintf(out1,"%E %E\n",z_l,conc(z_l, M, cosmological_params));
    z_l += 3.0/200.0;
  }
  fclose(out1);
//Graficos();

// stat_comp_addarcs();
//stat_bf_Z();
//stat_bf_M();
/*  double Om = 0.3;
  double Ok = 0.0;
  double Ol = 0.7;
  double cosmological_params[] = {Om, Ok, Ol};
  //printf("%.30E\n",Angular_Distance(0.0,0.3,cosmological_params));

  double M = 1.4248599303;
  double zl = 0.142704114318;
  double conc = 2.76073169708;

  for(int i=0;i<10;i++){
    double z = 3.0 / 999.0 * i;
    printf("%10.8f %20.15f\n",z,I(0.0, z, cosmological_params));
  }
*/
 // printf("%E\n",conc(0.5, 1,cosmological_params));
//  test();
  //calcula o argumento da integral f_arcs_II
 /* double zl = 0.5;
  double M = 1.0;
  double conc = 5;
  double Om=0.3;
  double Ol=0.0;
  double Ok=0.7;
  double m_lim = 24;

  double params[] = {zl,M,conc,m_lim,Om,Ol,Ok};

  CrossSection A;
  A.cross_section_init(0.4);

  int n_points = 100;
  double zs=zl;
  double zs_max = 5.0;
  double h = (zs_max-zs)/double(n_points);

  FILE *File  = fopen("farcsII-arg.txt","w+");
  for(int j = 0;j<=8;j++){

  for (int i=0;i<n_points;i++){
    fprintf(File,"%E\t%E\n",zs,f_arcs_II_arg(zs, params, A));
    zs+=h;
  }
  fprintf(File,"\n");
  zs=zl;
  M+=1;
  params[1] = M;
  }


  fclose(File);*/
 // Graficos();


  //stat_farcs_cosmology();


  //Conpute_Farcs_Cathalog_DC5();

  fim();
  return 0;
}
