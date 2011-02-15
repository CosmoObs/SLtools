#ifndef CROSSSECTION_H
#define CROSSSECTION_H

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#include <gsl/gsl_integration.h>

int DataFile(FILE *arquivo, int NumeroLinhas, int LinhaInicial, double *saida1, double *saida2, double *saida3, double *saida4, double *saida5)
{//le um arquivo de 5 colunas e passa para vetores, iniciando na linha 'LinhaInicial'

  char IN[500];
  double out[4];

  //printf("Lendo dados do arquivo.");


  for(int i=1;i<=LinhaInicial;i++){//anda no arquivo ate a linha LinhaInicial
    fgets(IN,500,arquivo);
    sscanf(IN,"%lg %lg %lg %lg %lg",
               &out[0] ,&out[1] ,&out[2] ,&out[3] ,&out[4]);
  }

  for(int i=1;i<(NumeroLinhas-LinhaInicial+1);i++)
  {
    fgets(IN,500,arquivo);
    sscanf(IN,"%lg %lg %lg %lg %lg",
               &out[0] ,&out[1] ,&out[2] ,&out[3] ,&out[4]);

    saida1[i-1] = out[0];
    saida2[i-1] = out[1];
    saida3[i-1] = out[2];
    saida4[i-1] = out[3];
    saida5[i-1] = out[4];
    //if (i<10) printf("%E, %E, %E, %E, %E\n",saida1[i-1],saida2[i-1],saida3[i-1],saida4[i-1],saida5[i-1]);
  }

  //printf("->Leitura finalizada.\n");
  rewind(arquivo);
  return 0;
}

int DataFile(FILE *arquivo, int NumeroLinhas, double *saida1, double *saida2, double *saida3, double *saida4, double *saida5)
{//le um arquivo de 5 colunas e passa para vetores

  char IN[500];
  double out[4];

  //printf("Lendo dados do arquivo.");

  for(int i=1;i<=NumeroLinhas;i++)
  {
    fgets(IN,500,arquivo);
    sscanf(IN,"%lg %lg %lg %lg %lg",
               &out[0] ,&out[1] ,&out[2] ,&out[3] ,&out[4]);

    saida1[i-1] = out[0];
    saida2[i-1] = out[1];
    saida3[i-1] = out[2];
    saida4[i-1] = out[3];
    saida5[i-1] = out[4];
    //if (i<10) printf("%E, %E, %E, %E, %E\n",saida1[i-1],saida2[i-1],saida3[i-1],saida4[i-1],saida5[i-1]);
  }

  //printf("->Leitura finalizada.\n");
  rewind(arquivo);
  return 0;
}

int DataFile(FILE *arquivo,int coluna, int NumeroLinhas, double *saida)
{//le um arquivo passa a coluna 'coluna' para um vetores

  char IN[500];
  double out[4];

  //printf("Lendo dados do arquivo.");

  for(int i=1;i<=NumeroLinhas;i++)
  {
    fgets(IN,500,arquivo);
    sscanf(IN,"%lg %lg %lg %lg %lg",
               &out[0] ,&out[1] ,&out[2] ,&out[3] ,&out[4]);

    saida[i-1] = out[coluna-1];
  }

  //printf("->Leitura finalizada.\n");
  rewind(arquivo);
  return 0;
}

class CrossSection{//retorna o valor da secao de choque e da magnificacao media lidos de um arquivo
  public:
    //Construtores
    CrossSection(){
      nlinhas = 5001;
    }

    double *ks_vec      ;// vetor a ser preenchido com dados do arquivo de entrerda
    double *q_vec       ;// vetor a ser preenchido com dados do arquivo de entrerda
    double *area_vec    ;// vetor a ser preenchido com dados do arquivo de entrerda
    double *mag_min_vec ;// vetor a ser preenchido com dados do arquivo de entrerda
    double *mag_med_vec ;// vetor a ser preenchido com dados do arquivo de entrerda

  //Funcoes
    void cross_section_init(double elip_in){//Le o arquivo de entrada e inicializa as funcoes de interpolacao

      int linha_inicial=0; //primeira linha onde ks != 0

      if(elip_in < 0.05){
        FileIn = fopen("q10_2.txt","r");
        //printf("elip_in = 0.0\n");
      }
      else if(elip_in >= 0.05 && elip_in < 0.15){
        FileIn = fopen("q09_2.txt","r");
        //printf("elip_in = 0.1\n");
      }
      else if(elip_in >= 0.15 && elip_in < 0.25){
        FileIn = fopen("q08_2.txt","r");
       // printf("elip_in = 0.2\n");
      }
      else if(elip_in >= 0.25 && elip_in < 0.35){
        FileIn = fopen("q07_2.txt","r");
        //printf("elip_in = 0.3\n");
      }
      else if(elip_in >= 0.35 && elip_in < 0.45){
        FileIn = fopen("q06_2.txt","r");
        //printf("elip_in = 0.4\n");
      }
      else if(elip_in >= 0.45 && elip_in < 0.55){
        FileIn = fopen("q05_2.txt","r");
        //printf("elip_in = 0.5\n");
      }
      else if(elip_in >= 0.55 && elip_in < 0.65){
        FileIn = fopen("q04_2.txt","r");
        //printf("elip_in = 0.6\n");
      }
      else if(elip_in >= 0.65 && elip_in < 0.75){
        FileIn = fopen("q03_2.txt","r");
        //printf("elip_in = 0.7\n");
      }
      else if(elip_in >= 0.75 && elip_in < 0.85){
        FileIn = fopen("q02_2.txt","r");
        //printf("elip_in = 0.8\n");
      }
      else{
        FileIn = fopen("q10_2.txt","r");
        printf("You choose elip_in = %f\nelip_in must have the following values: 0.0--0.85.\nNow setting elip_in = 0.0\n",elip_in);
      }

      ks_vec      = (double*) malloc(nlinhas*sizeof(double));// aloca o espaco de memoria para o vetor ks_vec
      q_vec       = (double*) malloc(nlinhas*sizeof(double));// aloca o espaco de memoria para o vetor q_vec
      area_vec    = (double*) malloc(nlinhas*sizeof(double));// aloca o espaco de memoria para o vetor area_vec
      mag_min_vec = (double*) malloc(nlinhas*sizeof(double));// aloca o espaco de memoria para o vetor mag_min_vec_vec
      mag_med_vec = (double*) malloc(nlinhas*sizeof(double));// aloca o espaco de memoria para o vetor mag_med_vec_vec

      DataFile(FileIn, nlinhas, ks_vec, q_vec, area_vec, mag_min_vec, mag_med_vec);
      //enche os vetors com os dados do arquivo de entrada, partindo da primeira linha

      for(int i=0;i<nlinhas;i++){//encontra em que linha aparece o primeiro ks!=0
        if(area_vec[i] > 0.0){
          linha_inicial = i-1;
          ks_min = ks_vec[i-1];
          //printf("%E\n",ks_min);
          break;
          }
      }

      free(area_vec); 

      ks_vec     =(double*) malloc( (nlinhas-linha_inicial)*sizeof(double));// aloca o espaco de memoria para o vetor ks_vec
      q_vec      =(double*) malloc( (nlinhas-linha_inicial)*sizeof(double));// aloca o espaco de memoria para o vetor q_vec
      area_vec   =(double*) malloc( (nlinhas-linha_inicial)*sizeof(double));// aloca o espaco de memoria para o vetor area_vec
      mag_min_vec=(double*) malloc( (nlinhas-linha_inicial)*sizeof(double));// aloca o espaco de memoria para o vetor mag_min_vec
      mag_med_vec=(double*) malloc( (nlinhas-linha_inicial)*sizeof(double));// aloca o espaco de memoria para o vetor mag_med_vec

      //printf("Lendo dados");      
      DataFile(FileIn, nlinhas, linha_inicial, ks_vec, q_vec, area_vec, mag_min_vec, mag_med_vec);
      //enche os vetors com os dados do arquivo de entrada, partindo da linha linha_inicial
      //printf(" -> Dados lidos\n");
      
      acc_area = gsl_interp_accel_alloc();
      spline_area = gsl_spline_alloc(gsl_interp_linear, nlinhas-linha_inicial);
      gsl_spline_init(spline_area, ks_vec, area_vec, nlinhas-linha_inicial);

      acc_mag_med = gsl_interp_accel_alloc();
      spline_mag_med = gsl_spline_alloc(gsl_interp_linear, nlinhas-linha_inicial);
      gsl_spline_init(spline_mag_med, ks_vec, mag_med_vec, nlinhas-linha_inicial);

      acc_mag_min = gsl_interp_accel_alloc();
      spline_mag_min = gsl_spline_alloc(gsl_interp_linear, nlinhas-linha_inicial);
      gsl_spline_init(spline_mag_min, ks_vec, mag_min_vec, nlinhas-linha_inicial);
    }


    double area(double ks_in){//retorna a area(secao de choque) interpolada 
      if(ks_in <= ks_min) return 0.0;
      return gsl_spline_eval(spline_area, ks_in, acc_area);
    }

    double mag_med(double ks_in){//retorna a magnificacao media interpolada
      if(ks_in <= ks_min) return 1.0;
      return gsl_spline_eval(spline_mag_med, ks_in, acc_mag_med);
    }

    double mag_min(double ks_in){//retorna a magnificacao minima interpolada
      if(ks_in <= ks_min) return 0.0;
      return gsl_spline_eval(spline_mag_min, ks_in, acc_mag_min);
    }




 // private:
    double elip, ks, zL; //elipticidade, ks e desvio para o vermelho da lente
    //double M, C; //Massa a parametro de concentracao
    int  nlinhas; //numero de linhas do arquivo de entrada

    double ks_min; //valor minimo de ks em que a secao de choque e considerada diferente de zero


    FILE *FileIn; //aqruivo de entrada
    
    gsl_interp_accel *acc_area; //acelerador para interpolar a funcao area(ks)
    gsl_spline *spline_area;

    gsl_interp_accel *acc_mag_med; //acelerador para interpolar a funcao mag_med(ks)
    gsl_spline *spline_mag_med;

    gsl_interp_accel *acc_mag_min; //acelerador para interpolar a funcao mag_min(ks)
    gsl_spline *spline_mag_min;




};

double DimCrossSection(double zl, double zs, double M, double conc, double cosmological_params[], CrossSection A){
  double XS = x_s(zl,M,conc,cosmological_params);
  double KS = kappa_s(zl,zs,M,conc,cosmological_params);
  //printf("%.20E\n",A.area(KS)*XS*XS);
  return A.area(KS)*XS*XS;
}

#endif
