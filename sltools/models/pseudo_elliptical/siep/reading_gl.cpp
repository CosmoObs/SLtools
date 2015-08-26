/*
 * Programa: reading.cpp
 * Descricao: Le de um arquivo, cujo nome sera fornecido pelo
 *            usuario, um conjunto de numeros reais e armazena
 *            em um vetor.
 *            O tamanho maximo do vetor e dado pela constante
 *            TAM_MAX.
 *            A quantidade de numeros no arquivo varia entre 0
 *            e TAM_MAX.
 * Autor:     Habib 
 * Data: 20/09/2004
*/
#include <stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

#define TAM_MAX 10000
int main () {

  // system("./cuting.sh");
	system("grep -v \\# crit.txt > lens_curves.dat");
	

  FILE *entrada1;
  entrada1=fopen("lens_curves.dat","r");
  FILE *saida1;
  saida1=fopen("tan_crit.dat","w");
  FILE *saida2;
  saida2=fopen("rad_crit.dat","w");
  FILE *saida3;
  saida3=fopen("tan_caust.dat","w");
  FILE *saida4;
  saida4=fopen("rad_caust.dat","w");
  
  int i = 0;
  float x1[TAM_MAX], y1[TAM_MAX],x2[TAM_MAX],y2[TAM_MAX];
  float u1[TAM_MAX], v1[TAM_MAX],u2[TAM_MAX],v2[TAM_MAX];

  if (entrada1== NULL) 
    {
      printf("Nao foi possivel abrir o arquivo.\n");
      exit(1);
     }
     while (! feof(entrada1) )
      {
        i++;
        fscanf(entrada1,"%f %f %f %f %f %f %f %f",&x1[i],&y1[i],&u1[i],&v1[i],&x2[i],&y2[i],&u2[i],&v2[i]);
      }
        int n_dats=i-1;
        printf("We have  %d data points\n", n_dats);

   double x1tan=x1[1];
   double x2tan=y1[1];
   int n_sep=0;

  for(int i1=2;i1<=n_dats;i1++){
//     if((fabs(x1[i1])==x1tan) && (fabs(y1[i1])==x2tan)){
    if((x2[i1]==x1tan) &&  (y2[i1]==x2tan)){
      n_sep=i1;
      printf("the key number is %i\n",n_sep);
      break;
      }
    }
  
    n_sep=n_sep+1;

  for(int i2=1;i2<=n_dats;i2++){
    if(i2<n_sep){
      fprintf(saida2,"%f %f\n",x1[i2],y1[i2]);
      fprintf(saida2,"%f %f\n", x2[i2],y2[i2]);
      fprintf(saida4,"%f %f\n",u1[i2],v1[i2]);
      fprintf(saida4,"%f %f\n", u2[i2],v2[i2]);
    }else{
      fprintf(saida1,"%f %f\n",x1[i2],y1[i2]);
      fprintf(saida1,"%f %f\n", x2[i2],y2[i2]);
      fprintf(saida3,"%f %f\n",u1[i2],v1[i2]);
      fprintf(saida3,"%f %f\n", u2[i2],v2[i2]);}
    }  
    fclose(saida1);
    fclose(saida2);
    fclose(saida3);
    fclose(saida4);
    system("xmgrace -view 0.15 0.15 0.85 0.85 tan_caust.dat rad_caust.dat -saveall source_plane.agr");
    system("xmgrace -view 0.15 0.15 0.85 0.85 tan_crit.dat  rad_crit.dat -saveall lens_plane.agr");
  return 0.0;
}
