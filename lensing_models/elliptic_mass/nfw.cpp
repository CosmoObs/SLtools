#include <cmath>          //funcoes matematicas
#include <ctime>          //funcoes de relogio
#include <fstream>        //funcoes de entrada e saida para arquivos
#include <cstdlib>         //system functions
#include <iostream>


#include "lente.h"
#include "integrais.h"
#include "propriedades.h"
#include "curvas.h"
#include "area.h"
//#include "graficos_xmgrace.h"


#define PI 3.1415926535897932384626433832795      //valor de PI

double area(double q, double ks, double RAZ,int iR, int iT, double* mag_min, double* x_mag_min, double* y_mag_min, double* mag_med);
double area(double q, double ks, double RAZ,int iR, int iT, double* mag_min, double* x_mag_min, double* y_mag_min, double* mag_med, double Mag_Min);

void Grid(int N_ks, int N_q, double ks_min, double ks_max, double q_min, double q_max, double RAZ, char file[]);

double elip(double e, double ks)
{return e*(2.272 + 0.4*ks - 0.29*pow(ks,2) + 0.117*pow(ks,3) - 0.024*pow(ks,4) + 0.002*pow(ks,5) - (0.06206896551724138*(0.6049790635362423 + ks)*(301.38784021310244 + ks*(28.928354269797094 + ks))*e)/((0.9221773473613604 + ks)*(6.139891618155882 + ks)));};

//int teste(char palavra[]);

int main()
{
  double q=0.5;
//  double M200 = 5;
  double RAZ = 7.5;
  double ks = 9.556262E-01;//1.283045E-02;
  
  int I = 1;
  int Npt = 200;

  char file1[] = "saida1.txt";
  char file2[] = "saida2.txt";
  char file3[] = "saida3.txt";

  graficosF(q,ks,RAZ, I, Npt, file1, file2, file3);

/*  double *mag_min = new double;
  double *mag_med = new double;
  printf("%E\n",AREA_2(0.973,0.671,10,mag_min,mag_med));
*/
 // graficosL(1.0, 0.01, 10, 1,100, "out1.txt", "out2.txt", "out3.txt");
//CCL(0.6 ,  1.916834E-02, 1, 200, "out1.txt", "out1.txt");

 // Grid(100, 0, 0, 2.5, 1.0, 1.0, 10, "e00.txt");

  /* CCL(0.9 , 1.0 , 1 , 350 , "CCL1.txt", "CCL2.txt");
   CAF(0.9 , 1.0 , 1 , 350 , "CAF1.txt", "CAF2.txt");*/

  //Ks minimo
  /*double q,ks,RAZ;
  double *mag_min = new double;
  double *mag_med = new double;
  double *x_mag_min = new double;
  double *y_mag_min = new double;
  double Mag_Min, Mag_Max;
  double Mag_Min0, h;
  double AREA;
  int iR, iT;
  int Npts;

  q=0.6;
  ks=1.0;
  RAZ=10.0;

  Grid(5000,0,0.0,1.0,0.2,0.2,10, "grid02.txt");*/

  /*Mag_Min = 1.0;
  Mag_Min0 = Mag_Min;
  Mag_Max = 400.0;
  Npts = 200;

  h = pow(Mag_Max/Mag_Min,1.0/double(Npts));

  iR=550;
  iT=550;*/

  /*  FILE *saida1 = fopen("saida1.txt","w+");
    for(int i=0;i<=400;i++)
    {
      AREA = AREA_2(q,ks,RAZ,mag_min,mag_med);
      printf("ks=%E  q=%E  Area=%E  mag_min=%E  mag_med=%E\n",ks,q,AREA,*mag_min,*mag_med);

      fprintf(saida1,"%.10E\t%.10E\n",RAZ,AREA);
      RAZ+=0.5;
    }*/

  /*FILE *saida1 = fopen("ks10q04.txt","w+");
  area(q, ks, RAZ, iR, iT, mag_min, x_mag_min, y_mag_min, mag_med, Mag_Min);
  fprintf(saida1,"ks=%E\tq=%E\tmag_min=%E\tmag_med=%E\n",ks,q,*mag_min,*mag_med);

  for(int i=1;i<=(Npts+1);i++){
  AREA = area(q, ks, RAZ, iR, iT, mag_min, x_mag_min, y_mag_min, mag_med, Mag_Min);
  printf("%E\t%E\t%E\t%E\t%i\n", AREA, *mag_min ,*mag_med,Mag_Min,i);
  fprintf(saida1,"%E\t%E\n",Mag_Min,AREA);
  Mag_Min = Mag_Min0*pow(h,i);
  }*/

  /*for(int i=0;i<=40;i++){
  AREA = area(q, ks, RAZ, iR, iT, mag_min, x_mag_min, y_mag_min, mag_med, Mag_Min);
  printf("%E\t%E\t%E\t%E\n", AREA, *mag_min ,*mag_med,Mag_Min);
  fprintf(saida1,"%E\t%E\n",Mag_Min,AREA);
  Mag_Min += 0.5;
  }
  printf("ParteII");
  for(int i=0;i<=30;i++){
  AREA = area(q, ks, RAZ, iR, iT, mag_min, x_mag_min, y_mag_min, mag_med, Mag_Min);
  printf("%E\t%E\t%E\t%E\n", AREA, *mag_min ,*mag_med,Mag_Min);
  fprintf(saida1,"%E\t%E\n",Mag_Min,AREA);
  Mag_Min += 15.0;
  }
  fclose(saida1);*/


  /*ks=1.2;
  RAZ=10.0;

  q=0.3;
  graficosL(q, ks, RAZ, 4, 100, "saida1-1.txt", "saida2-1.txt", "saida3-1.txt");
  OUT = AREA_2(q, ks, RAZ, mag_min, mag_med);
  printf("%E\n",OUT);

  q=1.0;
  graficosL(q, ks, RAZ, 1, 100, "saida1-2.txt", "saida2-2.txt", "saida3-2.txt");
  OUT = AREA_2(q, ks, RAZ, mag_min, mag_med);
  if(isnan(OUT)) printf("NAN->");

  printf("%E\n",OUT);*/

  //Calcula o valor de kappa em L1=0,L/W=-10,L/W=+10 e faz a transformacao da elipticidade
  /*double q, qsigma, ks, T;
  double rmed,rmax,rmin;
  double kmin,kmed,kmax;
  double kmin2,kmed2,kmax2;
  //double esigma,e;

  T = 0*PI/4.0;
  qsigma=1.0;
  ks=0.5;

  FILE *out = fopen("saida.txt","w+");
  for(int j=0;j<=3;j++){
    fprintf(out,"\n\nks=%.2f, T=%.1f\n",ks,T/PI*180);
    fprintf(out,"esigma\t\t\te\t\t\t\tk(L1=0)\t\t\tk(L/W=-10)\t\tk(L/W=10)\t\tX(L1=0)\t\t\tX(L/W=-10)\t\tX(L/W=10)\n");
    for(int i=0;i<=30;i++){
      q = 1.0-elip(1.0-qsigma,ks);

      rmed = L1ZERO(T, q, ks);
      rmax = L1OL2RAZ(rmed,L2ZERO(T,q,ks),T,q,ks,10.0,1);
      rmin = L1OL2RAZ(rmed,L2ZERO(T,q,ks),T,q,ks,10.0,-1);

      //printf("rmin =%E\trmed =%E\trmax =%E\n",rmin,rmed,rmax);

      kmed = kCSI2(rmed*cos(T),rmed*sin(T),q,ks);
      kmin = kCSI2(rmin*cos(T),rmin*sin(T),q,ks);
      kmax = kCSI2(rmax*cos(T),rmax*sin(T),q,ks);

      //printf("kmin =%E\tkmed =%E\tkmax =%E\n",kmin,kmed,kmax);

      kmed2 = (fiXX(rmed*cos(T),rmed*sin(T),q,ks)+fiYY(rmed*cos(T),rmed*sin(T),q,ks))/2.0;
      kmin2 = (fiXX(rmin*cos(T),rmin*sin(T),q,ks)+fiYY(rmin*cos(T),rmin*sin(T),q,ks))/2.0;
      kmax2 = (fiXX(rmax*cos(T),rmax*sin(T),q,ks)+fiYY(rmax*cos(T),rmax*sin(T),q,ks))/2.0;

      //printf("kmin2=%E\tkmed2=%E\tkmax2=%E\n",kmin2,kmed2,kmax2);
      //printf("  d\%1=%.10f\t  d\%2=%.10f\t  d\%3=%.10f\tq=%.2f\n\n",fabs((kmin-kmin2)/kmin)*100,fabs((kmed-kmed2)/kmed)*100,fabs((kmax-kmax2)/kmax)*100,q);

      fprintf(out,"%E\t%E\t%E\t%E\t%E\t",1.0-qsigma,1.0-q,kmed/ks,kmin/ks,kmax/ks);
      fprintf(out,"%E\t%E\t%E\n",rmed*cos(T),rmin*cos(T),rmax*cos(T));

      qsigma-=0.01;
    }
  qsigma=1.0;
  printf("%E\n",ks);
  ks+=0.5;
  }
  fclose(out);*/



//Gera uma saida igual a do gravlens para o comando plotdef1.
/*FILE *saida1 = fopen("saida.txt","w+");
double x,y;

fprintf(saida1,"\n\n\n\n\n");
x=0.0;
y=0.0;
for(int i=0;i<=250;i++){
  for(int j=0;j<=250;j++){
    fprintf(saida1,"%-+.4E %-+.4E %-+.4E %-+.4E %-+.4E %-+.4E %-+.4E %-+.4E %-+.4E\n",x,y,0.0,fiX(x,y,q,ks),fiY(x,y,q,ks),fiXX(x,y,q,ks),fiYY(x,y,q,ks),fiXY(x,y,q,ks),fiXY(x,y,q,ks));
    //printf("%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\n",x,y,0.0,fiX(x,y,q,ks),fiY(x,y,q,ks),fiXX(x,y,q,ks),fiYY(x,y,q,ks),fiXY(x,y,q,ks),fiXY(x,y,q,ks));
    y+=0.02;
  }
  y=0.0;
  x+=0.02;
}*/


    std::cout << "tempo decorrido: " << double(clock())/double(CLOCKS_PER_SEC) << " segundos" << std::endl;
    std::cout << "FIM!" << std::endl;
    //system("mpg321 /home/gbcaminha/Gabriel/MP3/CDs\\ Completos/Cake/Cake\\ -\\ B-Sides\\ And\\ Rarities/\\Cake\\ -\\ B-Sides\\ And\\ Rarities\\ -\\ 01\\ -\\ War\\ Pigs.mp3");
    return 0;
}

void Grid(int N_ks, int N_q, double ks_min, double ks_max, double q_min, double q_max, double RAZ, char file[])
{
  if(ks_min > ks_max){//troca ks_min por ks_max se ks_min > ks_max
  double swap=ks_max;
  ks_max=ks_min;
  ks_min=swap;
  }
  if(q_min > q_max){//troca q_min por q_max se q_min > q_max
  double swap=q_max;
  q_max=q_min;
  q_min=swap;
  }

  //if(ks_min < 0.05)std::cout<<"Grid->ks_min < 0.05, fazendo ks_min = "<<double(ks_min=0.05)<<std::endl;
  //if(q_min < 0.05)std::cout<<"Grid->q_min < 0.05, fazendo q_min = "<<double(q_min=0.05)<<std::endl;

  double dq=(q_max-q_min)/double(N_q);
  double dks=(ks_max-ks_min)/double(N_ks);

  double q=q_max;
  double ks=ks_max;
  double Area=0.0;

  double *mag_min = new double;
  double *mag_med = new double;

  FILE *saida1 = fopen(file,"w+");

  fprintf(saida1,"ks\t\t\t\tq\t\t\t\tArea\t\t\tmag_min\t\t\tmag_med\n");
  for(int j=0;j<=N_ks;j++)
  {
    for(int i=0;i<=N_q;i++)
    {
    Area = AREA_2(q,ks,RAZ,mag_min,mag_med);
    printf("ks=%E  q=%E  Area=%E  mag_min=%E  mag_med=%E\n",ks,q,Area,*mag_min,*mag_med);
    if(isnan(Area)) Area=0.0;
    if(isnan(*mag_min)) *mag_min=0.0;
    if(isnan(*mag_med)) *mag_med=0.0;

    fprintf(saida1,"%.10E\t%.10E\t%.10E\t%.10E\t%.10E\n",ks,q,Area,*mag_min,*mag_med);
    q-=dq;
    }
  q = q_max;
  ks-=dks;
  }
  fclose(saida1);
}


/******************************************************************************/
double area(double q, double ks, double RAZ,int iR, int iT, double* mag_min, double* x_mag_min, double* y_mag_min, double* mag_med)
{
    double dr,r,dT,T;
    double AREA_fonte, AREA_lente, dA_lente;
    double rmin,rmax,rmed, rINF;
    double X,Y,DX,DY;
	double Det;

    r = 0.0;
    T = 0.0;
    dT = PI/(2.0*double(iT));
    AREA_lente = 0.0;
	AREA_fonte = 0.0;

    rmed = L1ZERO(T,q,ks);
    rmin = rmed;
    rmax = rmed;

    rmax = L1OL2RAZ(rmax,10.0,T,q,ks,RAZ,+1);
	*mag_min = 1.0/DET(rmax,0.0,q,ks);

    for (int i=0;i<iT;i++)
    {
		rINF = L2ZERO(T,q,ks);
		rmed = L1ZERO(T,q,ks);
		rmin = L1OL2RAZ(rmed,rINF,T,q,ks,RAZ,-1);
		rmax = L1OL2RAZ(rmed,10.0,T,q,ks,RAZ,+1);
        r = rmin;
        dr = fabs(rmax - rmin)/double(iR);
        for(int i=0;i<iR;i++)
        {
            X = r*cos(T);
            Y = r*sin(T);
            DX = dr*cos(T);
            DY = dr*sin(T);

			Det = fabs(DET(X+DX/2.0,Y+DY/2.0,q,ks));
			dA_lente = r*dr*dT;

			AREA_lente += dA_lente;
            AREA_fonte += dA_lente*Det;

            r += dr;
        }
        T += dT;
    }
	*mag_med = AREA_lente/AREA_fonte;
    return 4.0*AREA_fonte;
}


/******************************************************************************/
double area(double q, double ks, double RAZ,int iR, int iT, double* mag_min, double* x_mag_min, double* y_mag_min, double* mag_med, double Mag_Min)
{
    double dr,r,dT,T;
    double AREA_fonte, AREA_lente, dA_lente;
    double rmin,rmax,rmed, rINF;
    double X,Y,DX,DY;
	double Det;

    r = 0.0;
    T = 0.0;
    dT = PI/(2.0*double(iT));
    AREA_lente = 0.0;
	AREA_fonte = 0.0;

    rmed = L1ZERO(T,q,ks);
    rmin = rmed;
    rmax = rmed;

    rmax = L1OL2RAZ(rmax,10.0,T,q,ks,RAZ,+1);
	*mag_min = 1.0/DET(rmax,0.0,q,ks);

    for (int i=0;i<iT;i++)
    {
		rINF = L2ZERO(T,q,ks);
		rmed = L1ZERO(T,q,ks);
		rmin = L1OL2RAZ(rmed,rINF,T,q,ks,RAZ,-1);
		rmax = L1OL2RAZ(rmed,10.0,T,q,ks,RAZ,+1);
        r = rmin;
        dr = fabs(rmax - rmin)/double(iR);
        for(int i=0;i<iR;i++)
        {
            X = r*cos(T);
            Y = r*sin(T);
            DX = dr*cos(T);
            DY = dr*sin(T);

			Det = fabs(DET(X+DX/2.0,Y+DY/2.0,q,ks));
			dA_lente = r*dr*dT;

			if(1.0/Det > Mag_Min) AREA_lente += dA_lente;
            if(1.0/Det > Mag_Min) AREA_fonte += dA_lente*Det;

            r += dr;
        }
        T += dT;
    }
	*mag_med = AREA_lente/AREA_fonte;
    return 4.0*AREA_fonte;
}


