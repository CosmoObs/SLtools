#ifndef GRAFICOS_XMGRACE_H
#define GRAFICOS_XMGRACE_H

#include <grace_np.h>

#define PI 3.1415926535897932384626433832795 

//Aplicacoes com xmgrace
void my_error_function(const char *msg){fprintf(stderr, "library message: \"%s\"\n", msg);}
int caf_graf_xmgrace(double x0[], double y0[], double x1[], double y1[],int Npt,double q, double ks)
{
    GraceRegisterErrorFunction(my_error_function);

    /* Start Grace with a buffer size of 2048 and open the pipe */
    if (GraceOpen(2048) == -1) {
        fprintf(stderr, "Can't run Grace. \n");
        exit(EXIT_FAILURE);
    }
    
    /* Send some initialization commands to Grace */
    GracePrintf("world xmax 1.5");
    GracePrintf("world ymax 1.5");
    GracePrintf("world ymin -1.5");
    GracePrintf("world xmin -1.5");
    //GracePrintf("xaxis tick major 1");
    //GracePrintf("xaxis tick minor 0.5");
    //GracePrintf("yaxis tick major 0.2");
   // GracePrintf("yaxis tick minor 0.1");
    
    GracePrintf("page size 690, 690"); //tamanho da pagina
    GracePrintf("view 0.1, 0.1, 0.95, 0.95"); //regiao ocupada pelo grafico

    GracePrintf("yaxis label \"Y\\s2\""); //label do eixo y
    GracePrintf("xaxis label \"Y\\s1\"");  //label do eixo x
    //GracePrintf("title \"Teste\"");  //titulo do grafico
    GracePrintf("subtitle \"ks=%.1f, q=%.1f\"",ks,q);  //sibtitulo do grafico

    GracePrintf("legend on"); //ativa legenda
    GracePrintf("legend loctype view");
    GracePrintf("legend 0.1425, 0.225");
    GracePrintf("legend box color 1");
    GracePrintf("legend box pattern 1");
    GracePrintf("legend box linewidth 1.0");
    GracePrintf("legend box linestyle 1");
    GracePrintf("legend box fill color 0");
    GracePrintf("legend box fill pattern 1");
    GracePrintf("legend font 0");
    GracePrintf("legend char size 0.750000");
    GracePrintf("legend color 1");
    GracePrintf("legend length 3");
    GracePrintf("legend vgap 1");
    GracePrintf("legend hgap 1");
    GracePrintf("legend invert false");

    GracePrintf("s0 line type 1");
    GracePrintf("s0 line linestyle 1");
    GracePrintf("s0 line linewidth 1.5");
    GracePrintf("s0 line color 1");
    GracePrintf("s0 line pattern 1");

    GracePrintf("s0 on");
    GracePrintf("s0 symbol 0");
    GracePrintf("s0 symbol size 0.3");
    GracePrintf("s0 symbol fill pattern 1");
    GracePrintf("s0 legend  \" \\xl\\f{}\\sr \\N\\f{}= 0\"");

    GracePrintf("s1 line type 1");
    GracePrintf("s1 line linestyle 2");
    GracePrintf("s1 line linewidth 1.5");
    GracePrintf("s1 line color 1");
    GracePrintf("s1 line pattern 1");

    GracePrintf("s1 on");
    GracePrintf("s1 symbol 0");
    GracePrintf("s1 symbol size 0.3");
    GracePrintf("s1 symbol fill pattern 1");
    GracePrintf("s1 legend  \" \\xl\\f{}\\st \\N\\f{}= 0\"");

    /* Display sample data */
    for (int i = 0; i <= Npt ; i++){
    GracePrintf("g0.s0 point %E, %E", x0[i],y0[i]);
    GracePrintf("g0.s1 point %E, %E", x1[i],y1[i]);
    }

    /*Formatacao final*/
    //GracePrintf("autoscale");
    GracePrintf("autoticks");
    GracePrintf("redraw");

    /*Final do programa*/
    std::cout<<"Press Enter to exit";
    //std::cin.get();     //aguarda a tecla Enter

    /*Exporta para ps*/
    GracePrintf("HARDCOPY DEVICE \"png\""); 
    GracePrintf("print to \"caf-q%.1f-ks%.1f.ps\"",q,ks);
    GracePrintf("print");
    if (GraceIsOpen()) {
        /* Tell Grace to save the data */
        GracePrintf("saveall \"sample.agr\"");

        /* Flush the output buffer and close Grace */
        GraceClose();

        /* We are done */
        //exit(EXIT_SUCCESS);
    } else {
        //exit(EXIT_FAILURE);
    }
  return 0;
}

double CAF_xmgrace(double q, double ks, int I, int Npt)
{
  if(ks<=0.2 ) printf("\"L2ZERO\"->ks<=0.20:Usando apenas o metodo da bisecao\n");
  if(ks<=0.08) printf("\"L1ZERO\"->ks<=0.08:Usando apenas o metodo da bisecao\n");

  double Rsup,Rinf,T;
  double Xi[Npt],Yi[Npt];
  double Xs[Npt],Ys[Npt];
  double Xif[Npt],Yif[Npt];
  double Xsf[Npt],Ysf[Npt];

  for(int i=0;i<=Npt;i++)
  {
    T = 2.0*PI*i/double(Npt)/double(I);
    Rsup = L1ZERO(T,q,ks);
    Rinf = L2ZERO(T,q,ks);

    Xi[i] = Rinf*cos(T);
    Yi[i] = Rinf*sin(T);
    Xs[i] = Rsup*cos(T);
    Ys[i] = Rsup*sin(T);

    Xif[i] = XF(Xi[i],Yi[i],q,ks);
    Yif[i] = YF(Xi[i],Yi[i],q,ks);

    Xsf[i] = XF(Xs[i],Ys[i],q,ks);
    Ysf[i] = YF(Xs[i],Ys[i],q,ks);
    //	fprintf(saida2,"%E\t%E\n", YF(Xi,Yi,q,ks) , XF(Xi,Yi,q,ks) );
    }

  caf_graf_xmgrace(Xif,Yif,Xsf,Ysf,Npt,q,ks);

  return 0;
}


int ccl_graf_xmgrace(double x0[], double y0[], double x1[], double y1[],int Npt,double q, double ks)
{
    GraceRegisterErrorFunction(my_error_function);

    /* Start Grace with a buffer size of 2048 and open the pipe */
    if (GraceOpen(2048) == -1) {
        fprintf(stderr, "Can't run Grace. \n");
        exit(EXIT_FAILURE);
    }
    
    /* Send some initialization commands to Grace */
    GracePrintf("world xmax 2");
    GracePrintf("world ymax 2");
    GracePrintf("world ymin -2");
    GracePrintf("world xmin -2");
    //GracePrintf("xaxis tick major 1");
    //GracePrintf("xaxis tick minor 0.5");
    //GracePrintf("yaxis tick major 0.2");
   // GracePrintf("yaxis tick minor 0.1");
    
    GracePrintf("page size 690, 690"); //tamanho da pagina
    GracePrintf("view 0.1, 0.1, 0.95, 0.95"); //regiao ocupada pelo grafico

    GracePrintf("yaxis label \"X\\s2\""); //label do eixo y
    GracePrintf("xaxis label \"X\\s1\"");  //label do eixo x
    //GracePrintf("title \"Teste\"");  //titulo do grafico
    GracePrintf("subtitle \"ks=%.1f, q=%.1f\"",ks,q);  //sibtitulo do grafico

    GracePrintf("legend on"); //ativa legenda
    GracePrintf("legend loctype view");
    GracePrintf("legend 0.1425, 0.225");
    GracePrintf("legend box color 1");
    GracePrintf("legend box pattern 1");
    GracePrintf("legend box linewidth 1.0");
    GracePrintf("legend box linestyle 1");
    GracePrintf("legend box fill color 0");
    GracePrintf("legend box fill pattern 1");
    GracePrintf("legend font 0");
    GracePrintf("legend char size 0.750000");
    GracePrintf("legend color 1");
    GracePrintf("legend length 3");
    GracePrintf("legend vgap 1");
    GracePrintf("legend hgap 1");
    GracePrintf("legend invert false");

    GracePrintf("s0 line type 1");
    GracePrintf("s0 line linestyle 1");
    GracePrintf("s0 line linewidth 1.5");
    GracePrintf("s0 line color 1");
    GracePrintf("s0 line pattern 1");

    GracePrintf("s0 on");
    GracePrintf("s0 symbol 0");
    GracePrintf("s0 symbol size 0.3");
    GracePrintf("s0 symbol fill pattern 1");
    GracePrintf("s0 legend  \" \\xl\\f{}\\sr \\N\\f{}= 0\"");

    GracePrintf("s1 line type 1");
    GracePrintf("s1 line linestyle 2");
    GracePrintf("s1 line linewidth 1.5");
    GracePrintf("s1 line color 1");
    GracePrintf("s1 line pattern 1");

    GracePrintf("s1 on");
    GracePrintf("s1 symbol 0");
    GracePrintf("s1 symbol size 0.3");
    GracePrintf("s1 symbol fill pattern 1");
    GracePrintf("s1 legend  \" \\xl\\f{}\\st \\N\\f{}= 0\"");

    /* Display sample data */
    for (int i = 0; i <= Npt ; i++){
    GracePrintf("g0.s0 point %E, %E", x0[i],y0[i]);
    GracePrintf("g0.s1 point %E, %E", x1[i],y1[i]);
    }

    /*Formatacao final*/
    //GracePrintf("autoscale");
    GracePrintf("autoticks");
    GracePrintf("redraw");

    /*Final do programa*/
    std::cout<<"Press Enter to exit";
    //std::cin.get();     //aguarda a tecla Enter

    /*Exporta para ps*/
    //GracePrintf("HARDCOPY DEVICE \"png\""); 
    GracePrintf("print to \"ccl-q%.1f-ks%.1f.ps\"",q,ks);
    GracePrintf("print");
    if (GraceIsOpen()) {
        /* Tell Grace to save the data */
        GracePrintf("saveall \"sample.agr\"");

        /* Flush the output buffer and close Grace */
        GraceClose();

        /* We are done */
        //exit(EXIT_SUCCESS);
    } else {
        //exit(EXIT_FAILURE);
    }
  return 0;
}

double CCL_xmgrace(double q, double ks, int I, int Npt)
{
  if(ks<=0.2 ) printf("\"L2ZERO\"->ks<=0.20:Usando apenas o metodo da bisecao\n");
  if(ks<=0.08) printf("\"L1ZERO\"->ks<=0.08:Usando apenas o metodo da bisecao\n");

  double Rsup,Rinf,T;
  double Xi[Npt+1],Yi[Npt+1];
  double Xs[Npt+1],Ys[Npt+1];
  //double Xif[Npt],Yif[Npt];
  //double Xsf[Npt],Ysf[Npt];

  for(int i=0;i<=Npt;i++)
  {
    T = 2.0*PI*i/double(Npt)/double(I);
    Rsup = L1ZERO(T,q,ks);
    Rinf = L2ZERO(T,q,ks);

    Xi[i] = Rinf*cos(T);
    Yi[i] = Rinf*sin(T);
    Xs[i] = Rsup*cos(T);
    Ys[i] = Rsup*sin(T);

    //Xif[i] = XF(Xi[i],Yi[i],q,ks);
    //Yif[i] = YF(Xi[i],Yi[i],q,ks);

    //Xsf[i] = XF(Xs[i],Ys[i],q,ks);
    //Ysf[i] = YF(Xs[i],Ys[i],q,ks);
    //	fprintf(saida2,"%E\t%E\n", YF(Xi,Yi,q,ks) , XF(Xi,Yi,q,ks) );
    }

  ccl_graf_xmgrace(Xi,Yi,Xs,Ys,Npt,q,ks);

  return 0;
}


int graficosf_graf_xmgrace(double x0[], double y0[], double x1[], double y1[],double x2[], double y2[], int Npt,double q, double ks, double RAZ)
{
    GraceRegisterErrorFunction(my_error_function);

    /* Start Grace with a buffer size of 2048 and open the pipe */
    if (GraceOpen(2048) == -1) {
        fprintf(stderr, "Can't run Grace. \n");
        exit(EXIT_FAILURE);
    }
    
    /* Send some initialization commands to Grace */
    GracePrintf("world xmax 0.7");
    GracePrintf("world ymax 0.7");
    GracePrintf("world ymin -0.7");
    GracePrintf("world xmin -0.7");
    //GracePrintf("xaxis tick major 1");
    //GracePrintf("xaxis tick minor 0.5");
    //GracePrintf("yaxis tick major 0.2");
   // GracePrintf("yaxis tick minor 0.1");
    
    GracePrintf("page size 690, 690"); //tamanho da pagina
    GracePrintf("view 0.1, 0.1, 0.95, 0.95"); //regiao ocupada pelo grafico

    GracePrintf("yaxis label \"X\\s2\""); //label do eixo y
    GracePrintf("xaxis label \"X\\s1\"");  //label do eixo x
    //GracePrintf("title \"Teste\"");  //titulo do grafico
    GracePrintf("subtitle \"ks=%.1f, q=%.1f\"",ks,q);  //sibtitulo do grafico

    GracePrintf("legend on"); //ativa legenda
    GracePrintf("legend loctype view");
    GracePrintf("legend 0.1425, 0.225");
    GracePrintf("legend box color 1");
    GracePrintf("legend box pattern 1");
    GracePrintf("legend box linewidth 1.0");
    GracePrintf("legend box linestyle 1");
    GracePrintf("legend box fill color 0");
    GracePrintf("legend box fill pattern 1");
    GracePrintf("legend font 0");
    GracePrintf("legend char size 0.750000");
    GracePrintf("legend color 1");
    GracePrintf("legend length 3");
    GracePrintf("legend vgap 1");
    GracePrintf("legend hgap 1");
    GracePrintf("legend invert false");

    GracePrintf("s0 line type 1");
    GracePrintf("s0 line linestyle 1");
    GracePrintf("s0 line linewidth 1.5");
    GracePrintf("s0 line color 1");
    GracePrintf("s0 line pattern 1");

    GracePrintf("s0 on");
    GracePrintf("s0 symbol 0");
    GracePrintf("s0 symbol size 0.3");
    GracePrintf("s0 symbol fill pattern 1");
    GracePrintf("s0 legend  \" \\xl\\f{}\\sr \\N\\f{}= 0\"");

    GracePrintf("s1 line type 1");
    GracePrintf("s1 line linestyle 3");
    GracePrintf("s1 line linewidth 1.5");
    GracePrintf("s1 line color 1");
    GracePrintf("s1 line pattern 1");

    GracePrintf("s1 on");
    GracePrintf("s1 symbol 0");
    GracePrintf("s1 symbol size 0.3");
    GracePrintf("s1 symbol fill pattern 1");
    //GracePrintf("s1 legend  \" C\\\\W = +%.0f\"",RAZ);

    GracePrintf("s2 line type 1");
    GracePrintf("s2 line linestyle 3");
    GracePrintf("s2 line linewidth 1.5");
    GracePrintf("s2 line color 1");
    GracePrintf("s2 line pattern 1");

    GracePrintf("s2 on");
    GracePrintf("s2 symbol 0");
    GracePrintf("s2 symbol size 0.3");
    GracePrintf("s2 symbol fill pattern 1");
    GracePrintf("s2 legend  \" |C\\\\W| = %.0f\"",RAZ);

    /* Display sample data */
    for (int i = 0; i <= Npt ; i++){
    GracePrintf("g0.s0 point %E, %E", x0[i],y0[i]);
    GracePrintf("g0.s1 point %E, %E", x1[i],y1[i]);
    GracePrintf("g0.s2 point %E, %E", x2[i],y2[i]);
    }

    /*Formatacao final*/
    //GracePrintf("autoscale");
    GracePrintf("autoticks");
    GracePrintf("redraw");

    /*Final do programa*/
    std::cout<<"Press Enter to exit";
    std::cin.get();     //aguarda a tecla Enter

    /*Exporta para ps*/
    //GracePrintf("HARDCOPY DEVICE \"png\""); 
    GracePrintf("print to \"graf-F-q%.1f-ks%.1f.ps\"",q,ks);
    GracePrintf("print");
    if (GraceIsOpen()) {
        /* Tell Grace to save the data */
        GracePrintf("saveall \"sample.agr\"");

        /* Flush the output buffer and close Grace */
        GraceClose();

        /* We are done */
        //exit(EXIT_SUCCESS);
    } else {
        //exit(EXIT_FAILURE);
    }
  return 0;
}

double graficosF_xmgrace(double q, double ks, double RAZ, int I, int Npt)
{

	double R,Rsup,Rinf;
	double rINF;
	double X[Npt],Y[Npt+1],Xs[Npt+1],Ys[Npt+1],Xi[Npt+1],Yi[Npt+1];
	double Xf[Npt],Yf[Npt+1],Xsf[Npt+1],Ysf[Npt+1],Xif[Npt+1],Yif[Npt+1];
	double T;

   // FILE *saida1 = fopen(file1,"w+");
   // FILE *saida2 = fopen(file2,"w+");
   // FILE *saida3 = fopen(file3,"w+");

    for(int i=0;i<=Npt;i++)
    {
		T = 2.0*PI*i/double(Npt)/double(I);
		rINF = L2ZERO(T,q,ks);
    	R = L1ZERO(T,q,ks);
    	Rsup = L1OL2RAZ(R,10.0,T,q,ks,RAZ,+1);
    	Rinf = L1OL2RAZ(R,rINF,T,q,ks,RAZ,-1);
		X[i] = R*cos(T);
		Y[i] = R*sin(T);
		Xs[i]= Rsup*cos(T);
		Ys[i]= Rsup*sin(T);
		Xi[i]= Rinf*cos(T);
		Yi[i]= Rinf*sin(T);
		
		Xf[i] = XF(X[i],Y[i],q,ks);
		Yf[i] = YF(X[i],Y[i],q,ks);
		Xif[i]= XF(Xi[i],Yi[i],q,ks);
		Yif[i]= YF(Xi[i],Yi[i],q,ks);
		Xsf[i]= XF(Xs[i],Ys[i],q,ks);
		Ysf[i]= YF(Xs[i],Ys[i],q,ks);
    	//fprintf(saida1,"%E\t%E\n", XF(X,Y,q,ks),YF(X,Y,q,ks));
    	//fprintf(saida2,"%E\t%E\n", XF(Xi,Yi,q,ks),YF(Xi,Yi,q,ks));
    	//fprintf(saida3,"%E\t%E\n", XF(Xs,Ys,q,ks),YF(Xs,Ys,q,ks));
    }

  graficosf_graf_xmgrace(Xf,Yf,Xsf,Ysf,Xif,Yif,Npt,q,ks,RAZ);
    //fclose(saida1);
    //fclose(saida2);
    //fclose(saida3);

	return 0.0;
}



int graficosl_graf_xmgrace(double x0[], double y0[], double x1[], double y1[],double x2[], double y2[], int Npt,double q, double ks, double RAZ)
{
    GraceRegisterErrorFunction(my_error_function);

    /* Start Grace with a buffer size of 2048 and open the pipe */
    if (GraceOpen(2048) == -1) {
        fprintf(stderr, "Can't run Grace. \n");
        exit(EXIT_FAILURE);
    }
    
    /* Send some initialization commands to Grace */
    GracePrintf("world xmax 2");
    GracePrintf("world ymax 2");
    GracePrintf("world ymin -2");
    GracePrintf("world xmin -2");
    //GracePrintf("xaxis tick major 1");
    //GracePrintf("xaxis tick minor 0.5");
    //GracePrintf("yaxis tick major 0.2");
   // GracePrintf("yaxis tick minor 0.1");
    
    GracePrintf("page size 690, 690"); //tamanho da pagina
    GracePrintf("view 0.1, 0.1, 0.95, 0.95"); //regiao ocupada pelo grafico

    GracePrintf("yaxis label \"X\\s2\""); //label do eixo y
    GracePrintf("xaxis label \"X\\s1\"");  //label do eixo x
    //GracePrintf("title \"Teste\"");  //titulo do grafico
    GracePrintf("subtitle \"ks=%.1f, q=%.1f\"",ks,q);  //sibtitulo do grafico

    GracePrintf("legend on"); //ativa legenda
    GracePrintf("legend loctype view");
    GracePrintf("legend 0.1425, 0.225");
    GracePrintf("legend box color 1");
    GracePrintf("legend box pattern 1");
    GracePrintf("legend box linewidth 1.0");
    GracePrintf("legend box linestyle 1");
    GracePrintf("legend box fill color 0");
    GracePrintf("legend box fill pattern 1");
    GracePrintf("legend font 0");
    GracePrintf("legend char size 0.750000");
    GracePrintf("legend color 1");
    GracePrintf("legend length 3");
    GracePrintf("legend vgap 1");
    GracePrintf("legend hgap 1");
    GracePrintf("legend invert false");

    GracePrintf("s0 line type 1");
    GracePrintf("s0 line linestyle 1");
    GracePrintf("s0 line linewidth 1.5");
    GracePrintf("s0 line color 1");
    GracePrintf("s0 line pattern 1");

    GracePrintf("s0 on");
    GracePrintf("s0 symbol 0");
    GracePrintf("s0 symbol size 0.3");
    GracePrintf("s0 symbol fill pattern 1");
    GracePrintf("s0 legend  \" \\xl\\f{}\\sr \\N\\f{}= 0\"");

    GracePrintf("s1 line type 1");
    GracePrintf("s1 line linestyle 2");
    GracePrintf("s1 line linewidth 1.5");
    GracePrintf("s1 line color 1");
    GracePrintf("s1 line pattern 1");

    GracePrintf("s1 on");
    GracePrintf("s1 symbol 0");
    GracePrintf("s1 symbol size 0.3");
    GracePrintf("s1 symbol fill pattern 1");
    //GracePrintf("s1 legend  \" C\\\\W = +%.0f\"",RAZ);

    GracePrintf("s2 line type 1");
    GracePrintf("s2 line linestyle 2");
    GracePrintf("s2 line linewidth 1.5");
    GracePrintf("s2 line color 1");
    GracePrintf("s2 line pattern 1");

    GracePrintf("s2 on");
    GracePrintf("s2 symbol 0");
    GracePrintf("s2 symbol size 0.3");
    GracePrintf("s2 symbol fill pattern 1");
    GracePrintf("s2 legend  \" |C\\\\W| = %.0f\"",RAZ);

    /* Display sample data */
    for (int i = 0; i <= Npt ; i++){
    GracePrintf("g0.s0 point %E, %E", x0[i],y0[i]);
    GracePrintf("g0.s1 point %E, %E", x1[i],y1[i]);
    GracePrintf("g0.s2 point %E, %E", x2[i],y2[i]);
    }

    /*Formatacao final*/
    //GracePrintf("autoscale");
    GracePrintf("autoticks");
    GracePrintf("redraw");

    /*Final do programa*/
    std::cout<<"Press Enter to exit";
    std::cin.get();     //aguarda a tecla Enter

    /*Exporta para ps*/
    //GracePrintf("HARDCOPY DEVICE \"png\""); 
    GracePrintf("print to \"graf-L-q%.1f-ks%.1f.ps\"",q,ks);
    GracePrintf("print");
    if (GraceIsOpen()) {
        /* Tell Grace to save the data */
        GracePrintf("saveall \"sample.agr\"");

        /* Flush the output buffer and close Grace */
        GraceClose();

        /* We are done */
        //exit(EXIT_SUCCESS);
    } else {
        //exit(EXIT_FAILURE);
    }
  return 0;
}

double graficosL_xmgrace(double q, double ks, double RAZ, int I, int Npt)
{

	double R,Rsup,Rinf;
	double rINF;
	double X[Npt],Y[Npt+1],Xs[Npt+1],Ys[Npt+1],Xi[Npt+1],Yi[Npt+1];
	double T;

   // FILE *saida1 = fopen(file1,"w+");
   // FILE *saida2 = fopen(file2,"w+");
   // FILE *saida3 = fopen(file3,"w+");

    for(int i=0;i<=Npt;i++)
    {
		T = 2.0*PI*i/double(Npt)/double(I);
		rINF = L2ZERO(T,q,ks);
    	R = L1ZERO(T,q,ks);
    	Rsup = L1OL2RAZ(R,10.0,T,q,ks,RAZ,+1);
    	Rinf = L1OL2RAZ(R,rINF,T,q,ks,RAZ,-1);
		X[i] = R*cos(T);
		Y[i] = R*sin(T);
		Xs[i]= Rsup*cos(T);
		Ys[i]= Rsup*sin(T);
		Xi[i]= Rinf*cos(T);
		Yi[i]= Rinf*sin(T);
    	//fprintf(saida1,"%E\t%E\n", XF(X,Y,q,ks),YF(X,Y,q,ks));
    	//fprintf(saida2,"%E\t%E\n", XF(Xi,Yi,q,ks),YF(Xi,Yi,q,ks));
    	//fprintf(saida3,"%E\t%E\n", XF(Xs,Ys,q,ks),YF(Xs,Ys,q,ks));
    }

  graficosl_graf_xmgrace(X,Y,Xs,Ys,Xi,Yi,Npt,q,ks,RAZ);
    //fclose(saida1);
    //fclose(saida2);
    //fclose(saida3);

	return 0.0;
}
#endif
