/** @file
* Example of doxygen documentation for C functions FIXME. 
*/

/** @package curvas
*  Package.
*
*  Detailed descrition FIXME
*
*/

#ifndef CURVAS_H
#define CURVAS_H

#define PI 3.1415926535897932384626433832795 


double NRrtbis(double (*func)(double, double, double, double), double r1, double r2, double xacc, double q, double ks, double T, double* r_sup, double* r_inf)
{
	const int JMAX=100;
	int j;
	double dx,f,fmid,xmid,rtb;
	double COST,SINT;

	COST = cos(T);
	SINT = sin(T);

	f=func(r1*COST,r1*SINT,q,ks);
	fmid=func(r2*COST,r2*SINT,q,ks);

	if (f*fmid >= 0.0 && f+fmid<100.0) printf("NR->Root must be bracketed for bisection in rtbis\n");

	rtb = f < 0.0 ? (dx=r2-r1,r1) : (dx=r1-r2,r2);
	for (j=0;j<JMAX;j++) {
		xmid=rtb+(dx *= 0.5);
		fmid=func(xmid*COST,xmid*SINT,q,ks);
		if (fmid <= 0.0) rtb=xmid;
		if (fabs(dx) < xacc || fmid == 0.0){
		*r_sup = xmid+dx;
		*r_inf = xmid-dx;
		//printf("%f\t%f\t%f\n",rtb,xmid+dx,xmid-dx);
		return rtb;
		}
	}
	printf("NR->Too many bisections in rtbis\n");
	return 0.0;
}
double NRrtbis(double (*func)(double, double, double, double), double r1, double r2, double xacc, double q, double ks, double T, double* r_sup, double* r_inf, double RAZ)
{
	const int JMAX=100;
	int j;
	double dx,f,fmid,xmid,rtb;
	double COST,SINT;

	COST = cos(T);
	SINT = sin(T);

	f=func(r1*COST,r1*SINT,q,ks) - RAZ;
	fmid=func(r2*COST,r2*SINT,q,ks) - RAZ;

	if (f*fmid >= 0.0 && f+fmid<100.0) printf("NR->Root must be bracketed for bisection in rtbis\n");

	rtb = f < 0.0 ? (dx=r2-r1,r1) : (dx=r1-r2,r2);
	for (j=0;j<JMAX;j++) {
		xmid=rtb+(dx *= 0.5);
		fmid=func(xmid*COST,xmid*SINT,q,ks) - RAZ;
		if (fmid <= 0.0) rtb=xmid;
		if (fabs(dx) < xacc || fmid == 0.0){
		*r_sup = xmid+dx;
		*r_inf = xmid-dx;
		//printf("%f\t%f\t%f\n",rtb,xmid+dx,xmid-dx);
		return rtb;
		}
	}
	printf("NR->Too many bisections in rtbis\n");
	return 0.0;
}



double NRrtflsp(double (*func)(double, double, double, double), double r1, double r2, double xacc, double q, double ks, double T)
{
	const int MAXIT=1000;
	int j;
	double fl,fh,xl,xh,swap,dx,del,f,rtf;
	double COST,SINT;
	
	COST = cos(T);
	SINT = sin(T);

	fl=func(r1*COST,r1*SINT,q,ks);
	fh=func(r2*COST,r2*SINT,q,ks);
	if (fl*fh > 0.0) printf("NR->Root must be bracketed in rtflsp\n");
	if (fl < 0.0) {
		xl=r1;
		xh=r2;
	} else {
		xl=r2;
		xh=r1;
		swap=fl;
		fl=fh;
		fh=swap;
	}
	dx=xh-xl;
	for (j=0;j<MAXIT;j++) {
		rtf=xl+dx*fl/(fl-fh);
		f=func(rtf*COST,rtf*SINT,q,ks);
		if (f < 0.0) {
			del=xl-rtf;
			xl=rtf;
			fl=f;
		} else {
			del=xh-rtf;
			xh=rtf;
			fh=f;
		}
		dx=xh-xl;
		//printf("%i\n",j);
		if (fabs(del) < xacc || f == 0.0) return rtf;
	}
	printf("NR->Maximum number of iterations exceeded in rtflsp\n");
	return 0.0;
}
double NRrtflsp(double (*func)(double, double, double, double), double r1, double r2, double xacc, double q, double ks, double T, double RAZ)
{
	const int MAXIT=1000;
	int j;
	double fl,fh,xl,xh,swap,dx,del,f,rtf;
	double COST,SINT;
	
	COST = cos(T);
	SINT = sin(T);

	fl=func(r1*COST,r1*SINT,q,ks) - RAZ;
	fh=func(r2*COST,r2*SINT,q,ks) - RAZ;
	if (fl*fh > 0.0) printf("NR->Root must be bracketed in rtflsp\n");
	if (fl < 0.0) {
		xl=r1;
		xh=r2;
	} else {
		xl=r2;
		xh=r1;
		swap=fl;
		fl=fh;
		fh=swap;
	}
	dx=xh-xl;
	for (j=0;j<MAXIT;j++) {
		rtf=xl+dx*fl/(fl-fh);
		f=func(rtf*COST,rtf*SINT,q,ks) - RAZ;
		if (f < 0.0) {
			del=xl-rtf;
			xl=rtf;
			fl=f;
		} else {
			del=xh-rtf;
			xh=rtf;
			fh=f;
		}
		dx=xh-xl;
		//printf("%i\n",j);
		if (fabs(del) < xacc || f == 0.0) return rtf;
	}
	printf("NR->Maximum number of iterations exceeded in rtflsp\n");
	return 0.0;
}




/*---ACHANDO OS ZEROS---*******************************************************/
double L1ZERO(double T, double q, double ks)
{
	double *r_sup = new double;
	double *r_inf = new double;

	if(ks<=0.08) return NRrtbis(L1, 1.0, 0.0, 1E-8, q, ks, T, r_sup, r_inf);

	NRrtbis(L1, 10.0, 0.0001, 1E-1, q, ks, T, r_sup, r_inf);
	return NRrtflsp(L1, *r_sup, *r_inf, 1E-8, q, ks, T);
}
/*---ACHANDO OS ZEROS---*******************************************************/
double L2ZERO(double T, double q, double ks)
{

	double *r_sup = new double;
	double *r_inf = new double;

	if(ks<=0.2) return NRrtbis(L2, 1.0, 0.0, 1E-8, q, ks, T, r_sup, r_inf);

	NRrtbis(L2, 10.0, 0.0001, 1E-1, q, ks, T, r_sup, r_inf);
    return NRrtflsp(L2, *r_sup, *r_inf, 1E-8, q, ks, T);
}

/*---ACHANDO OS ZEROS DAS RAZOES DOS AUTOVALORES---*****************************************************************************/
double L2OL1RAZ(double rIN, double T, double q, double ks, double RAZ,int S)
{
    int i;
    double F,F1;
    double r,dr;
    double X,Y;

    
    F1 = 0.0;
    r = rIN;
    dr = S;
    
       for (i=0;i<=5000;i++)
       {
            X = r*cos(T);
            Y = r*sin(T);
            F = fabs(L2OL1(X,Y,q,ks)) - RAZ;
            if(F < 0) r = r + dr;
            if(F > 0){  r = r - dr;  dr = dr/5.0;  }
            if(fabs(F1-F) < 1E-7) break;
            F1 = F;
       }
   return r + 5.0*dr;
}

/******************************************************************************/
double L1OL2RAZ(double rMED, double rINF, double T, double q, double ks, double RAZ,int S)
{

	double *r_sup = new double;
	double *r_inf = new double;
	
	if(S<0) return NRrtbis(L1OL2, rMED, rINF, 1E-8, q, ks, T, r_sup, r_inf, -1.0/RAZ);

	if(S>0) return NRrtbis(L1OL2, rMED, 10.0, 1E-8, q, ks, T, r_sup, r_inf,  1.0/RAZ);

	printf("\"L1OL2RAZ\"-> S=0");
	return 0.0;
}

/*---GRAFICOS---***************************************************************/
double CAF(double q, double ks, int I, int Npt, char file1[], char file2[])
{
	if(ks<=0.2 ) printf("\"L2ZERO\"->ks<=0.20:Usando apenas o metodo da bisecao\n");
	if(ks<=0.08) printf("\"L1ZERO\"->ks<=0.08:Usando apenas o metodo da bisecao\n");

	double Rsup,Rinf,T;
	double Xi,Yi;
	double Xs,Ys;

    FILE *saida1 = fopen(file1,"w+");
    FILE *saida2 = fopen(file2,"w+");

    for(int i=0;i<=Npt;i++)
    {
		T = 2.0*PI*i/double(Npt)/double(I);
    	Rsup = L1ZERO(T,q,ks);
    	Rinf = L2ZERO(T,q,ks);

		Xi = Rinf*cos(T);
		Yi = Rinf*sin(T);
		Xs = Rsup*cos(T);
		Ys = Rsup*sin(T);

    	fprintf(saida1,"%E\t%E\n", YF(Xs,Ys,q,ks) , XF(Xs,Ys,q,ks) );
    	fprintf(saida2,"%E\t%E\n", YF(Xi,Yi,q,ks) , XF(Xi,Yi,q,ks) );
    }

    fclose(saida1);
    fclose(saida2);

	return 0;
}

double CCL(double q, double ks, int I, int Npt, char file1[], char file2[])
{
	if(ks<=0.2 ) printf("\"L2ZERO\"->ks<=0.20:Usando apenas o metodo da bisecao\n");
	if(ks<=0.08) printf("\"L1ZERO\"->ks<=0.08:Usando apenas o metodo da bisecao\n");

	double Rsup,Rinf,T;

    FILE *saida1 = fopen(file1,"w+");
    FILE *saida2 = fopen(file2,"w+");

    for(int i=0;i<=Npt;i++)
    {
		T = 2.0*PI*i/double(Npt)/double(I);

    	Rsup = L1ZERO(T,q,ks);
    	Rinf = L2ZERO(T,q,ks);
    	fprintf(saida1,"%E\t%E\n", Rsup*cos(T), Rsup*sin(T) );
    	fprintf(saida2,"%E\t%E\n", Rinf*cos(T), Rinf*sin(T) );
    }

    fclose(saida1);
    fclose(saida2);

	return 0;
}

double graficosL(double q, double ks, double RAZ, int I, int Npt, char file1[], char file2[], char file3[])
{

	double R,Rsup,Rinf;
	double rINF;
	double T;

    FILE *saida1 = fopen(file1,"w+");
    FILE *saida2 = fopen(file2,"w+");
    FILE *saida3 = fopen(file3,"w+");


    for(int i=0;i<=Npt;i++)
    {
	T = 2.0*PI*i/double(Npt)/double(I);
	rINF = L2ZERO(T,q,ks);

    	R = L1ZERO(T,q,ks);
    	Rsup = L1OL2RAZ(R,10.0,T,q,ks,RAZ,+1);
    	Rinf = L1OL2RAZ(R,rINF,T,q,ks,RAZ,-1);
    	fprintf(saida1,"%E\t%E\n", R*cos(T),R*sin(T));
    	fprintf(saida2,"%E\t%E\n", Rinf*cos(T),Rinf*sin(T) );
    	fprintf(saida3,"%E\t%E\n", Rsup*cos(T),Rsup*sin(T) );
    }

    fclose(saida1);
    fclose(saida2);
    fclose(saida3);

	return 0.0;
}


double graficosF(double q, double ks, double RAZ, int I, int Npt, char file1[], char file2[], char file3[])
{

	double R,Rsup,Rinf;
	double rINF;
	double X,Y,Xs,Ys,Xi,Yi;
	double T;

    FILE *saida1 = fopen(file1,"w+");
    FILE *saida2 = fopen(file2,"w+");
    FILE *saida3 = fopen(file3,"w+");

    for(int i=0;i<=Npt;i++)
    {
		T = 2.0*PI*i/double(Npt)/double(I);
		rINF = L2ZERO(T,q,ks);
    	R = L1ZERO(T,q,ks);
    	Rsup = L1OL2RAZ(R,10.0,T,q,ks,RAZ,+1);
    	Rinf = L1OL2RAZ(R,rINF,T,q,ks,RAZ,-1);
		X = R*cos(T);
		Y = R*sin(T);
		Xs= Rsup*cos(T);
		Ys= Rsup*sin(T);
		Xi= Rinf*cos(T);
		Yi= Rinf*sin(T);
    	fprintf(saida1,"%E\t%E\n", XF(X,Y,q,ks),YF(X,Y,q,ks));
    	fprintf(saida2,"%E\t%E\n", XF(Xi,Yi,q,ks),YF(Xi,Yi,q,ks));
    	fprintf(saida3,"%E\t%E\n", XF(Xs,Ys,q,ks),YF(Xs,Ys,q,ks));
    }

    fclose(saida1);
    fclose(saida2);
    fclose(saida3);

	return 0.0;
}





#endif


/******************************************************************************/
/*double L1OL2RAZ(double rIN, double T, double q, double ks, double RAZ,int S)
{
    int i;
    double F,F1;
    double r,dr;
    double COST,SINT;
    double X,Y;
    
    F1 = 0.0;
    dr = double(S)/10.0;
    r = rIN;// + S*dr;
    COST = cos(T);
    SINT = sin(T);
    
    
       for (i=0;i<=5000;i++)
       {
            X = r*COST;
            Y = r*SINT;
            F = fabs(L1OL2(X,Y,q,ks));// - RAZ;
            if(F < RAZ){ r = r + dr; }
            if(F > RAZ){  r = r - dr;  dr = dr/5.0;}
            if(fabs(F1-F) < 1.0E-12) break;
            F1 = F;
       }
    return r + 5.0*dr;
}*/
/******************************************************************************/
