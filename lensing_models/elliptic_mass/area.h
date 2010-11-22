#ifndef AREA_H
#define AREA_H


#define PI 3.1415926535897932384626433832795 


double AREA_(double q, double ks, double RAZ, double T, double* Intl)
{//Calcula a parte radial da itegral da area

    double dr,r;
    double rmin,rmax,rmed;
	double rINF;
    double X,Y,DX,DY;
	double N=50;    
    double h;
    double INTF,INTL;
    double p0,p1,p2;

    r = 0.0;
    DX = dr*cos(T);
    DY = dr*sin(T); 

    rmed = L1ZERO(T,q,ks);
	rINF = L2ZERO(T,q,ks);
    rmin = rmed;
    rmax = rmed;
    rmin = L1OL2RAZ(rmed,rINF,T,q,ks,RAZ,-1);
    rmax = L1OL2RAZ(rmed,10.0,T,q,ks,RAZ,+1);
    
    h = (rmax-rmin)/double(N);
    r = rmin;
    INTF = 0.0;
	INTL = 0.0;
    
    p0 = 0.3750000000000;
    p1 = 1.1666666666667;
    p2 = 0.9583333333333;
    

    X = r*cos(T);
    Y = r*sin(T);
    INTF += p0*r*fabs(DET(X,Y,q,ks));
	INTL += p0*r;
    r += h;

    X = r*cos(T);
    Y = r*sin(T);
    INTF += p1*r*fabs(DET(X,Y,q,ks));
	INTL += p1*r;
    r += h;

    X = r*cos(T);
    Y = r*sin(T);
    INTF += p2*r*fabs(DET(X,Y,q,ks));
	INTL += p2*r;
    r += h;
    
    for(int i=3;i<=N-3;i++)
    {
 	    X = r*cos(T);
 	    Y = r*sin(T);
        INTF += r*fabs(DET(X,Y,q,ks));
		INTL += r;
        r += h;
    }

    X = r*cos(T);
    Y = r*sin(T);
    INTF += p2*r*fabs(DET(X,Y,q,ks));
	INTL += p2*r;
    r += h;

    X = r*cos(T);
    Y = r*sin(T);
    INTF += p1*r*fabs(DET(X,Y,q,ks));
	INTL += p1*r;
    r += h;

    X = r*cos(T);
    Y = r*sin(T);
    INTF += p0*r*fabs(DET(X,Y,q,ks));
    INTL += p0*r;

    *Intl = h*INTL;

    return h*INTF;   
}


double AREA_2(double q, double ks, double RAZ,double* mag_min, double* mag_med)
{//Calcula a parte angular da itegral da area

    double INTF=0.0;
	double INTL=0.0;
	double* Intl = new double;
	double rmed,rmax;
    
    int D=25;
    double x[D+1],w[D+1];


x[0]=0.0;                   w[0]=0.0;
x[1]=3.489547766462153E-03;	w[1]=8.948668416825708E-03;
x[2]=1.832811322659356E-02;	w[2]=2.069915808381061E-02;
x[3]=4.478766702371562E-02;	w[3]=3.215353848424626E-02;
x[4]=8.247309200843789E-02;	w[4]=4.312204738131632E-02;
x[5]=1.308138534715479E-01;	w[5]=5.343718241684769E-02;
x[6]=1.890773326654075E-01;	w[6]=6.294235885649400E-02;
x[7]=2.563803746720973E-01;	w[7]=7.149342977868132E-02;
x[8]=3.317027589232035E-01;	w[8]=7.896074975268089E-02;
x[9]=4.139026868380073E-01;	w[9]=8.523111375466741E-02;
x[10]=5.017340977447335E-01;	w[10]=9.020946578407005E-02;
x[11]=5.938655608095249E-01;	w[11]=9.382033728824524E-02;
x[12]=6.889004594746043E-01;	w[12]=9.600899021380682E-02;
x[13]=7.853981633974483E-01;	w[13]=9.674224637150773E-02;
x[14]=8.818958673202922E-01;	w[14]=9.600899021380682E-02;
x[15]=9.769307659853717E-01;	w[15]=9.382033728824524E-02;
x[16]=1.069062229050163E+00;	w[16]=9.020946578407005E-02;
x[17]=1.156893639956889E+00;	w[17]=8.523111375466741E-02;
x[18]=1.239093567871693E+00;	w[18]=7.896074975268089E-02;
x[19]=1.314415952122799E+00;	w[19]=7.149342977868132E-02;
x[20]=1.381718994129489E+00;	w[20]=6.294235885649400E-02;
x[21]=1.439982473323349E+00;	w[21]=5.343718241684769E-02;
x[22]=1.488323234786459E+00;	w[22]=4.312204738131632E-02;
x[23]=1.526008659771181E+00;	w[23]=3.215353848424626E-02;
x[24]=1.552468213568303E+00;	w[24]=2.069915808381061E-02;
x[25]=1.567306779028435E+00;	w[25]=8.948668416825708E-03;


//Calculo da magnificacao minima
	rmed = L1ZERO(0.0,q,ks);
    rmax = L1OL2RAZ(rmed,10.0,0.0,q,ks,RAZ,+1);
	*mag_min = 1.0/DET(rmax,0.0,q,ks);
//Fim do calculo da magnificacao minima


//Calculo da parte angular da integral da area
    for(int i=1;i<=D;i++)
    {
    INTF += w[i]*AREA_(q,ks,RAZ,x[i],Intl);
	INTL += w[i]*(*Intl);
    }
//Fim de calculo da parte angular da integral da area

    *mag_med = INTL/INTF;
    return 4.0*INTF;
}


/******************************************************************************/


double AREA_(double q, double ks, double RAZ, double T, double* Intl, double MAG_MIN)
{//Calcula a parte radial da itegral da area

    double dr,r;
    double rmin,rmax,rmed;
	double rINF;
    double X,Y,DX,DY;
	double N=50;    
    double h;
    double INTF,INTL;
    double p0,p1,p2;

    r = 0.0;
    DX = dr*cos(T);
    DY = dr*sin(T); 

    rmed = L1ZERO(T,q,ks);
	rINF = L2ZERO(T,q,ks);
    rmin = rmed;
    rmax = rmed;
    rmin = L1OL2RAZ(rmed,rINF,T,q,ks,RAZ,-1);
    rmax = L1OL2RAZ(rmed,10.0,T,q,ks,RAZ,+1);

    
    h = (rmax-rmin)/double(N);
    r = rmin;
    INTF = 0.0;
	INTL = 0.0;
    
    p0 = 0.3750000000000;
    p1 = 1.1666666666667;
    p2 = 0.9583333333333;
    

    X = r*cos(T);
    Y = r*sin(T);
    INTF += p0*r*fabs(DET(X,Y,q,ks));
	INTL += p0*r;
    r += h;

    X = r*cos(T);
    Y = r*sin(T);
    INTF += p1*r*fabs(DET(X,Y,q,ks));
	INTL += p1*r;
    r += h;

    X = r*cos(T);
    Y = r*sin(T);
    INTF += p2*r*fabs(DET(X,Y,q,ks));
	INTL += p2*r;
    r += h;
    
    for(int i=3;i<=N-3;i++)
    {
 	    X = r*cos(T);
 	    Y = r*sin(T);
        INTF += r*fabs(DET(X,Y,q,ks));
		INTL += r;
        r += h;
    }

    X = r*cos(T);
    Y = r*sin(T);
    INTF += p2*r*fabs(DET(X,Y,q,ks));
	INTL += p2*r;
    r += h;

    X = r*cos(T);
    Y = r*sin(T);
    INTF += p1*r*fabs(DET(X,Y,q,ks));
	INTL += p1*r;
    r += h;

    X = r*cos(T);
    Y = r*sin(T);
    INTF += p0*r*fabs(DET(X,Y,q,ks));
	INTL += p0*r;
    
	*Intl = h*INTL;

    return h*INTF;   
}
#endif
