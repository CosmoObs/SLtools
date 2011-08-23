/** @file
* Example of doxygen documentation for C functions FIXME. 
*/

/** @package propriedades
*  Package.
*
*  Detailed descrition FIXME
*
*/

#ifndef PROPRIEDADES_H
#define PROPRIEDADES_H


#define PI 3.1415926535897932384626433832795

/*---CALCULANDO AS DERIVADAS DO POTENCIAL DA LENTE*****************************/
double fiX(double X, double Y, double q, double ks){return X*J(0,X,Y,q,ks);}
double fiY(double X, double Y, double q, double ks){return Y*J(1,X,Y,q,ks);}

double fiXX(double X, double Y, double q, double ks){return 2.0*X*X*K(0,X,Y,q,ks) + J(0,X,Y,q,ks);}
double fiYY(double X, double Y, double q, double ks){return 2.0*Y*Y*K(2,X,Y,q,ks) + J(1,X,Y,q,ks);}
double fiXY(double X, double Y, double q, double ks){return 2.0*X*Y*K(1,X,Y,q,ks);}
/******************************************************************************/

/*---EQUACAO DA LENTE---*******************************************************/
double XF(double X, double Y, double q, double ks){return X*(1.0 - J(0,X,Y,q,ks));}
double YF(double X, double Y, double q, double ks){return Y*(1.0 - J(1,X,Y,q,ks));}
/******************************************************************************/

/*---CALCULANDO OS AUTOVALORES E ALGUMAS RELACOES---***************************/
double L1(double X, double Y, double q, double ks)
{
    double Phixx, Phiyy, Phixy;
    double PhixxMPhiyy;
    
    Phixx = fiXX(X,Y,q,ks);
    Phiyy = fiYY(X,Y,q,ks);
    Phixy = fiXY(X,Y,q,ks);
    PhixxMPhiyy = Phixx-Phiyy;
    
    return (2.0 - Phixx - Phiyy - sqrt(4.0*Phixy*Phixy + PhixxMPhiyy*PhixxMPhiyy) )/2.0;
}
/******************************************************************************/
double L2(double X, double Y, double q, double ks)
{
    double Phixx, Phiyy, Phixy;
    double PhixxMPhiyy;
    
    Phixx = fiXX(X,Y,q,ks);
    Phiyy = fiYY(X,Y,q,ks);
    Phixy = fiXY(X,Y,q,ks);
    PhixxMPhiyy = Phixx-Phiyy;
    
    return (2.0 - Phixx - Phiyy + sqrt(4.0*Phixy*Phixy + PhixxMPhiyy*PhixxMPhiyy) )/2.0;
}
/******************************************************************************/
double L2OL1(double X, double Y, double q, double ks)
{
    double Phixx, Phiyy, Phixy;
    double PhixxMPhiyy;
    
    Phixx = fiXX(X,Y,q,ks);
    Phiyy = fiYY(X,Y,q,ks);
    Phixy = fiXY(X,Y,q,ks);
    PhixxMPhiyy = Phixx-Phiyy;
    
    return (2.0 - Phixx - Phiyy + sqrt(4.0*Phixy*Phixy + PhixxMPhiyy*PhixxMPhiyy) )/
           (2.0 - Phixx - Phiyy - sqrt(4.0*Phixy*Phixy + PhixxMPhiyy*PhixxMPhiyy) );
}
/******************************************************************************/
double L1OL2(double X, double Y, double q, double ks)
{
    double Phixx, Phiyy, Phixy;
    double PhixxMPhiyy;
    
    Phixx = fiXX(X,Y,q,ks);
    Phiyy = fiYY(X,Y,q,ks);
    Phixy = fiXY(X,Y,q,ks);
    PhixxMPhiyy = Phixx-Phiyy;
    
    return (2.0 - Phixx - Phiyy - sqrt(4.0*Phixy*Phixy + PhixxMPhiyy*PhixxMPhiyy) )/
           (2.0 - Phixx - Phiyy + sqrt(4.0*Phixy*Phixy + PhixxMPhiyy*PhixxMPhiyy) );
}
/******************************************************************************/
double DET(double X, double Y, double q, double ks)
{
    double Phixx, Phiyy, Phixy;
    
    Phixx = fiXX(X,Y,q,ks);
    Phiyy = fiYY(X,Y,q,ks);
    Phixy = fiXY(X,Y,q,ks);

    return (1.0-Phixx)*(1.0-Phiyy) - Phixy*Phixy;
}
/******************************************************************************/

#endif
