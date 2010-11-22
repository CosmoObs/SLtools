/** \file file.h
 * A brief file description.
 * A more elaborated file description.
 */


#ifndef LENTE_H
#define LENTE_H



#define PI 3.1415926535897932384626433832795


/*******************************************************************************
FUNCOES RELACIONADAS AO MODELO DA LENTE
*******************************************************************************/
//! A function used to compute the spherical NFW convergence
/*!
  \param X radial coordinate in units of r_s (r_s is the NFW profile parameter)
  \return \f$ \frac{{\rm ArcTanH}(\sqrt{x^2-1})}{\sqrt{x^2-1}} \f$
  \sa k()
*/
double F(double X)
{
    double f;
    double x2;
    double arg;
    
    x2 = X*X;
    
    if (X < 1.0 - 1e-8)
    {
        arg = sqrt(1.0-x2);
        
        f = 1.0/arg * atanh(arg);
    }
    else if (X > 1.0 + 1e-8)
    {
        arg = sqrt(x2 - 1.0);
        
        f = 1.0/arg * atan(arg);
    }
    else
    {
        f = 1.0;
    }
    
    return f;
}

/******************************************************************************/
//! Spherical NFW convergence
/*!
  \param X radial coordinate in units of r_s (r_s is the NFW profile parameter)
  \param ks NFW profile parameter
  \return the convergence, given by \f$ 2\kappa_s\frac{1-F(X)}{x^2-1} \f$
  \sa F()
*/
double k(double X, double ks)
{
    double K;
    
    if(X == 1.0)
    {
        return 2.0/3.0*ks;
    }
    else
    {
        K = 2.0*ks*(1.0 - F(X))/(X*X - 1.0);
        return K;
    }    
}
/******************************************************************************/
//! Derivative of F
/*!
  \param X radial coordinate in units of r_s (r_s is the NFW profile parameter)
  \return \f$ \frac{dF(x)}{dx} = \frac{1-x^2 F(x)}{x(x^2-1)}\f$
  \sa F()
*/
double Flinha(double X)
{
    double x2;
    double arg;
    double flinha;
    
    x2 = X*X;
    arg = x2 - 1.0;
    
    flinha = (1.0 - x2*F(X))/(X*arg);
    
    
    return flinha;
}
/******************************************************************************/
//! Derivative of K
/*!
  \param X radial coordinate in units of r_s (r_s is the NFW profile parameter)
  \param ks NFW profile parameter
  \return \f$ \frac{dk(x)}{d\left(x^2\right)}\f$
  \sa F()
*/
double klinha(double X, double ks)
{
    double X2=X*X;
	double UMX2=X2-1.0;
	
	if( fabs(X-1.0) < 1e-8) return -0.4*ks;
    return ks*( 2.0*(F(fabs(X))-1.0) - UMX2*Flinha(fabs(X))/X )/(UMX2*UMX2);
}


#endif
