//gcc -g -shared -Wall -fPIC -I/home/grasi/Tools/Cuba libVDIntegral.c -o libVDIntegral.so -lm -L/home/grasi/Tools/Cuba -lcuba -lgsl -lgslcblas
#include "VDIntegral.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cuba.h"
#include <float.h>
#include <gsl/gsl_sf_gamma.h> // flags to GSL -lgsl -lgslcblas



int Integrand(const int *ndim, const double xx[],
              const int *ncomp, double ff[], void *userdata) {


  double b = ((double*)userdata)[0];
  double alpha = ((double*)userdata)[1];
  //double beta = ((double*)userdata)[2];
  double delta = ((double*)userdata)[2];
  double ra = ((double*)userdata)[3];

  double xi = alpha + delta - 2.0;

  double v = xx[0];
  double w = xx[1];

  double R = w/(1 - w);
  double r = R/v;
  double x  = R/b;
  double weight_term = exp(-(pow(x,2.0))/2.0);
  double projection_term = sqrt(1.0 - pow(v,2.0));
  
  double Jeans_term = pow(r, 2.0- xi);

  double hyper_cube_term = 1/pow(1.0 - w, 2.0);


  //double w_term = pow(R, 2.0 - xi)/pow(1.0 - w, 2.0);
  //double v_term = pow(v, xi - 2.0);
  
  //----------anisotropy----------------
  // ---------beta cte -----------------
  //double anisotropy = 1.0 - (beta*pow(v,2.0));
  //
  //---------- mamon's beta ---------
  double y = ra/(r + ra);
  double z = pow(R,2.0)/r;
  double anisotropy = xi - (xi*z) - (y) + (z*y);



  ff[0] = 2.0*weight_term*Jeans_term*hyper_cube_term*anisotropy/projection_term;

  return 0;
}
#define STATEFILE NULL
#define SPIN NULL

double Integral_func (double b, double alpha, double delta, double ra){

    double params [] = {b, alpha, delta, ra};

    int NDIM = 2;
    int NCOMP = 1;
    int NVEC = 1;
    double EPSREL = 1e-5;
    double EPSABS = 1e-20;
    int VERBOSE  = 0;
    int MINEVAL = 0;
    int MAXEVAL = 50000;
    int KEY =  13;

    int nregions, neval, fail;
    double integral_result[NCOMP], error[NCOMP], prob[NCOMP];

    Cuhre(NDIM, NCOMP, Integrand, params, NVEC, EPSREL, EPSABS, VERBOSE,
        MINEVAL, MAXEVAL, KEY, STATEFILE, SPIN,
        &nregions, &neval, &fail, integral_result, error, prob);

    return integral_result[0];

}


