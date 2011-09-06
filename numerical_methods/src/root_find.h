/** @file
* Example of doxygen documentation for C functions FIXME. 
*/

/** @package general_methods
*  Package with general functions, root finders, integrators etc.
*
*  Detailed descrition FIXME
*
*/

#ifndef ROOT_FIND_H
#define ROOT_FIND_H

#include <cmath> 

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

//!  Root finder that use the gsl library, (brent method)
//!
//!  For a detailed description see: http://www.gnu.org/software/gsl/manual/html_node/One-dimensional-Root_002dFinding.html. To compile: g++ -Wall  main.cpp `pkg-config gsl --cflags --libs`
/*!
  \param f_in function in which the root will be found
  \param pot_params function parameters
  \param est_err vector where the estimated error will be written
  \param x_lo lower value to be considered by the root finder 
  \param x_hi upper value to be considered by the root finder
  \param max_iter maximum number of interactions
  \param relative_error relative error
  \param v !=0 for verbose mode
  \return the root
  \return est_err_out
*/
double root_find(double f_in(double,void*), void *pot_params, double *est_err_out=NULL, 
                 double x_lo = 1E-4, double x_hi = 10.0, int max_iter = 100, double relative_error = 1E-4, int v=0){
  double r = 0.0;

  gsl_function F;

  F.function = f_in;
  F.params = pot_params;

  int status;
  int iter = 0;

  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc (T);
  gsl_root_fsolver_set (s, &F, x_lo, x_hi);

  if(v) printf ("using %s method\n",gsl_root_fsolver_name (s));
  if(v) printf ("%5s [%9s, %9s] %9s %9s\n","iter", "lower", "upper", "root", "err(est)");

  do{
    iter++;
    status = gsl_root_fsolver_iterate (s);
    r = gsl_root_fsolver_root (s);
    x_lo = gsl_root_fsolver_x_lower (s);
    x_hi = gsl_root_fsolver_x_upper (s);
    status = gsl_root_test_interval (x_lo, x_hi,0, relative_error);

    if (v && status == GSL_SUCCESS) printf("Converged:\n");
        if(v) printf ("%5d [%.7f, %.7f] %.7f %.7f\n",iter, x_lo, x_hi, r, x_hi - x_lo);
  }while (status == GSL_CONTINUE && iter < max_iter);

  if(est_err_out != NULL) *est_err_out = x_hi - x_lo;

  return r;
}

#endif
