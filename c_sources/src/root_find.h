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

#include <math.h> 

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
double root_find(double f_in(double,void*), void *pot_params, double *est_err_out, 
                 double x_lo, double x_hi, int max_iter, double relative_error, int v);

#endif
