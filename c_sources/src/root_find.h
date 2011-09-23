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




double root_find(double f_in(double,void*), void *pot_params, double *est_err_out, 
                 double x_lo, double x_hi, int max_iter, double relative_error, int v);

#endif
