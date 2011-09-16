/** @file
* Example of doxygen documentation for C functions FIXME. 
*/

/** @package perturbative_method
*  Package to generate the arcs
*
*  Detailed descrition FIXME
*
*/


#ifndef GENERATE_ARCS_H
#define GENERATE_ARCS_H

#include <stdio.h>

#include "perturbative_method.h"
#include "theta_find.h"


void plot_arcs(elliptical_source source_in, f_type f1_in, f_type Df0Dtheta_in, double pert_params[], double kappa2, double _r_e, int npts, FILE *file_in);

void plot_arcs_sep(elliptical_source source_in, f_type f1_in, f_type Df0Dtheta_in, double pert_params[], double kappa2, double _r_e, int npts, FILE *file_in);

#endif
