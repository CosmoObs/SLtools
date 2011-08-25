/** @file
* Example of doxygen documentation for C functions FIXME. 
*/

/** @package theta_find
*  Short description FIXME 
*   
*
*  Detailed description FIXME
*
*/

#ifndef THETA_FIND_H
#define THETA_FIND_H

#include <stdio.h>
#include <stdlib.h>

#include "perturbative_method.h"

void theta_find(f_type Df0Dtheta_in, elliptical_source source_in,  double pert_params[], double _r_e, int N);

void theta_find_2(f_type Df0Dtheta_in, elliptical_source source_in,  double pert_params[], double _r_e, double out[10][2], int N);
#endif
