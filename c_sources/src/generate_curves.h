/** @file
* Module to generate curves useful for the perturbative approach, as
* 
*  Tangential Critical Curve and  Tangential Caustic. This function need
* - main field of the perturbative approach, i.e, \f$ f_n(\theta)\f$ and their derivatives
* - pert_params[] :a vector that contains all the perturbation parameters
* - npts: number of theta divisions 
* - Einstein Radius: \f$ R_{_{\mathrm{E}}}\f$ 
* - two files to write the tangential curves 
*
*
*  Elliptical source, contour plot. It need
*  - structure corresponding to the source: size, ellipticity, orientation,etc
*  - npts: number of theta divisions
*  - a file to write the plot in a text file
*
*
*  Curves of Constant Distortion. This function need
* - main field of the perturbative approach, i.e, \f$ f_n(\theta)\f$ and their derivatives
* - pert_params[] :a vector that contains all the perturbation parameters
* - npts: number of theta divisions 
* - Einstein Radius: \f$ R_{_{\mathrm{E}}}\f$ 
* - A POSITIVE VALUE of the length-to-width ratio threshold, i.e \f$ R_{\rm th} \f$
* - two files to write theses curves in a text file (with the following format x1, x2, y1, y2), each for both value of \f$ R_{\rm th}\f$
*
* This function computes automatically both positive and negative constant distortion curves
*
*/

/** @package perturbative_method
* 
* To plot the critical and caustic curves (model independent)
* 
* To plot of an elliptical source, not aligned to the main axis, in its corresponding plane 
*/
#ifndef GENERATE_CUVES_H
#define GENERATE_CURVES_H

#include <stdio.h> 
#include <stdlib.h> 
#include "perturbative_method.h"

void plot_curves(f_type f1_in, f_type Df0Dtheta_in, f_type D2f0Dtheta2_in, double pert_params[], double kappa2, double _r_e, int npts, FILE *file_in1);

void plot_sources(elliptical_source source_in, int npts, FILE *file_in);

void plot_rth_curves(f_type f1_in, f_type Df0Dtheta_in, f_type D2f0Dtheta2_in, double pert_params[], double kappa2, double _r_e, double r_th, int npts, FILE *file_in1, FILE *file_in2);

#endif
