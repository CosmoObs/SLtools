/** @file
* Example of doxygen documentation for C functions FIXME. 
*/

/** @package general_methods
*  Package with general functions, root finders, integrators etc.
*
*  Detailed descrition FIXME
*
*/

#ifndef SOLVING_SYSTEM_H
#define SOLVING_SYSTEM_H

#include <stdlib.h>
#include <stdio.h>
//#include <fstream>        //funcoes de entrada e saida para arquivos
//#include <iostream>
#include <math.h> 

double det_matrix (double A[3][3]);
void solving_system(double A[3][3], double B[3], double X[3]);

#endif
