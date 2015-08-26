/** @file
* Example of doxygen documentation for C functions FIXME. 
*/

/** @package general_methods
*  Package with general functions, root finders, integrators etc.
*
*  Detailed descrition FIXME
*
*/


//!  Function to calculate the determinant of some \f$ 3\times 3 \f$ matrix
/*!
  \param A[][] matrix array
  \return det_mx determinant
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h> 

double det_matrix (double A[3][3]){
  double det_pos,det_neg,det_mtx;
  double m_11=A[0][0],m_12=A[0][1],m_13=A[0][2];
  double m_21=A[1][0],m_22=A[1][1],m_23=A[1][2];
  double m_31=A[2][0],m_32=A[2][1],m_33=A[2][2];

  det_pos=m_11*m_22*m_33 + m_12*m_23*m_31 + m_13*m_21*m_32;

  det_neg=m_13*m_22*m_31+ m_12*m_21*m_33 + m_11*m_23*m_32;

  det_mtx=det_pos-det_neg;
  
  return det_mtx;

}

//!  Function to solve a equation systems \f$ A*X=B \f$, using the Cramer's Rule. 
//!
//! A is a matrix of \f$ 3\times 3 \f$, i.e, \f$ A_{[3\times3]}\f$
//!
//! X is a vector solution, i.e \f$ X_{[3\times1]}\f$ 
//!  
/*!
  \param A[][] bidimensional array
  \param B[]  array containing some things.
  \return X[]  array containing the vector solution
*/

void solving_system(double A[3][3], double B[3], double X[3]){
  double trp_x1[3][3],trp_x2[3][3],trp_x3[3][3];
  double det_a, det_x1, det_x2, det_x3; 
//   double X[3];

  trp_x1[0][0]=B[0],trp_x1[0][1]=A[0][1],trp_x1[0][2]=A[0][2];
  trp_x1[1][0]=B[1],trp_x1[1][1]=A[1][1],trp_x1[1][2]=A[1][2];
  trp_x1[2][0]=B[2],trp_x1[2][1]=A[2][1],trp_x1[2][2]=A[2][2];

  det_x1=det_matrix(trp_x1);

  trp_x2[0][0]=A[0][0],trp_x2[0][1]=B[0],trp_x2[0][2]=A[0][2];
  trp_x2[1][0]=A[1][0],trp_x2[1][1]=B[1],trp_x2[1][2]=A[1][2];
  trp_x2[2][0]=A[2][0],trp_x2[2][1]=B[2],trp_x2[2][2]=A[2][2];

  det_x2=det_matrix(trp_x2);

  trp_x3[0][0]=A[0][0],trp_x3[0][1]=A[0][1],trp_x3[0][2]=B[0];
  trp_x3[1][0]=A[1][0],trp_x3[1][1]=A[1][1],trp_x3[1][2]=B[1];
  trp_x3[2][0]=A[2][0],trp_x3[2][1]=A[2][1],trp_x3[2][2]=B[2];

  det_x3=det_matrix(trp_x3);

  det_a=det_matrix(A);
//   printf("\n The main determinant is %f\n",det_a);

  if(det_a!=0.0)
    X[0]=det_x1/det_a,  
    X[1]=det_x2/det_a,
    X[2]=det_x3/det_a;
//     printf("The solution is x[1] = %f,\tx[2] = %f,\tx[3] = %f\n",X[0],X[1],X[2]);  
  else
   printf("\nThe system have infinite solutions\n"),
   exit(1);
  
// return 0.0;
}
