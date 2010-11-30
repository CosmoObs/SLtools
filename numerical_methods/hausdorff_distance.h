/** @file
* Example of doxygen documentation for C functions FIXME. 
*/

/** @package hausdorff_distance.h
*  Package to compute distances between points and curves. See http://twiki.linea.gov.br/bin/view/StrongLensing/HausdorffDistance
*
*  Detailed descrition FIXME
*
*/

#ifndef HAUSDORFF_DISTANCE_H
#define HAUSDORFF_DISTANCE_H

#include <cstdio>
#include <cmath>
#include <algorithm>
#include <cstdlib>


// to compile make: g++ -Wall -o exec_h hausdorff.cpp `pkg-config gsl --cflags --libs`


/** \brief Compute the distance between two points
*
*  \param x1 x coordinate for first point
*  \param y1 y coordinate for first point
*  \param x2 x coordinate for second point
*  \param y2 y coordinate for second point
*  \return \f$ \sqrt{ (x_1-x_2)^2 +  (y_1-y_2)^2} \f$
*/
double dist2pts(double x1, double y1, double x2, double y2){return sqrt( pow(x1-x2,2) + pow(y1-y2,2)  );}


/** \brief Minimum distance between a points and a curve
*
*  \param ptx x coordinate for point
*  \param pty y coordinate for point
*  \param cx[] x coordinates values that define the curve
*  \param cy[] y coordinates values that define the curve
*  \param cnpts number of points
*  \return Minimum distance between (ptx,pty) and the curve defined by {(cx[],cy[])} 
*/
double dist_pt_curve(double ptx, double pty, double cx[], double cy[], int cnpts){
  double dist[cnpts];
  for(int i=0;i<cnpts;i++) dist[i] = dist2pts(ptx, pty, cx[i], cy[i]);
  return *std::min_element(dist,dist+cnpts);
}


/** \brief Maximum distance between a points and a curve
*
*  \param ptx x coordinate for point
*  \param pty y coordinate for point
*  \param cx[] x coordinates values that define the curve
*  \param cy[] y coordinates values that define the curve
*  \param cnpts number of points
*  \return Maximum distance between (ptx,pty) and the curve defined by {(cx[],cy[])} 
*/
double dist_max_pt_curve(double ptx, double pty, double cx[], double cy[], int cnpts){
  double dist[cnpts];
  for(int i=0;i<cnpts;i++) dist[i] = dist2pts(ptx, pty, cx[i], cy[i]);
  return *std::max_element(dist,dist+cnpts);
}


/** \brief Function to define the distance between two curves.
*
*  \param cx1[] x coordinates values that define the curve 1
*  \param cy1[] y coordinates values that define the curve 1
*  \param c1npts number of points of curve 1
*  \param cx2[] x coordinates values that define the curve 2
*  \param cy2[] y coordinates values that define the curve 2
*  \param c2npts number of points of curve 2
*  \param type specify the type of distance, see http://twiki.linea.gov.br/bin/view/StrongLensing/HausdorffDistance
*  \return distance between the two curves
*/
double dist2curves(double c1x[], double c1y[], int c1npts, double c2x[], double c2y[], int c2npts, int type = 3){

  double dist_vec[c1npts];
  double dist_out=0.0; 
  for(int i=0;i<c1npts;i++) dist_vec[i] = dist_pt_curve(c1x[i], c1y[i], c2x, c2y, c2npts);

  if(type == 1){
    dist_out = *std::min_element(dist_vec,dist_vec+c1npts);;
  } else if(type == 2) {
    dist_out = *std::max_element(dist_vec,dist_vec+c1npts);;
  } else if(type == 3){
    for(int i=0;i<c1npts;i++) dist_out += dist_vec[i]/c1npts;
  } else {
    printf("Using definition d3 in http://twiki.linea.gov.br/bin/view/StrongLensing/HausdorffDistance\n");
    for(int i=0;i<c1npts;i++) dist_out += dist_vec[i]/c1npts;
  }

  return dist_out;
}


/** \brief Compute the maximum distance between two curves
*
*  \param cx1[] x coordinates values that define the curve 1
*  \param cy1[] y coordinates values that define the curve 1
*  \param c1npts number of points of curve 1
*  \param cx2[] x coordinates values that define the curve 2
*  \param cy2[] y coordinates values that define the curve 2
*  \param c2npts number of points of curve 2
*  \return Maximum distance between the two curves
*/
double dist2curves_max(double c1x[], double c1y[], int c1npts, double c2x[], double c2y[], int c2npts){

  double dist_vec[c1npts];
  double dist_out=0.0; 
  for(int i=0;i<c1npts;i++) dist_vec[i] = dist_max_pt_curve(c1x[i], c1y[i], c2x, c2y, c2npts);

  dist_out = *std::max_element(dist_vec,dist_vec+c1npts);;

  return dist_out;
}

/** \brief Compute the f_4 modified Hausdorff distance between two curves. (see http://twiki.linea.gov.br/bin/view/StrongLensing/HausdorffDistance)
*
*  \param cx1[] x coordinates values that define the curve 1
*  \param cy1[] y coordinates values that define the curve 1
*  \param c1npts number of points of curve 1
*  \param cx2[] x coordinates values that define the curve 2
*  \param cy2[] y coordinates values that define the curve 2
*  \param c2npts number of points of curve 2
*  \param type specify the type of distance, see http://twiki.linea.gov.br/bin/view/StrongLensing/HausdorffDistance
*  \return f_4 modified Hausdorff distance
*/
double haus_dist_f4(double c1x[], double c1y[], int c1npts, double c2x[], double c2y[], int c2npts, int type = 3){
  
  double distAB = dist2curves(c1x, c1y, c1npts, c2x, c2y, c2npts, type);
  double distBA = dist2curves(c2x, c2y, c2npts, c1x, c1y, c1npts, type);

  return (double(c1npts)*distAB + double(c2npts)*distBA)/double(c1npts+c2npts);

}


#endif
