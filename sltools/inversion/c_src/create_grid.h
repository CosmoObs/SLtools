#ifndef CREATE_GRID_H
#define CREATE_GRID_H

int create_grid (double x_min, double x_max, int x_pts,
                 double y_min, double y_max, int y_pts,
                 std::vector <double> &vec_x, std::vector <double> &vec_y,
                 bool begin_edge = true );
#include "create_grid.cpp"
//template <class LensModelTemplate> //should have the function kappa(x_in, y_in)
//void kappa_grid( LensModelTemplate lens_in )


#endif /*create_grid.h*/
