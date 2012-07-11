#include <vector>
#include <fstream>
#include <iostream>

#include "create_grid.h"
#include "utils.h"

int create_grid (double x_min, double x_max, int x_pts,
                 double y_min, double y_max, int y_pts,
                 std::vector <double> &vec_x, std::vector <double> &vec_y,
                 bool begin_edge)
{
    double x_step = (x_max - x_min)/double( begin_edge ? x_pts - 1 : x_pts );
    double y_step = (y_max - y_min)/double( begin_edge ? y_pts - 1 : y_pts );

    double x_var = ( begin_edge ? x_min : x_min + x_step/2.0 );
    double y_var = ( begin_edge ? y_min : y_min + y_step/2.0 );

    if (! vec_x.empty ())
    {
        std::cout << "Vector X is not empty\n";
        vec_x.erase (vec_x.begin (), vec_x.end ());
    }
    if (! vec_y.empty ())
    {
        std::cout << "Vector Y is not empty\n";
        vec_y.erase (vec_y.begin (), vec_y.end ());
    }

    for (int i_x = 0; i_x < x_pts; i_x++)
    {
        for (int i_y = 0; i_y < y_pts; i_y++)
        {
            vec_x.push_back (x_var);
            vec_y.push_back (y_var);
            y_var += y_step;
        }
        y_var = ( begin_edge ? y_min : y_min + y_step/2.0 );
        x_var += x_step;
    }
    return 0;
}
template <class LensModelTemplate> //should have the function kappa(x_in, y_in)
void kappa_grid( LensModelTemplate lens_in )
{
    std::vector <double> vec_x;
    std::vector <double> vec_y;
    std::vector <double> vec_func;

    std::vector<std::vector<double>> grid_fits;
    std::vector<double> vec_fits;

    double x_min = -1.0;
    double x_max = 1.0;

    long double y_min = -1.0;
    long double y_max = 1.0;

    int x_step = 201;
    int y_step = 201;

    create_grid (x_min, x_max, x_step, y_min, y_max, y_step, vec_x, vec_y,
                 false);

    std::ofstream file_out;
    file_out.open("kappa_grid.dat");

    file_out << std::scientific;
    file_out << "# x, y, kappa(x, y)\n";
    for (unsigned int i = 0; i < vec_x.size(); i++)
    {
        vec_func.push_back( lens_in.kappa(vec_x[i], vec_y[i]) );
        file_out << vec_x[i] << "\t";
        file_out << vec_y[i] << "\t";
        file_out << vec_func[i] << "\n";
    }

    file_out.close();
}


/*class MakeGrid
{
    MakeGrid()
    {
        std::cout << "Initializing MakeGrid\n";
    }
};*/
