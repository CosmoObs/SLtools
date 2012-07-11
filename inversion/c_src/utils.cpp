// ==================================
// Author: Gabriel Bartosch Caminha - gbcaminha@gmail.com
// ==================================

//#include "utils.h"

#include <vector>
#include <iostream>
#include <string>


#include "fitsio.h"

template <class T>
bool check_matrix(std::vector<std::vector<T>> data_in)
{
    std::vector<unsigned long int> lin_dim;
    unsigned long int col_dim;

    col_dim = data_in.size();
    for (unsigned long int i = 0; i < col_dim; i++)
    {
        lin_dim.push_back(data_in[i].size());
        if ( lin_dim[0] != lin_dim[i] )
        {
            std::cout << "error: check_matrix: the matrix has multi"
                      << " dimensional lines\n";
            return false;
        }
    }
    std::cout << "check_matrix: The matrix has " << col_dim << " colunms and " 
              << lin_dim[0] << " lines\n";
    return true;
}


bool write_fits_image(std::vector<std::vector<double>> data_in,
                      char fits_name[])
{
    if ( !check_matrix(data_in) )
    {
        std::cout << "error: write_fits_image: data has is not well suited. "
                  << "I will stop here!\n";
        return false;
    }

    int status;
    fitsfile *fptr;
    fits_create_file(&fptr, fits_name, &status);

    unsigned long int dimension = 2, fpixel = 1;
    unsigned long int nelements = data_in.size() *  data_in[0].size();
    long int shape[2] = {data_in.size(), data_in[0].size()};

    double data_c[shape[1]][shape[0]];

    fits_create_img(fptr, DOUBLE_IMG, dimension, shape, &status);

    for (unsigned int i = 0; i < data_in[0].size(); i++)
    {
        for (unsigned long int j = 0; j < data_in.size(); j++)
        {
            //std::cout << std::scientific << data_in[j][i] << "  ";
            data_c[i][j] = data_in[j][i];
        }
        //std::cout << "\n";
    }

    fits_write_img(fptr, TDOUBLE, fpixel, nelements, data_c[0], &status);
    fits_close_file(fptr, &status);
    fits_report_error(stderr, status);

    return true;
}



