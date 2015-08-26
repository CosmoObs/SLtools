// ==================================
// Author: Gabriel Bartosch Caminha - gbcaminha@gmail.com
// ==================================

#include <iostream>
#include <string>

#include "radial_models.h"

class SisLens
{
    public:
    SisLens(double sis_par_in)
        : lens_name ("SIS"),
          _sis_par (sis_par_in)

    {
        std::cout << "SisLens:Initializing SisLens\n";
        std::cout << "    SisLens:Sis parameter = " << this->_sis_par << "\n";
    }
    /**************************************************************************/
    double kappa(double radial)
    {
        return this->_sis_par/2.0/radial;
    }
    /**************************************************************************/
    double gamma(double radial)
    {
        return this->_sis_par/2.0/radial;
    }
    /**************************************************************************/
    double alpha()
    {
        return this->_sis_par;
    }
    /**************************************************************************/
    int set_sis_par(double sis_par_in)
    {
        this->_sis_par = sis_par_in;
        std::cout << "SisLens: setting sis_par = " << this->_sis_par << "\n";
        return 0;
    }
    /**************************************************************************/
    double get_sis_par()
    {
        return this->_sis_par;
    }
    /**************************************************************************/
    double test(double a_in)
    {
        std::cout << "SisLens: test: " << a_in << "\n";
        return a_in*a_in;
    }


    std::string lens_name;
    private:
    double _sis_par;
};

class TestClass
{
    public:
    TestClass(double a_in)
        : lens_name ("test_class")
    {
        std::cout << "TestClass: Initializing: " << a_in << "\n";
    }
    double test(double a_in)
    {
        std::cout << "TestClass: test: " << a_in << "\n";
        return a_in*a_in;
    }
    double kappa(double radial)
    {
        return 0.5/radial;
    }
    double gamma(double radial)
    {
        return 1.0/radial;
    }

    std::string lens_name;
};

