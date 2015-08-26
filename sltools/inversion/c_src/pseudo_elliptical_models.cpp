// ==================================
// Author: Gabriel Bartosch Caminha - gbcaminha@gmail.com
// ==================================
#include <iostream>
#include <cmath>

#include "pseudo_elliptical_models.h"

template <class RadialModel>
class PseudoEllipticalLens
{
    public:
    PseudoEllipticalLens(RadialModel radial_lens_in, double ell_in,
                         std::string parametrization)
        : radial_lens (radial_lens_in),
          _ell (ell_in),
          _parametrization (parametrization),
          _a1 (this->_set_a1()),
          _a2 (this->_set_a2()),
          _sqrt_a2oa1 (sqrt(this->_a2/this->_a1)),
          _a_cap ( 0.5*(this->_a1 + this->_a2) ),
          _b_cap ( 0.5*(this->_a1 - this->_a2) )

    {
        std::cout << "PseudoEllipticalLens:Initializing PseudoEllipricalLens"
                  << "\n";
        std::cout << "    PseudoEllipticalLens: considering radial model as: "
                  << this->radial_lens.lens_name << "\n";
        std::cout << "    PseudoEllipticalLens: input ellipticity: "
                  << this->_ell << "\n";
        std::cout << "    PseudoEllipticalLens: using parametrization: "
                  << this->_parametrization << "\n";
    }
    double kappa(double x_in, double y_in)
    {
        double phie = atan2(y_in*this->_sqrt_a2oa1, x_in);
        double xe = sqrt(this->_a1*x_in*x_in + this->_a2*y_in*y_in);

        double out = this->_a_cap*this->radial_lens.kappa(xe) -
                     this->_b_cap*this->radial_lens.gamma(xe)*cos(2.0*phie);

        if (isnan(out))
        {
            std::cout << "Warning: PseudoEllipticalLens: kappa_ell(" << x_in
                      << ", " << y_in << ") is not a number." << "\n";
        }
        return out;
    }
    double gamma(double x_in, double y_in)
    {
        double phie = atan2(y_in*this->_sqrt_a2oa1, x_in);
        double xe = sqrt(this->_a1*x_in*x_in + this->_a2*y_in*y_in);

        double kappa_rad = this->radial_lens.kappa(xe);
        double gamma_rad = this->radial_lens.gamma(xe);

        double cons_1 = pow(this->_a_cap*gamma_rad, 2);
        double cons_2 = 2.0*this->_a_cap*this->_b_cap*kappa_rad*gamma_rad*
                        cos(2.0*phie);
        double cons_3 = pow(this->_b_cap, 2);
        double cons_4 = pow(kappa_rad, 2) - pow(sin(2.0*phie)*gamma_rad, 2);


        return cons_1 - cons_2 + cons_3*cons_4;
    }
    double gamma_1(double x_in, double y_in)
    {
        double phie = atan2(y_in*this->_sqrt_a2oa1, x_in);
        double xe = sqrt(this->_a1*x_in*x_in + this->_a2*y_in*y_in);
        return this->_b_cap*this->radial_lens.kappa(xe)
                     - this->_a_cap*this->radial_lens.gamma(xe)*cos(2.0*phie);
    }
    double gamma_2(double x_in, double y_in)
    {
        double phie = atan2(y_in*this->_sqrt_a2oa1, x_in);
        double xe = sqrt(this->_a1*x_in*x_in + this->_a2*y_in*y_in);
        return -sqrt( pow(this->_a_cap, 2) + pow(this->_b_cap, 2)  )
               *this->radial_lens.gamma(xe)*sin(2.0*phie);
    }


    RadialModel radial_lens;

    private:

    double _set_a1()
    {
        if ( this->_parametrization.compare("G&K") == 0 )
        {
            return 1.0 - this->_ell;
        }
        else
        {
            std::cout << "error: PseudoEllipticalLens: _set_a1: "
                      << "parametrization not defined, assuming G&K\n";
            return 1.0 - this->_ell;
        }
    }
    double _set_a2()
    {
        if ( this->_parametrization.compare("G&K") == 0 )
        {
            return 1.0 + this->_ell;
        }
        else
        {
            std::cout << "error: PseudoEllipticalLens: _set_a2: "
                      << "parametrization not defined, assuming G&K\n";
            return 1.0 + this->_ell;
        }

    }
    double _ell;
    std::string _parametrization;
    double _a1, _a2;
    double _sqrt_a2oa1;
    double _a_cap, _b_cap;
};

/*class PellSis
{
    public:
    PellSis(SisLens sis_in)
        : psis (sis_in)
    {

    }
    PseudoEllipticalLens<SisLens> psis;
};*/
