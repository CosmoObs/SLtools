// ==================================
// Author: Gabriel Bartosch Caminha - gbcaminha@gmail.com
// ==================================

#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>

#include "composit.h"

#include "radial_models.cpp"
#include "pseudo_elliptical_models.cpp"

class CompositeLens
{
    public:
    CompositeLens(std::vector<SisLens> sis_vec_in,
                  std::vector<TestClass> test_vec_in)
        : sis_vec  ( sis_vec_in ),
          test_vec ( test_vec_in ),
          psis_vec (),
          ptest_vec (),
          n_sis  ( this->sis_vec.size() ),
          n_test ( this->test_vec.size() )
    {
        std::cout << "CompositeLens: Initializing\n";
        std::cout << "    CompositeLens: number of sis lenses: " << this->n_sis
                  << "\n";
        std::cout << "    CompositeLens: number of test lenses: "
                  << this->n_test << "\n";
    }
    void set_ell_sis(std::vector<double> ell_in)
    {
        for (unsigned int i = 0; i < this->n_sis; i++)
            this->psis_vec.push_back(PseudoEllipticalLens<SisLens>
                                     (this->sis_vec[i], ell_in[i], "G&K" ));
    }
    double kappa(double x_in, double y_in)
    {
        double kappa_out = 0.0;
        double radial = sqrt(pow(x_in, 2) + pow(y_in, 2));
        for (auto i : this->sis_vec)  kappa_out += i.kappa(radial);
        for (auto i : this->test_vec) kappa_out += i.kappa(radial);
        return kappa_out;
    }

    std::vector<SisLens> sis_vec;
    std::vector<TestClass> test_vec;

    std::vector< PseudoEllipticalLens<SisLens> > psis_vec;
    std::vector< PseudoEllipticalLens<TestClass> > ptest_vec;
    unsigned int n_sis, n_test;
};

class CompositePsis
{
    public:
    CompositePsis( std::vector<double> sis_par_in, std::vector<double> ell_in,
                   std::vector<double> theta_in,   std::vector<double> x0_in,
                   std::vector<double> y0_in )
        : sis_par (sis_par_in),
          ell (ell_in),
          theta (theta_in),
          x0 (x0_in),
          y0 (y0_in),
          psis_vec (),
          n_lens (this->sis_par.size())
    {
        init_psis();
    }
    void init_psis()
    {
        for (unsigned int i = 0; i < this->n_lens; i++)
        {
            this->psis_vec.push_back(PseudoEllipticalLens<SisLens>
                                     (SisLens(this->sis_par[i]), this->ell[i],
                                      "G&K") );
        }
    }
    void lens_info()
    {
        std::cout << "CompositePsis: Number of lenses: " << this->n_lens
                  << "\n";
        std::cout << "Lens      sis_par          ell        theta           x0"
                  << "           y0\n";
        for (unsigned int i = 0; i < this->n_lens; i++)
        {
            std::cout << "   " << i << std::scientific << " "
                      << this->sis_par[i] << " " << this->ell[i] << " "
                      << this->theta[i] << " " << this->x0[i] << " "
                      << this->y0[i] << "\n";
        }
    }

    double kappa(double x_in, double y_in)
    {
        double kappa_out = 0.0;
        double x_eff = 0.0;
        double y_eff = 0.0;
        double theta_ef = 0.0;

        for (unsigned int i = 0; i < this->n_lens; i++)
        {
            theta_ef = this->theta[i]/180.0 * M_PIl;
            x_eff = (x_in - this->x0[i])*cos(theta_ef) -
                    (y_in - this->y0[i])*sin(theta_ef);
            y_eff = (x_in - this->x0[i])*sin(theta_ef) +
                    (y_in - this->y0[i])*cos(theta_ef);
            kappa_out += this->psis_vec[i].kappa(x_eff, y_eff);
        }
        return kappa_out;
    }
    double gamma(double x_in, double y_in)
    {
        double gamma_1 = 0.0;
        double gamma_2 = 0.0;
        double x_eff = 0.0;
        double y_eff = 0.0;
        double theta_ef = 0.0;
        for (unsigned int i = 0; i < this->n_lens; i++)
        {
            theta_ef = this->theta[i]/180.0 * M_PIl;
            x_eff = (x_in - this->x0[i])*cos(theta_ef) -
                    (y_in - this->y0[i])*sin(theta_ef);
            y_eff = (x_in - this->x0[i])*sin(theta_ef) +
                    (y_in - this->y0[i])*cos(theta_ef);
            gamma_1 += this->psis_vec[i].gamma_1(x_eff, y_eff);
            gamma_2 += this->psis_vec[i].gamma_2(x_eff, y_eff);
        }
        return sqrt( pow(gamma_1, 2) + pow(gamma_2, 2) );
    }

    std::vector<double> sis_par, ell, theta, x0, y0;
    std::vector< PseudoEllipticalLens<SisLens> > psis_vec;
    unsigned int n_lens;
};
