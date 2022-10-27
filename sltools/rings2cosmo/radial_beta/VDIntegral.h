#ifndef _VDINTEGRAL_H_

#define _VDINTEGRAL_H_

int Integrand(const int *ndim, const double xx[],
              const int *ncomp, double ff[], void *userdata);
#
double Integral_func (double b, double alpha, double beta, double delta);

#endif