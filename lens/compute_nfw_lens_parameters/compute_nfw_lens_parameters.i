%module compute_nfw_lens_parameters 
%{
extern double kappa_s(double zl, double zs, double M, double C, double Om0, double Ol0);
extern double x_s(double zl, double M, double C, double Om0, double Ol0);
extern int main(int argc, char **argv);
%}

extern double kappa_s(double zl, double zs, double M, double C, double Om0, double Ol0);
extern double x_s(double zl, double M, double C, double Om0, double Ol0);
extern int main(int argc, char **argv);
