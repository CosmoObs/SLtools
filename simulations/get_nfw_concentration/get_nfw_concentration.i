%module get_nfw_concentration
%{
#include"./rand.h"
extern void read_config(char* get_nfw_concentration_conf);
extern float c(float m, float z, long int start);
extern float mean_log10c(float m, float z);
extern float c(float mass, float redshift, long int start);
extern float sigma_log10c(float m, float z);
extern void help(void);

%}
extern void read_config(char* get_nfw_concentration_conf);
extern float mean_log10c(float m, float z);
extern float c(float mass, float redshift, long int start);
extern float sigma_log10c(float m, float z);
extern void help(void);
