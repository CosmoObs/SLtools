%module get_nfw_concentration
%{
#include"./rand.h"
extern void read_config(char* get_nfw_concentration_conf);
extern float c(float mass, float redshift, long int start);
extern void help(void);

%}
extern void read_config(char* get_nfw_concentration_conf);
extern float c(float mass, float redshift, long int start);
extern void help(void);
