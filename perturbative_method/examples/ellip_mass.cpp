#include <cstdio>
#include <cstdlib>

#include <glib.h>

#include "nfw_circular_model.h" //import functions for nfw model
#include "potential_derivatives.h" //import functions to compute the lens_potential derivatives


static gboolean compute_enfw_potential_derivatives_in = FALSE;
static gdouble r_s_in = 1.0;
static gdouble kappa_s_in = 1.0;

static gdouble X_lens_min_in = -1.0;
static gdouble X_lens_max_in = 1.0;
static gint X_divisions_in = 10;

static gdouble Y_lens_min_in = -1.0;
static gdouble Y_lens_max_in = 1.0;
static gint Y_divisions_in = 10;


typedef struct _Lens_Plane_Grid_Entries
{
  gdouble teste1;
  gdouble teste2;
} Lens_Plane_Grid_Entries;



GOptionGroup *
grid_group (Lens_Plane_Grid_Entries *grid)
{
  GOptionEntry grid_entries[] =
  {
    { "r_s", NULL, 0, G_OPTION_ARG_DOUBLE, &grid->teste1, "scale radius for nfw model", "1.0" },
    { "r_s", NULL, 0, G_OPTION_ARG_DOUBLE, &grid->teste2, "scale radius for nfw model", "1.0" },
    { NULL }
  };
  GOptionGroup *grid_group = g_option_group_new ("test", " - test options", "Show help options related to test", NULL, NULL);
  g_option_group_add_entries (grid_group, grid_entries);

  return grid_group;
}

static GOptionEntry entries[] = 
{
  { "compute_enfw_pot", 'p', 0, G_OPTION_ARG_NONE, &compute_enfw_potential_derivatives_in, "compute the potential derivative for the enfw model", NULL },
  { "kappa_s", NULL, 0, G_OPTION_ARG_DOUBLE, &kappa_s_in, "characteristic convergence for nfw model", "1.0" },
  { "r_s", NULL, 0, G_OPTION_ARG_DOUBLE, &r_s_in, "scale radius for nfw model", "1.0" },

  { "X_lens_min", NULL, 0, G_OPTION_ARG_DOUBLE, &X_lens_min_in, "lower X value for the grid", "-1.0" },
  { "X_lens_max", NULL, 0, G_OPTION_ARG_DOUBLE, &X_lens_max_in, "upper X value for the grid", "1.0" },
  { "X_divisions", NULL, 0, G_OPTION_ARG_INT, &X_divisions_in, "number of divisions for the x axis", "10" },

  { "Y_lens_min", NULL, 0, G_OPTION_ARG_DOUBLE, &Y_lens_min_in, "lower X value for the grid", "-1.0" },
  { "Y_lens_max", NULL, 0, G_OPTION_ARG_DOUBLE, &Y_lens_max_in, "upper X value for the grid", "1.0" },
  { "Y_divisions", NULL, 0, G_OPTION_ARG_INT, &Y_divisions_in, "number of divisions for the y axis", "10" },


  { NULL }
};


void potential_derivatives_enfw(double kappa_s, double r_s, double a, double b, double X_lens_min, double X_lens_max, int X_divisions, double Y_lens_min, double Y_lens_max, int Y_divisions){

  double conv_params[] = {kappa_s, r_s};
  double X_step = fabs(X_lens_max-X_lens_min)/double(X_divisions);
  double Y_step = fabs(Y_lens_max-Y_lens_min)/double(Y_divisions);
  double X = X_lens_min;
  double Y = Y_lens_min;

  double p1, p2, p11, p22, p12;

  for(int ix=0;ix<=X_divisions;ix++){
    for(int iy=0;iy<=Y_divisions;iy++){
      p1  = potential_1 (conv_nfw_circ, conv_params, X, Y, a,b);
      p2  = potential_2 (conv_nfw_circ, conv_params, X, Y, a,b);
      p11 = potential_11(conv_nfw_circ, conv_nfw_circ_prime, conv_params, X, Y, a,b);
      p22 = potential_22(conv_nfw_circ, conv_nfw_circ_prime, conv_params, X, Y, a,b);
      p12 = potential_12(conv_nfw_circ_prime, conv_params, X, Y, a,b);

      p1 = isnan(p1) || isinf(p1) ? 0.0 : p1;
      p2 = isnan(p2) || isinf(p2) ? 0.0 : p2;
      p11 = isnan(p11) || isinf(p11) ? 0.0 : p11;
      p22 = isnan(p22) || isinf(p22) ? 0.0 : p22;
      p12 = isnan(p12) || isinf(p12) ? 0.0 : p12;

      printf("%E %E %E %E %E %E %E\n", X, Y, p1, p2, p11, p22, p12);

      Y += Y_step;
    }
    X += X_step;
    Y = Y_lens_min;
  }
}

int main(int argc, char *argv[]){

  GError *error = NULL;
  GOptionContext *context;

  context = g_option_context_new ("- compute quantities related to the elliptical models");

  g_option_context_add_main_entries (context, entries, NULL);

  g_option_context_add_group (context,  grid_group ( &Lens_Plane_Grid_Entries ) );


 // g_option_context_add_main_entries (context, grid_group, NULL);
  g_option_context_parse (context, &argc, &argv, &error);

  printf("Running ellip_mass\n");
  double a = 1.0;
  double b = 0.5;




  if(compute_enfw_potential_derivatives_in){
    printf("NFW parameters: kappa_s=%f r_s=%f\n", kappa_s_in, r_s_in);
    potential_derivatives_enfw(kappa_s_in, r_s_in, a, b, X_lens_min_in , X_lens_max_in , X_divisions_in, Y_lens_min_in, Y_lens_max_in, Y_divisions_in);
  }


  return 0;
}
