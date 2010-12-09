#include <cstdio>
#include <cmath>
#include <algorithm>
#include <cstdlib>

#include "../../numerical_methods/hausdorff_distance.h"


#include "../perturbative_method.h"
#include "../pnfw_model.h"
#include "../nfw_circular_model.h"

// to compile make: g++ -Wall -o exec_h hausdorff.cpp `pkg-config gsl --cflags --libs`

double generate_curves_f(double ks,double rs,double elp, int npt, double cc_x[], double cc_y[], double ca_x[], double ca_y[]){
//Escreve no arquivo in_penf_par.txt os imputs para o programa em fortran.
  FILE * h_in = fopen ("../../lensing_models/pnfw/in_pnfw_par.txt","w");

  double lw_u = 10.0;
  int iflag=1;

  fprintf(h_in,"%f %f %f %f %i %i",ks, rs, elp, lw_u, iflag, npt); 
  fclose(h_in);

  system("sh hausdorff.sh");


//Read the output from fortran program and put in cc_x[] and cc_y[]
  FILE *cc_file = fopen("../../lensing_models/pnfw/fort.65","r");
  FILE *ca_file = fopen("../../lensing_models/pnfw/fort.75","r");
  int number_columns = 500;  //number of coluns (characteres)
  char buffer[number_columns]; //variavel a receber as linhas do arquivo


  int i=0;
  while (!feof(cc_file)){
    fgets (buffer , number_columns , cc_file);
    sscanf(buffer ,"%lg %lg", &cc_x[i] ,&cc_y[i]);
    fgets (buffer , number_columns , ca_file);
    sscanf(buffer ,"%lg %lg", &ca_x[i] ,&ca_y[i]);
//    printf("%i %f %f %f %f \n",i,cc_x[i],cc_y[i],ca_x[i],ca_y[i]);
    i++;
  }

  fclose(cc_file);
  fclose(ca_file);

  return 0.0;
}


double generate_curves_p(double ks,double rs,double elp, int npt, double cc_x_per[], double cc_y_per[], double ca_x_per[], double ca_y_per[]){

  double cc_r_per[4*npt];

  double pot_params[] = {ks,rs};
  double _r_e_nfw= r_e_nfw_find(pot_params);
  double kappa_2=kappa2_nfw(pot_params,_r_e_nfw);  

  double pert_params[] = {ks,rs,elp};


  double theta = 0.0;
  for(int i=0;i<4*npt;i++){

    cc_r_per[i] = r_crit(f1_pnfw, D2f0Dtheta2_pnfw, kappa_2, theta, pert_params, _r_e_nfw);

    cc_x_per[i] = cc_r_per[i]*cos(theta);
    cc_y_per[i] = cc_r_per[i]*sin(theta);


    ca_x_per[i] = caustic_y1(Df0Dtheta_pnfw, D2f0Dtheta2_pnfw, theta, pert_params, _r_e_nfw);
    ca_y_per[i] = caustic_y2(Df0Dtheta_pnfw, D2f0Dtheta2_pnfw, theta, pert_params, _r_e_nfw);


    //printf("%E %E\n",theta, atan(ca_y_per[i]/ca_x_per[i]));
    theta += 6.2831853/double(4*npt);
  }
  return 0.0;
}


double mean_rad_vec(double vec_x[], double vec_y[], int npt){
  double out = 0.0;
  for(int i=0;i<npt;i++) out+=sqrt( pow(vec_x[i],2) + pow(vec_y[i],2) )/double(npt);
  return out;
}

void invar_test(){
  double ks = 1.1;
  double rs = 1.2;
  double elp= 0.1;

  int npt = 200;

  double cc_x[4*npt], cc_y[4*npt];
  double ca_x[4*npt], ca_y[4*npt];

  double cc_x_per[4*npt], cc_y_per[4*npt];
  double ca_x_per[4*npt], ca_y_per[4*npt];

  double r_cc_mean = 0.0;
  double r_ca_mean = 0.0;

  double r_cc_mean_per = 0.0;
  double r_ca_mean_per = 0.0;

  generate_curves_f(ks, rs, elp, npt, cc_x    , cc_y    , ca_x    , ca_y);
  generate_curves_p(ks, rs, elp, npt, cc_x_per, cc_y_per, ca_x_per, ca_y_per);
  
  r_cc_mean = mean_rad_vec(cc_x,cc_y,4*npt);
  r_ca_mean = mean_rad_vec(ca_x,ca_y,4*npt);;
  r_cc_mean_per = mean_rad_vec(cc_x_per,cc_y_per,4*npt);;
  r_ca_mean_per = mean_rad_vec(ca_x_per,ca_y_per,4*npt);;

  printf("\n\ndist2curves for rs=%.2f:\n",rs);
  printf("1- %E 1/rs- %E\n",dist2curves(cc_x, cc_y, npt, cc_x_per, cc_y_per, npt,1), dist2curves(cc_x, cc_y, npt, cc_x_per, cc_y_per, npt,1)/rs);
  printf("2- %E 2/rs- %E\n",dist2curves(cc_x, cc_y, npt, cc_x_per, cc_y_per, npt,2), dist2curves(cc_x, cc_y, npt, cc_x_per, cc_y_per, npt,2)/rs);
  printf("3- %E 3/rs- %E\n\n",dist2curves(cc_x, cc_y, npt, cc_x_per, cc_y_per, npt,3), dist2curves(cc_x, cc_y, npt, cc_x_per, cc_y_per, npt,3)/rs);

  printf("hd1- %E hd1/rs- %E\n",haus_dist_f4(cc_x, cc_y, npt, cc_x_per, cc_y_per, npt,1), haus_dist_f4(cc_x, cc_y, npt, cc_x_per, cc_y_per, npt,1)/rs);
  printf("hd2- %E hd2/rs- %E\n",haus_dist_f4(cc_x, cc_y, npt, cc_x_per, cc_y_per, npt,2), haus_dist_f4(cc_x, cc_y, npt, cc_x_per, cc_y_per, npt,2)/rs);
  printf("hd3- %E hd3/rs- %E\n\n",haus_dist_f4(cc_x, cc_y, npt, cc_x_per, cc_y_per, npt,3), haus_dist_f4(cc_x, cc_y, npt, cc_x_per, cc_y_per, npt,3)/rs);

  printf("<r_cc>-     %E <r_cc>/rs-     %E\n",r_cc_mean,r_cc_mean/rs);
  printf("<r_ca>-     %E <r_cc>/rs-     %E\n",r_ca_mean,r_ca_mean/rs);
  printf("<r_cc_per>- %E <r_cc_per>/rs- %E\n",r_cc_mean_per,r_cc_mean_per/rs);
  printf("<r_ca_per>- %E <r_cc_per>/rs- %E\n",r_ca_mean_per,r_ca_mean_per/rs);


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  rs = 0.17;
  generate_curves_f(ks, rs, elp, npt, cc_x    , cc_y    , ca_x    , ca_y);
  generate_curves_p(ks, rs, elp, npt, cc_x_per, cc_y_per, ca_x_per, ca_y_per);

  r_cc_mean = mean_rad_vec(cc_x,cc_y,4*npt);
  r_ca_mean = mean_rad_vec(ca_x,ca_y,4*npt);;
  r_cc_mean_per = mean_rad_vec(cc_x_per,cc_y_per,4*npt);;
  r_ca_mean_per = mean_rad_vec(ca_x_per,ca_y_per,4*npt);;

  printf("\n\ndist2curves for rs=%.2f:\n",rs);
  printf("1- %E 1/rs- %E\n",dist2curves(cc_x, cc_y, npt, cc_x_per, cc_y_per, npt,1), dist2curves(cc_x, cc_y, npt, cc_x_per, cc_y_per, npt,1)/rs);
  printf("2- %E 2/rs- %E\n",dist2curves(cc_x, cc_y, npt, cc_x_per, cc_y_per, npt,2), dist2curves(cc_x, cc_y, npt, cc_x_per, cc_y_per, npt,2)/rs);
  printf("3- %E 3/rs- %E\n\n",dist2curves(cc_x, cc_y, npt, cc_x_per, cc_y_per, npt,3), dist2curves(cc_x, cc_y, npt, cc_x_per, cc_y_per, npt,3)/rs);

  printf("hd1- %E hd1/rs- %E\n",haus_dist_f4(cc_x, cc_y, npt, cc_x_per, cc_y_per, npt,1), haus_dist_f4(cc_x, cc_y, npt, cc_x_per, cc_y_per, npt,1)/rs);
  printf("hd2- %E hd2/rs- %E\n",haus_dist_f4(cc_x, cc_y, npt, cc_x_per, cc_y_per, npt,2), haus_dist_f4(cc_x, cc_y, npt, cc_x_per, cc_y_per, npt,2)/rs);
  printf("hd3- %E hd3/rs- %E\n\n",haus_dist_f4(cc_x, cc_y, npt, cc_x_per, cc_y_per, npt,3), haus_dist_f4(cc_x, cc_y, npt, cc_x_per, cc_y_per, npt,3)/rs);

  printf("<r_cc>-     %E <r_cc>/rs-     %E\n",r_cc_mean,r_cc_mean/rs);
  printf("<r_ca>-     %E <r_cc>/rs-     %E\n",r_ca_mean,r_ca_mean/rs);
  printf("<r_cc_per>- %E <r_cc_per>/rs- %E\n",r_cc_mean_per,r_cc_mean_per/rs);
  printf("<r_ca_per>- %E <r_cc_per>/rs- %E\n",r_ca_mean_per,r_ca_mean_per/rs);


}



int main(){
  //invar_test();

//Escreve no arquivo in_penf_par.txt os imputs para o programa em fortran.
  //FILE * h_in = fopen ("../../lensing_models/pnfw/in_pnfw_par.txt","w");
  double ks = 0.08;
  double rs = 1.;
  double elp= 0.7;
  int npt = 200;

  double cc_x[4*npt], cc_y[4*npt];
  double ca_x[4*npt], ca_y[4*npt];

  double cc_x_per[4*npt], cc_y_per[4*npt];
  double ca_x_per[4*npt], ca_y_per[4*npt];

//  generate_curves_f(ks, rs, elp, npt, cc_x    , cc_y    , ca_x    , ca_y);
//  generate_curves_p(ks, rs, elp, npt, cc_x_per, cc_y_per, ca_x_per, ca_y_per);


//Mean radius for CC and CA from the fortran code
  double r_cc_mean = mean_rad_vec(cc_x,cc_y,4*npt);
  double r_ca_mean = mean_rad_vec(ca_x,ca_y,4*npt);;

  double r_cc_mean_per = mean_rad_vec(cc_x_per,cc_y_per,4*npt);;
  double r_ca_mean_per = mean_rad_vec(ca_x_per,ca_y_per,4*npt);;



  FILE *h_dist_cc = fopen("h_cc_rs1_elp07-3.dat","w");
  FILE *h_dist_ca = fopen("h_ca_rs1_elp07-3.dat","w");

  for(int i=0;i<50;i++){
    generate_curves_f(ks, rs, elp, npt, cc_x    , cc_y    , ca_x    , ca_y);
    generate_curves_p(ks, rs, elp, npt, cc_x_per, cc_y_per, ca_x_per, ca_y_per);

    r_cc_mean = mean_rad_vec(cc_x,cc_y,4*npt);
    r_ca_mean = mean_rad_vec(ca_x,ca_y,4*npt);;

    r_cc_mean_per = mean_rad_vec(cc_x_per,cc_y_per,4*npt);
    r_ca_mean_per = mean_rad_vec(ca_x_per,ca_y_per,4*npt);

    fprintf(h_dist_cc,"%E %E\n",ks,haus_dist_f4(cc_x_per, cc_y_per, 4*npt, cc_x, cc_y, 4*npt,3)/(r_cc_mean) );
    fprintf(h_dist_ca,"%E %E\n",ks,haus_dist_f4(ca_x_per, ca_y_per, 4*npt, ca_x, ca_y, 4*npt,3)/(r_ca_mean) );
    ks+=1.0/50.0;
  }

  fclose(h_dist_cc);
  fclose(h_dist_ca);

  printf("r_cc_mean=%E r_ca_mean=%E\n",r_cc_mean,r_ca_mean);
  printf("r_cc_mean_per=%E r_ca_mean_per=%E\n",r_cc_mean_per,r_ca_mean_per);

  return 0;
}
