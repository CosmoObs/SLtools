#include <cstdio>
#include <cmath>
#include <algorithm>
#include <cstdlib>

double dist2pts(double x1, double y1, double x2, double y2){return sqrt( pow(x1-x2,2) + pow(y1-y2,2)  );}

double dist_pt_curve(double ptx, double pty, double cx[], double cy[], int cnpts){
  double dist[cnpts];
  for(int i=0;i<cnpts;i++) dist[i] = dist2pts(ptx, pty, cx[i], cy[i]);
  return *std::min_element(dist,dist+cnpts);
}

double dist_max_pt_curve(double ptx, double pty, double cx[], double cy[], int cnpts){
  double dist[cnpts];
  for(int i=0;i<cnpts;i++) dist[i] = dist2pts(ptx, pty, cx[i], cy[i]);
  return *std::max_element(dist,dist+cnpts);
}

//type = 0 -> d3, type = 1 -> d1 , else -> d2
double dist2curves(double c1x[], double c1y[], int c1npts, double c2x[], double c2y[], int c2npts, int type = 0){

  double dist_vec[c1npts];
  double dist_out=0.0; 
  for(int i=0;i<c1npts;i++) dist_vec[i] = dist_pt_curve(c1x[i], c1y[i], c2x, c2y, c2npts);

  if(type == 0){
    for(int i=0;i<c1npts;i++) dist_out += dist_vec[i]/c1npts;
  } else if(type == 1) {
    dist_out = *std::min_element(dist_vec,dist_vec+c1npts);;
  } else {
    dist_out = *std::max_element(dist_vec,dist_vec+c1npts);;
  }

  return dist_out;
}


double dist2curves_max(double c1x[], double c1y[], int c1npts, double c2x[], double c2y[], int c2npts){

  double dist_vec[c1npts];
  double dist_out=0.0; 
  for(int i=0;i<c1npts;i++) dist_vec[i] = dist_max_pt_curve(c1x[i], c1y[i], c2x, c2y, c2npts);

  dist_out = *std::max_element(dist_vec,dist_vec+c1npts);;


  return dist_out;
}
int main(){


  int npts1 = 1000; //number of poits for ellipse 1
  int npts2 = 1000; //number of poits for ellipse 2
  double c1x[npts1],c1y[npts1];//vector for curve 1
  double c2x[npts2],c2y[npts2];//vector for curve 2
  double a1,e1; //major semi axis and ellipticity for ellipse 1
  double a2,e2; //major semi axis and ellipticity for ellipse 2
  double x02,y02;//center of the ellipse 2

  a1=2.0;
  e1 = 0.;

  a2 =2.25;
  e2 = 0.5;
  x02 = 0.1;
  y02 = 0.2;

  double t=0;
  for(int i=0;i<npts1;i++){//defining curve 1
    c1x[i] = a1*cos(t);
    c1y[i] = a1*sqrt(1.0-e1)*sin(t);
    t+=6.28319/double(npts1);
  }
  t=0.0;
  for(int i=0;i<npts2;i++){//defining curve 2
    c2x[i] = x02 + a2*cos(t);
    c2y[i] = y02 + a2*sqrt(1.0-e2)*sin(t);
    t+=6.28319/double(npts2);
  }


/********************************************************************************************************************/
  printf("Distancies from curve A to curve B\n");
  printf("  d1    : %f\n", dist2curves(c1x,c1y,npts1,c2x,c2y,npts2,1) );
  printf("  d2    : %f\n", dist2curves(c1x,c1y,npts1,c2x,c2y,npts2,2) );
  printf("  d3    : %f\n", dist2curves(c1x,c1y,npts1,c2x,c2y,npts2,0) );
  printf("  d4    : %f\n", dist2curves_max(c1x,c1y,npts1,c2x,c2y,npts2) );
  printf("  d3/d2 : %f\n", dist2curves(c1x,c1y,npts1,c2x,c2y,npts2,0)/dist2curves(c1x,c1y,npts1,c2x,c2y,npts2,2) );
  printf("  d3/d4 : %f\n", dist2curves(c1x,c1y,npts1,c2x,c2y,npts2,0)/dist2curves_max(c1x,c1y,npts1,c2x,c2y,npts2) );

  printf("Distancies from curve B to curve A\n");
  printf("  d1    : %f\n", dist2curves(c2x,c2y,npts2,c1x,c1y,npts1,1) );
  printf("  d2    : %f\n", dist2curves(c2x,c2y,npts2,c1x,c1y,npts1,2) );
  printf("  d3    : %f\n", dist2curves(c2x,c2y,npts2,c1x,c1y,npts1,0) );
  printf("  d4    : %f\n", dist2curves_max(c1x,c1y,npts1,c2x,c2y,npts2) );
  printf("  d3/d2 : %f\n", dist2curves(c1x,c1y,npts1,c2x,c2y,npts2,0)/dist2curves(c1x,c1y,npts1,c2x,c2y,npts2,2) );
  printf("  d3/d4 : %f\n", dist2curves(c1x,c1y,npts1,c2x,c2y,npts2,0)/dist2curves_max(c1x,c1y,npts1,c2x,c2y,npts2) );



  FILE * c1File = fopen ("c1.dat","w");
  for(int i=0;i<npts1;i++){ fprintf(c1File,"%f  %f\n",c1x[i], c1y[i]); }
  FILE * c2File = fopen ("c2.dat","w");
  for(int i=0;i<npts2;i++){ fprintf(c2File,"%f  %f\n",c2x[i], c2y[i]); }
  fclose(c1File);
  fclose(c2File);
  system("xmgrace c1.dat c2.dat &");
/********************************************************************************************************************/






/********************************************************************************************************************/
  double factor = 2.0;
  printf("\n\nEnlarging curves by a factor %f\n",factor);
  for(int i=0;i<npts1;i++){ c1x[i]=factor*c1x[i]; c1y[i]=factor*c1y[i]; }
  for(int i=0;i<npts2;i++){ c2x[i]=factor*c2x[i]; c2y[i]=factor*c2y[i]; }
  printf("Distancies from curve A to curve B\n");
  printf("  d1    : %f\n", dist2curves(c1x,c1y,npts1,c2x,c2y,npts2,1) );
  printf("  d2    : %f\n", dist2curves(c1x,c1y,npts1,c2x,c2y,npts2,2) );
  printf("  d3    : %f\n", dist2curves(c1x,c1y,npts1,c2x,c2y,npts2,0) );
  printf("  d4    : %f\n", dist2curves_max(c1x,c1y,npts1,c2x,c2y,npts2) );
  printf("  d3/d2 : %f\n", dist2curves(c1x,c1y,npts1,c2x,c2y,npts2,0)/dist2curves(c1x,c1y,npts1,c2x,c2y,npts2,2) );
  printf("  d3/d4 : %f\n", dist2curves(c1x,c1y,npts1,c2x,c2y,npts2,0)/dist2curves_max(c1x,c1y,npts1,c2x,c2y,npts2) );

  printf("Distancies from curve B to curve A\n");
  printf("  d1    : %f\n", dist2curves(c2x,c2y,npts2,c1x,c1y,npts1,1) );
  printf("  d2    : %f\n", dist2curves(c2x,c2y,npts2,c1x,c1y,npts1,2) );
  printf("  d3    : %f\n", dist2curves(c2x,c2y,npts2,c1x,c1y,npts1,0) );
  printf("  d4    : %f\n", dist2curves_max(c1x,c1y,npts1,c2x,c2y,npts2) );
  printf("  d3/d2 : %f\n", dist2curves(c1x,c1y,npts1,c2x,c2y,npts2,0)/dist2curves(c1x,c1y,npts1,c2x,c2y,npts2,2) );
  printf("  d3/d4 : %f\n", dist2curves(c1x,c1y,npts1,c2x,c2y,npts2,0)/dist2curves_max(c1x,c1y,npts1,c2x,c2y,npts2) );
  FILE * c12File = fopen ("c1.dat","w");
  for(int i=0;i<npts1;i++){ fprintf(c12File,"%f  %f\n",c1x[i], c1y[i]); }
  FILE * c22File = fopen ("c2.dat","w");
  for(int i=0;i<npts2;i++){ fprintf(c22File,"%f  %f\n",c2x[i], c2y[i]); }
  fclose(c12File);
  fclose(c22File);
  system("xmgrace c1.dat c2.dat &");
/********************************************************************************************************************/

  return 0;
}
