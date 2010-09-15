/* Control file to gravlens source (to considerer SIS model) 
Started by H.Dumet and MSS Gill on  Sept 2, 2010
Commented and small changes by MSSG
 */
// Needed headers, in the g++ vs. gcc style
#include <cstdio>
#include <cmath> 
#include <ctime>
#include <iostream>

double fpar[10]; 
double par_var[10];
double ccpar[10];
double x0_l,y0_l,x0_s,y0_s; // coordinate
double a_lo,a_hi;
double eta_s,PA_s,R0; // source characteristics 
//**********************************  Begin main code ******************************/
//! Function to make gravlens config files
/*!
 \param ccpar[1]  gridlo1's value
 \return NOTHING
*/

double nullfunction();

int main(){

ccpar[1]=1.E-5; //  gridlo1's value
ccpar[2]=20.0;  //  gridhi1's value
ccpar[3]=1.E-8; // xtol's value
ccpar[4]=1.E-8; // crittol's value
ccpar[5]=1.E-8 ;// inttol's value

ccpar[6]=6; // maxlev's value
ccpar[7]=6 ;// imglev's value


x0_l=0.0, y0_l=0.0; //change later (lens coordinates)

/////////// Print to screen, just to check
for(int i=1;i<=10;i++){
    fpar[i]=0.;
    par_var[i]=0;
// printf("%i %E %E\n",i,fpar[i],par_var[i]);
  }
fpar[1]=1.0, fpar[2]=x0_l, fpar[3]=y0_l;
fpar[8]=0.0,fpar[10]=1.0;
for(int i=1;i<=10;i++){
printf("%i %E %E %E\n",i,fpar[i],par_var[i],ccpar[i]);
}

/////////////// Open lens (potential-only) config file for output
FILE *outfile = fopen ("lens_file.txt" , "w");

// creating the control file for the lens
fprintf(outfile," gridmode  %i\n",2);
fprintf(outfile," set ngrid1 = %i\n",25);
fprintf(outfile," set ngrid2 = %i\n",75);
fprintf(outfile," set gridlo1 = %E\n",ccpar[1]); 
fprintf(outfile," set gridhi1 = %E\n",ccpar[2]); 
fprintf(outfile," set xtol = %E\n",ccpar[3]); 
fprintf(outfile," set crittol = %E\n",ccpar[4]); 
fprintf(outfile," set inttol = %E\n",ccpar[5]); 
fprintf(outfile," set maxlev = %E\n",ccpar[6]); 
fprintf(outfile," set imglev = %E\n",ccpar[7]); 
fprintf(outfile," startup  %i %i\n",1,1); 
fprintf(outfile," alpha  %E %E %E %E %E %E %E %E %E %E\n",fpar[1],fpar[2],fpar[3],fpar[4],fpar[5],fpar[6],fpar[7],fpar[8],fpar[9],fpar[10]); 
fprintf(outfile,  "%E %E %E %E %E %E %E %E %E %E\n",par_var[1],par_var[2],par_var[3],par_var[4],par_var[5],par_var[6],par_var[7],par_var[8],par_var[9],par_var[10]); 
fprintf(outfile," plotcrit sis_curves.dat \n"); 


/////////////// Open source config file for output (to make the arcs)
int steps=1, nang=400;
FILE *outfile2 = fopen ("source_file.txt" , "w");
x0_s = 1.0/10.0, y0_s = 0 ; 
R0 = 1.0/15.0, eta_s=0.0, PA_s=0.0 ;
a_lo=R0,a_hi=R0;

// creating the source file
fprintf(outfile2," gridmode  %i\n",1);
fprintf(outfile2," set ngrid1 = %i\n",25);
fprintf(outfile2," set ngrid2 = %i\n",75);
fprintf(outfile2," set gridlo1 = %E\n",ccpar[1]); 
fprintf(outfile2," set gridhi1 = %E\n",ccpar[2]); 
fprintf(outfile2," set xtol = %E\n",ccpar[3]); 
fprintf(outfile2," set crittol = %E\n",ccpar[4]); 
fprintf(outfile2," set inttol = %E\n",ccpar[5]); 
fprintf(outfile2," set maxlev = %E\n",ccpar[6]); 
fprintf(outfile2," set imglev = %E\n",ccpar[7]); 
fprintf(outfile2," startup  %i %i\n",1,1); 
fprintf(outfile2," alpha  %E %E %E %E %E %E %E %E %E %E\n",fpar[1],fpar[2],fpar[3],fpar[4],fpar[5],fpar[6],fpar[7],fpar[8],fpar[9],fpar[10]); 
fprintf(outfile2,  "%E %E %E %E %E %E %E %E %E %E\n",par_var[1],par_var[2],par_var[3],par_var[4],par_var[5],par_var[6],par_var[7],par_var[8],par_var[9],par_var[10]);
fprintf(outfile2," ellsrc src_img_plot.txt  %E %E %E %E %E %E %i %i\n",x0_s,y0_s,eta_s,PA_s,a_lo,a_hi,steps,nang);
fprintf(outfile2," findimg  %E %E\n", x0_s,y0_s,"src_im_pos.txt");


// creating the file for the source

return 0;

}
