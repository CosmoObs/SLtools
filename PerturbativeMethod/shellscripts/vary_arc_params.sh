#!/bin/sh

# MSSG
# Sept 2010

# To vary various arc params, then output the results for 20
# variations into 20 files under the arcdata subdir

# Compile code, in case we need to
g++ ../examples/perturbed.PE_SIS_with_ellipticalSource.cpp -o arcmaker
           
# Start loop over 20 variations
for i  in $(seq 20);
do

# This is how to make an integer into a floating pt num in a shellscript -- i think this must be done in bash shell
ivar=$(echo "scale=1; $i/10" |bc)  
# i2=$(printf '%.3f \n' $i1)
echo $i $ivar
#           npts  x0    y0 R0   eta_s   theta_s       mp  rp     theta_p eta_central
 ./arcmaker 1000  0.2   0  0.05 0.1     0             0.5 $ivar  0       0
mv arcs_sis.dat  arcdata/arcs$i.txt 
mv src_plot.dat  arcdata/srcplot$i.txt 
mv tang_caust.dat  arcdata/caustic$i.txt 
mv tang_crit.dat  arcdata/critcurve$i.txt 


done
