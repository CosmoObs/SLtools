#!/bin/sh

# MSSG
# Sept 2010

# To vary various arc params, then output the results for 20
# variations into 20 files under the arcdata subdir

# Compile code, in case we need to
g++ perturbedSIS_with_ellipticalSource.cpp -o make.arcdata.output
           
# Start loop over 20 variations
for i  in $(seq 20);
do

# This is how to make an integer into a floating pt num in a shellscript
ivar=$(echo "scale=1; $i/10" |bc)
# i2=$(printf '%.3f \n' $i1)
echo $i $ivar
#                      npts  x0  y0 R0  eta   theta_s    mp  rp theta_p
 ./make.arcdata.output 1000  1   0  0.1 0.1  0           0.5 0.5 $ivar    > arcdata/arcs$i.txt 

done
