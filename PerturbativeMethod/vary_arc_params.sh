#!/bin/sh

g++ perturbedSIS_with_ellipticalSource.cpp 
           
for i  in $(seq 20);
do

ivar=$(echo "scale=1; $i/10" |bc)
# i2=$(printf '%.3f \n' $i1)
echo $i $ivar
#        npts  x0  y0 R0  eta theta_s    mp  rp theta_p
 ./a.out 1000  0.1 0  0.1 0 0.0        0.6 $ivar 0   > arcs$i.txt 

done
