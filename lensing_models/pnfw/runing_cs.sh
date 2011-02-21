#!/bin/bash
# this script perform the calculation for the arc cross section
# compiling the code
echo "compiling the code (waiting please)"
g77-3.4 -Wall -o dcs_pnfw arc_cross_section.f root_finding.f lensing_functions.f arc_cross_subroutine.f 
echo "the compilation is done"


#read(1,*) iflag, ks, e, lw_u, muth
iflag=4

# echo "enter the characterstic convergence"
# read ks1

echo "enter the value of the ellipticity"
read el1

ks1=0.5
lwu=10
muth=1.0
#el1=0.0
eflag=3
echo "enter the minimum magnification is :" $muth, "length-to-width ratio :" $lwu

echo $iflag $ks1 $el1 $lwu $muth $eflag> ifile

echo "scaling the arc cross section with ellipticity"
qumin=1.1
qumax=25.1
dqu=100
echo $qumin $qumax $dqu > ifile2

./dcs_pnfw

