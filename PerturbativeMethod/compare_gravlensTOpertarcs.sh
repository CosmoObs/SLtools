#!/bin/bash

################################### Shellscript to execute all steps to compare Gravlens output to Pertarcs (Alard method) output

########################## Compile code
echo "Compiling the code"
######## First the code to make the gravlens config file
g++ -o glconfigmaker comp_gravlens/configfile_maker_gravlens.cpp
######## Then the pertarcs code 
g++ -o pertarcs makearcs_pertmethod.cpp

######################### Run the pertarcs code to make the output text files which contain the arc points in the image plane
### (This contains an image of the critical curve at R_E automatically)
echo "Running Gabriel´s  code"
./pertarcs > pert_arcs.txt

####### Optionally look at these, separately
# xmgrace -view 0.15 0.15 0.85 0.85 pert_arcs.txt &

######################## Run the code to make the gravlens config file
echo "Running the Gravlens software "
./glconfigmaker

######################### Run gravlens first over just the potential file to get the critical curve at R_E in the image plane
./gravlens lens_file.txt

######################### Run gravlens next over just the source file to get the arcs made in the image plane
./gravlens source_file.txt

######################### Load both the gravlens and pertarcs images into grace to compare them
xmgrace -view 0.15 0.15 0.85 0.85 sis_curves.dat src_img_plot.txt pert_arcs.txt