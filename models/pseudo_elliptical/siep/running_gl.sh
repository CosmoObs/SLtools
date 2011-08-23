#/bin/bash

# script to runing the gravlens and separate its output file for the SIEP.

g++ -o cgr configfile_maker_gravlens.cpp
g++ -o rgl reading_gl.cpp

./cgr

gravlens lens_file.txt

./rgl


