# Contents of slcode/PerturbativeMethod dir
# Started by MSSG, Sept 8, 2010

---- makearcs_pertmethod.cpp -- first code that made arcs successfully
 to compare vs. MMakler Math'ca file

To get this to work and show the arcs, compile with g++, run
executable and dump values to a text output file, then plot these
(without connecting lines is best) in either e.g. Grace or Gnuplot.

---- compare_gravlensTOpertarcs.sh - compiles the two codes and runs
them, then loads the output files into Grace; gravlens installed in
the same directory for now,


---- configfile_maker_gravlens.cpp*--  Code to compare gravlens to pertarcs output.


---- Recent code added by Gabriel, to do modifications of actually
     perturbed potential (changed versions of above):

simple_example_old.cpp*
simple_example.cpp*
perturbative_method.h*
test.cpp*
theta.h
