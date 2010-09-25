# Contents of slcode/PerturbativeMethod dir
# Started by MSSG, Sept 8, 2010

[For making Doxygen files in the html and latex subdirs, remember to
first generate Doxyfile with doxygen -g in the dir, then run doxygen
in the dir]

---- The below is a copy of this page:

http://twiki.on.br/bin/save/StrongLensing/PerturbativeArcContoursv01


---------------------------------------------------------------------------------------------------------
---------------------------------------------------------------------------------------------------------


--------------------- Pert Method code in dir: .../repositories/sltools/PerturbativeMethod

---- makearcs_pertmethod.cpp -- unDoxygenated, deprecated -- first
 code that made arcs successfully to compare vs. MMakler Math'ca file
 (start: 8/2010)


Just a hardcoded SIS central potential with an off-axis circular source.

To get this to work and show the arcs, compile with g++, run
executable and dump values to a text output file, then plot these
(without connecting lines is best) in either e.g. Grace or Gnuplot.

---- perturbedSIS_with_ellipticalSource.cpp* -- Doxygenated -- this
     incorporates an off-axis perturber, and takes 9 arguments in the
     command line corresponding to source and perturber properties
     (start: 9/2010)

This is driven by vary_arc_params.sh* , a shellscript which will
output a bunch (20 by default) of ascii files with differing input
params, then you can load these one by one into gnuplot, with the
script plotarcs.gp (e.g. just run from command line: gnuplot
plotarcs.gp)


--------------------- Comparing to Gravlens

---- compare_gravlensTOpertarcs.sh - compiles:

 configfile_maker_gravlens.cpp -- Doxygenated -- which makes the config file for gravlens

 makearcs_pertmethod.cpp -- code that makes arcs in pert method, as discussed above

Then runs these, then loads the output files into Grace.  Needs
gravlens installed in the same directory for now, though this is not
committed to the repo (start: 8/2010)

------------------------ Modularized code

----- Based on our initial codes above, we began modularizing pieces
      as we expanded.  To wit:

---- lens_models.h* -- Doxygenated -- This takes out the f1,df0dt
     functions so they can be accessed directly from any code here

---- theta_find.h* -- unDoxygenated yet -- code to find the extrema of
     arcs by checking for sign changes of the radicand in the x
     solution (the distance to arc from origin in lens plane).  Writes
     the theta values for start and end of arc into an output txt
     file.

---- pnfw_test.cpp* -- values of shear, conv, and needed deriv values
     in the PNFW model

----- perturbative_method.h* -- general accessor utility code for many needed functions:
 - central potential phi
 - perturber potential psi 
 - f function of needed derivs
 - elliptical source params
 - f1bar
 - df0dthetabar
 - radicand  (arg_sqrt)
 - x+, x- solns ( dr_plus,  dr_minus)
 - rcrit
 - Cartesian coords for caustic
 - Cartesian coords for source 


--------------------- Writeups subdir

Tex notes for Alard 2007 and 2008 papers, and a note for PNFW models.


---------------------- Examples subdir:

Codes that make use of the more general codes in the above dir to make arcs 
