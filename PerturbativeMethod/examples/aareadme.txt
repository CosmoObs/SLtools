# Sept 27, 2010
# Current: Oct 11, 2010
# MSSG

----------- This is the examples subdir, for code to make arcs,
            caustics, crit curves in the forward direction for now,
            from all the primary code in the main dir

------------- Subdir Contains:

-external_shear_teste.cpp*
-makearcs_pertmethod.cpp*
-pnfw_test.cpp*
-simple_example_old.cpp*
-theta_example.cpp*

-configfile_maker_gravlens.cpp*
-perturbed.PE_SIS_with_ellipticalSource.cpp*

-plotarcsandcaustics.gp*




------------- Code Purposes:

-- perturbed.PE_SIS_with_ellipticalSource.cpp - Deprecated, partly Doxygenated - To make
            arcs, caustics, crit curves for the following config:
            Central pseudo-elliptical SIS pot with one circular SIS
            perturber, and accepts an elliptical source with
            orientation along x-axis.  This code is now superceded by
            sis_sub_model.h* in the main dir, but was the original one
            i (MSSG) used for making all these curves, so kept in the
            dir for now.

-- external_shear_teste.cpp* - unDoxygenated - Code by HSDM to test the external shr model

-- pnfw_test.cpp* - unDoxygenated - Code by HSDM to test the PNFW model

-- makearcs_pertmethod.cpp* - Deprecated, unDoxygenated - Very orig code to make the
   arcs with a simple central SIS pot, successfully compared vs. MMakler Math'ca file.
   Just a hardcoded SIS central potential with an off-axis circular source.

-- configfile_maker_gravlens.cpp* - unDoxygenated - Orig file, makes
   the config file for gravlens.  To run see
   compare_gravlensTOpertarcs.sh in shellscripts subdir

-- simple_example_old.cpp* - Deprecated, unDoxygenated - Slightly more complex
   version of makearcs_pertmethod.cpp* to make the arcs with a simple
   central SIS pot and an elliptical source

-- theta_example.cpp* - unDoxygenated - Code by GBC (?) to exercise the theta_find
  function in main dir that gives the minima and maxima in theta of
  the arcs made by a certain potential config




-plotarcsandcaustics.gp* - Gnuplot macro to plot up the various created files


--------------------------- To use:  


-------------- perturbed.PE_SIS_with_ellipticalSource.cpp

------- Compile with:  g++ perturbed.PE_SIS_with_ellipticalSource.cpp -o arcmaker


--- Then to make caustic curves of Alard 2008, Fig 8, use these params

 First is num of pts, next 5 are source params: x0   y0  R0 (size)  eta_s (ellip)   theta_s (ang. position)
 Next 3 are perturber properties:               mp   rp  theta_p
 Last is central lens ellip:                    eta 


   ./arcmaker 5000   0 0 0 0 0   0 0 0 0.1
   mv tang_caust.dat tang_caustic.EllipCentralLens.NoPerturber.dat
		
   ./arcmaker 5000   0 0 0 0 0   0.03 1.3 0 0.1
   mv tang_caust.dat tang_caustic.EllipCentralLens.PerturberWithTheta_pEqualsZero.dat
  
  ./arcmaker 5000    0 0 0 0 0    0.03 1.3 .31415 0.1
  mv tang_caust.dat tang_caustic.EllipCentralLens.PerturberWithTheta_pEqualsPiOverTen.dat


Then run the correct line of  plotarcsandcaustics.gp*

--- To make arcs of Fig 9 do:


./arcmaker 5000 0.15 0 0.025 0 0 0.03 1.3 0 0.1
	
mv arcs_sis.dat arcs.EllipCentralLens.PerturberWithTheta_pEqualsZero.SourceWith.Xeq0.15.Yeq0.R_Seq0.025.SourceEllipEqZero.dat

--- See also  vary_arc_params.sh* in shellscripts subdir


------------------ makearcs_pertmethod.cpp 



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

