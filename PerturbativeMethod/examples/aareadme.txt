# Sept 27, 2010

---------------------------- Code contents of this dir:



---- perturbedSIS_with_ellipticalSource.cpp* -- Doxygenated -- this
     incorporates an off-axis perturber, and takes 9 arguments in the
     command line corresponding to source and perturber properties
     (start: 9/2010)



---------------------------- How to run things:

--- To make caustic curves of Fig 8, do:

   ./a.out 5000 0 0 0 0 0 0 0 0 0.1
   mv tang_caust.dat tang_caustic.EllipCentralLens.NoPerturber.dat
		
   ./a.out 5000 0 0 0 0 0 0.03 1.3 0 0.1
   mv tang_caust.dat tang_caustic.EllipCentralLens.PerturberWithTheta_pEqualsZero.dat
  
  ./a.out 5000 0 0 0 0 0 0.03 1.3 .31415 0.1
  mv tang_caust.dat tang_caustic.EllipCentralLens.PerturberWithTheta_pEqualsPiOverTen.dat


Then run  plotarcsandcaustics.gp*


--- To make arcs of Fig 9 do:


./a.out 5000 0.15 0 0.025 0 0 0.03 1.3 0 0.1
	v arcs_sis.dat arcs.EllipCentralLens.PerturberWithTheta_pEqualsZero.SourceWith.Xeq0.15.Yeq0.R_Seq0.025.SourceEllipEqZero.dat
