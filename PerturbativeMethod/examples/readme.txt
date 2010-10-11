# Sept 27, 2010
# Current: Oct 5, 2010

----------- This is the examples subdir, for code to make arcs,
            caustics, crit curves in the forward direction

----------- Using:  perturbedSIE_with_ellipticalSource.cpp

------- Compile with:  g++ perturbedSIE_with_ellipticalSource.cpp

--- Then to make caustic curves of Alard 2008, Fig 8, do:

(First is num of pts, next 5 are source params, next 3 are perturber properties, last is central lens ellipp) 

   ./a.out 5000   0 0 0 0 0   0 0 0 0.1
   mv tang_caust.dat tang_caustic.EllipCentralLens.NoPerturber.dat
		
   ./a.out 5000   0 0 0 0 0   0.03 1.3 0 0.1
   mv tang_caust.dat tang_caustic.EllipCentralLens.PerturberWithTheta_pEqualsZero.dat
  
  ./a.out 5000    0 0 0 0 0    0.03 1.3 .31415 0.1
  mv tang_caust.dat tang_caustic.EllipCentralLens.PerturberWithTheta_pEqualsPiOverTen.dat


Then run the correct line of  plotarcsandcaustics.gp*


--- To make arcs of Fig 9 do:


./a.out 5000 0.15 0 0.025 0 0 0.03 1.3 0 0.1
	
mv arcs_sis.dat arcs.EllipCentralLens.PerturberWithTheta_pEqualsZero.SourceWith.Xeq0.15.Yeq0.R_Seq0.025.SourceEllipEqZero.dat
