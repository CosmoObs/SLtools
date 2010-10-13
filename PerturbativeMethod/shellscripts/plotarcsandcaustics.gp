# MSSG
# Sept 2010

# set terminal postscript color
set pointsize 0.4

##### Alard 2008, Fig 8, Caustics:

# set out 'arcs.forAn.EllipCentralLens.and.OnePerturber.ps'

# plot 'arcs.EllipCentralLens.PerturberWithTheta_pEqualsZero.SourceWith.Xeq0.15.Yeq0.R_Seq0.025.SourceEllipEqZero.dat' u 1:2

pause -1

# set out 'tangential.caustics.forAn.EllipCentralLens.and.OnePerturber.ps'

# plot  "tang_caustic.EllipCentralLens.NoPerturber.dat" u 1:2,"tang_caustic.EllipCentralLens.PerturberWithTheta_pEqualsZero.dat" u 1:2,"tang_caustic.EllipCentralLens.PerturberWithTheta_pEqualsPiOverTen.dat" u 1:2


##### Alard 2008, Fig 9, Arcs: 

# set out 'arcs.forAn.EllipCentralLens.and.OnePerturber.ps'

# plot "tang_caust.dat" u 1:2, "tang_crit.dat" u 1:2, "src_plot.dat" u 1:2, "arcs.EllipCentralLens.PerturberWithTheta_pEqualsZero.SourceWith.Xeq0.15.Yeq0.R_Seq0.025.SourceEllipEqZero.dat" u 1:2 ti "Lensed Arcs for Ellip Lens and Perturber With m_p = 0.03 of central lens, Source at (0.15,0), R_S=0.025", "./tang_crit.dat" u 1:2 ti "Tangential Crit Curve For This Potential"
 
# pause -1


# plot "./arcs_sis.dat" u 1:2

 set out 'arcs.forAn.EllipCentralLens.and.OnePerturber.TwoSources.ps'

 plot "./arcs_sis.centerleft.dat" u 1:2,"./arcs_sis.centerright.dat" u 1:2

# plot  "tang_caustic.EllipCentralLens.NoPerturber.dat" u 1:2

 pause -1

# plot "tang_caust.rp0.8.dat" u 1:2, "tang_caust.rp0.9.dat" u 1:2,"tang_caust.rp1.0.dat" u 1:2,"tang_caust.rp1.1.dat" u 1:2

# pause -1

 set terminal X11
