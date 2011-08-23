# Sept 29,2010
# Last changes: MSSG

New code:

arcs.agr*
arquivosparacomparaocomogravlens.zip*
gracefile.cfg*
reading_gl.cpp*

------------------ Making arcs 

--- vary_arc_params.sh* , a shellscript which will output a bunch (20
by default) of ascii files with differing input params, then you can
load these one by one into gnuplot, with the script plotarcs.gp
(e.g. just run from command line: gnuplot plotarcs.gp)

Other gnuplot scripts:

plotarcsandcaustics.gp*
plotcaustic.gp*
plotcritcurves.gp*


--------------------- Comparing to Gravlens

---- compare_gravlensTOpertarcs.sh - compiles (in the dir above):

 configfile_maker_gravlens.cpp -- Doxygenated -- which makes the config file for gravlens

 makearcs_pertmethod.cpp -- code that makes arcs in pert method, as discussed above

Then runs these, then loads the output files into Grace.  Needs
gravlens installed in the same directory for now, though this is not
committed to the repo (start: 8/2010)
