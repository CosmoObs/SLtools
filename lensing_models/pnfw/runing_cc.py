# -*- coding: utf-8 -*-
import os 
from sltools.gravlens.find_CC_new import run_find_CC as run

ks=2.0
rs=1.0
el=0.5
thetal=90
#	
lw=10
iflag=3
eflag=3
npt=251

outpar = open('in_pnfw_par.txt','w')
outpar.write(" %f %f %f %f %i %i %i\n " %(ks,rs,el,lw,iflag,npt,eflag))
outpar.close()


print 'compiling the code (waiting please)'
os.system("gfortran-4.4 -o c_pnfw pnfw_contours.f root_finding.f lensing_functions.f arc_cross_subroutine.f") 
print "the compilation is done"

print "runing run_find_CC "

run('nfwpot',ks,rs,0.0,0.0,e_L=el,theta_L=thetal,show_plot=1,write_to_file=1)

print "runing the fortran code"
os.system("./c_pnfw")
print " files fort.65 and fort.75 correspond to tangential curves (SOLID LINES)"
print " files fort.25 and fort.35 correspond to radial curves (SOLID LINES)"

os.system("xmgrace -view 0.15 0.15 0.85 0.85 fort.65 lens_curves_tan.dat fort.25 lens_curves_rad.dat -batch edit_critical -saveall lens_plane.agr")
os.system("xmgrace -view 0.15 0.15 0.85 0.85 fort.75 -block lens_curves_tan.dat -bxy 3:4 fort.35 -block lens_curves_rad.dat -bxy 3:4 -batch edit_caustic -saveall src_plane.agr")

os.system("rm fort.*1 fort.*2 fort.*3 fort.*4")

