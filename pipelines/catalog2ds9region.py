#!/usr/bin/env python
# ====================================================
# Authors:
# Angelo Fausti Neto - angelofausti@gmail.com
# Cristina Furlanetto - furlanetto.cristina@gmail.com
# ====================================================


"""Module to create SAO DS9 region files from FITS catalogs"""

##@package catalog2ds9region
#
# This package creates SAO DS9 region files for objects in FITS catalogs.
# It has four mandatory arguments: the (input) FITS catalog, the reference center 
# coordinates (RA and DEC) and the search radius. It searches the objects in the 
# catalog that are inside a circle with this radius and centered at this reference 
# position.
# For each object it uses the catalog information to create an ellipse. It uses A_IMAGE 
# and B_IMAGE from the catalog as the ellipse semi-major axis and semi-minor axis, 
# respectively. It also uses THETA_IMAGE from the catalog as the ellipse orientation.
# The ellipse semi-major and semi-minor axes can be multiplied by a scale
# factor to enlarge the size of the object region. 
#
# Executable package: YES
#
# To see the package help message, type:
#
# > python catalog2ds9region.py --help
#
# To run the code, just type:
#
# > python catalog2ds9region.py -i input_file.fits --ra 10 --dec 20

import sys;
from numpy import *;
from math import cos, radians;
from catalog.fits_data import *;


def catalog2ds9region(input_file, ra_0, dec_0, radius, pixel_scale, scale_factor, output_file):
	''' 
	Main function to create SAO DS9 region files for objects in FITS catalogs. 
	It searches objects in the catalog that are inside a circle with a  
	given radius and centered at a reference position (RA and DEC). 
	For each object it uses the catalog information to create an ellipse. The ellipse 
	semi-major and semi-minor axes can be multiplied by a scale factor to enlarge the 
	size of the object region.

	Input:
	- input_file <str>: Input FITS catalog
	- ra_0 <float>: RA value of reference (deg)
	- dec_0 <float>: DEC value of reference (deg)
	- radius <float>: Search objetcs inside the radius (arcsec) 
	- pixel_scale <float>: Angular size of the pixel
	- scale_factor <float>: Scale the ellipse size
	- output_file <str>: Output SAO DS9 region file
	
	Output:
	- 0 success or 1 fail
	'''

	ra_0 = float(ra_0)

	if ra_0 < 0:	
		ra_0 = (360. + ra_0)

        dec_0 = float(dec_0)
        
	radius = float(radius) 				
	pixel_scale = float(pixel_scale) 
	scale_factor = float(scale_factor)

	tbhdu = pyfits.open(input_file)[1]
	
	data = dict_from_tbHDU(tbhdu, "RA", "DEC", "A_IMAGE", "B_IMAGE", "THETA_IMAGE")
    	
	search_radius = ((data['RA'] - ra_0)**2)*(math.cos(math.radians(dec_0))**2) + (data['DEC'] - dec_0)**2

	select = search_radius < (radius/3600.)**2

	index = where(select)[0]
	
	if len(index) == 0:
		print "No objects found."

	output = open(output_file,"w")
	
	# DS9 region file header
	output.write("# Region file format: DS9 version 4.1 \n")
	output.write("# Filename: %s \n" %input_file) 
	output.write("global color=green dashlist=8 3 width=1 font=\"helvetica 10 normal\" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n")
	output.write("fk5\n")

	for i in index:
		output.write("ellipse(%f,%f,%f\",%f\",%f) \n" %(data['RA'][i],data['DEC'][i], pixel_scale*scale_factor*data['A_IMAGE'][i],pixel_scale*scale_factor*data['B_IMAGE'][i], data['THETA_IMAGE'][i]))


#\cond	
if __name__ == "__main__" :

	from optparse import OptionParser;

	parser = OptionParser();

	parser.add_option('-i','--input_file',
			  dest='input_file', default=None,
			  help='Input FITS catalog');
	parser.add_option('-o','--output_file',
			  dest='output_file', default='ds9.reg',
			  help='Output ds9 region file [ds9.reg]');
	parser.add_option('--ra',
			  dest='ra_0', default=None,
			  help='Right ascension value of reference (deg)');
	parser.add_option('--dec',
			  dest='dec_0', default=None,
			  help='Declination value of reference (deg)');
	parser.add_option('-r','--radius',
			  dest='radius', default=30,
			  help='Search objetcs inside the radius (arcsec) [30]');
	parser.add_option('-p','--pixel_scale',
			  dest='pixel_scale', default=0.27,
			  help='Angular size of the pixel [0.27]');
	parser.add_option('-s','--scale_fator',
			  dest='scale_factor', default=2.5,
			  help='Scale the ellipse size [2.5]');

	(opts,args) = parser.parse_args();

	input_file = opts.input_file;
	output_file = opts.output_file;
	ra_0 = opts.ra_0;
	dec_0 = opts.dec_0;
	radius = opts.radius;
	pixel_scale = opts.pixel_scale;
	scale_factor = opts.scale_factor;

	if ( opts=={} ):
		print"";
		parser.print_help();
		print "";
		sys.exit(1);
	if ( input_file==None or ra_0 == None or dec_0 == None ):
		print "";
		parser.print_help();
		print "";
		sys.exit(1);

	catalog2ds9region(input_file,ra_0,dec_0,radius,pixel_scale,scale_factor,output_file)

	print "Created %s file." % output_file

	sys.exit(0);

#\endcond 


