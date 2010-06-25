#!/usr/bin/env python
#-*- coding:utf-8 -*-

"""@package catalog2ds9region

Module to create SAO DS9 region files from DES Data Challenge FITS catalogs

Executable package: YES

To see the package help message, type:

> python catalog2ds9region.py --help

 'catalog2ds9region' has three mandatory arguments: the (input) FITS catalog and
 the reference center coordinates (RA and DEC)

> python catalog2ds9region.py -i input_file.fits --ra 10 --dec 20

"""

from numpy import *
import math
import sys

# TODO: remove halos/get_catalog_data, fix select_halos.py and select_halos_by_radec to use the current module

# move get_catalog_data to tools/catalog and use  Intra-package References (http://docs.python.org/tutorial/modules.html)
# the solution below didn't work

#sys.path.append(os.path.abspath(os.path.abspath("../../"))
#from tools.catalog.get_catalog_data import get_catalog_data

from get_catalog_data import get_catalog_data

def catalog2ds9region(input_file, ra_0, dec_0, radius, pixel_scale, scale_factor, output_file):
''' Main function to create DS9 region file 

	@param input_file : Input FITS catalog
	@param ra_0, dec_0 : RA and DEC values of reference (deg)
	@param radius : Search objetcs inside the radius (arcsec) 
	@param pixel_scale : Angular size of the pixel
	@param scale_factor : Scale the ellipse size
	@param output_file : Output ds9 region file

	@return 0 success or 1 fail
'''

	ra_0 = float(ra_0)

	if ra_0 < 0:	
		ra_0 = (360. + ra_0)

        dec_0 = float(dec_0)
	radius = float(radius) 				
	pixel_scale = float(pixel_scale) 
	scale_factor = float(scale_factor)

	data = get_catalog_data(input_file, "RA", "DEC", "A_IMAGE", "B_IMAGE", "THETA_IMAGE", "TILENAME")
    	
	search_radius = ((data['RA'] - ra_0)**2)*(math.cos(math.radians(dec_0))**2) + (data['DEC'] - dec_0)**2

	select = search_radius < (radius/3600.)**2

	index = where(select)[0]
	
	if len(index) == 0:
		print "No objects found."

	output = open(output_file,"w")
	
	# DS9 region file header
	output.write("# Region file format: DS9 version 4.1 \n")
	output.write("# Filename: DES2223.fits \n")
	output.write("global color=green dashlist=8 3 width=1 font=\"helvetica 10 normal\" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n")
	output.write("fk5\n")

        # Circle entities
	for i in index:
		output.write("ellipse(%f,%f,%f\",%f\",%f) \n" %(data['RA'][i],data['DEC'][i], pixel_scale*scale_factor*data['A_IMAGE'][i],pixel_scale*scale_factor*data['B_IMAGE'][i], data['THETA_IMAGE'][i]))


#@cond	
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

#@endcond 


