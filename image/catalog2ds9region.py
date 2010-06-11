#!/usr/bin/python
#
from numpy import *
import pyfits
import os

def get_parameters_from_catalog(cat_name,ra_0,dec_0,radius_arcsec,output_file):

	#cat_root_name, cat_extension = os.path.splitext(cat_name)
	
	hdulist = pyfits.open(cat_name,ignore_missing_end = True)
	tbdata = hdulist[1].data

	ra = tbdata.field("RA")
	dec = tbdata.field("DEC")
	a = tbdata.field("A_IMAGE")
	b = tbdata.field("B_IMAGE")
	tilename = tbdata.field("TILENAME")
	
	r = sqrt(((ra - ra_0)**2)*(math.cos(math.radians(dec_0)))**2 + (dec - dec_0)**2)
	
	select = r < radius_arcsec / 3600.

	index = where(select)[0]

	pix_scale = 0.27

	output = open(output_file,"w")
	
	output.write("# Region file format: DS9 version 4.1 \n")
	output.write("# Filename: DES2223.fits \n")
	output.write("global color=green dashlist=8 3 width=1 font=\"helvetica 10 normal\" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n")
	output.write("fk5\n")


	for i in index:
		output.write("ellipse(%f,%f,%f\",%f\",0.) \n" %(ra[i],dec[i], 2.5*a[i]*pix_scale,2.5*b[i]*pix_scale))
	


get_parameters_from_catalog("DC5B_FOOTPRINT_3.fits",334.92784,-42.40495,50.,"tmp.reg") 

