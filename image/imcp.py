"""@package imcp

Module to deal with (fits) image copy/cuts

Funtions:
   cutout( FITS image file ) -> (np.array, header)

"""

# -*- coding: utf-8 -*-
#!/usr/bin/env python
import sys;

import os
import pyfits
import string
import commands
import numpy as np;

#============================================================================================================================
def cutout( input_file, RAo=None, DECo=None, Xo=None, Yo=None, x_cut_size=1000, y_cut_size=1000, output_file=None ):

	"""-> (np.array,header)

	Do a (square) cut out of given FITS image

	Input:
	 input_file: Input FITS image file
	 RAo, DECo: Centre of output (cut) image in sky coordinates
	 Xo, Yo: Centre of output (cut) image in pixel position
	 x_cut_size: Horizontal size (in pixels) of output image
	 y_cut_size: Vertical size (in pixels) of output image

	Output:
	 output_file: Output FITS (cut) image file. If *None*, the np.array and header are returned.

	"""

	hdr = pyfits.getheader( input_file );

	NAXIS1 = hdr['NAXIS1'];
	NAXIS2 = hdr['NAXIS2'];

	CRPIX1 = hdr['CRPIX1'];
	CRPIX2 = hdr['CRPIX2'];

	LTV1 = None;
	LTV2 = None;

	if (RAo != None and DECo != None) and (Xo == None and Yo == None):
		_out = commands.getoutput( 'sky2xy %s %s %s' % (input_file,RAo,DECo) );
		_out = string.split(_out);
		x_halo = int(_out[-2])-1;
		y_halo = int(_out[-1])-1;
	elif (Xo != None and Yo != None) and (RAo == None and DECo == None):
		x_halo = int(Xo)-1;
		y_halo = int(Yo)-1;
	elif (Xo != None and Yo != None and RAo != None and DECo != None):
		x_halo = int(Xo)-1;
		y_halo = int(Yo)-1;
		print >> sys.stderr, "Warning: You should be passing just x,y or RA,DEC values.";
		print >> sys.stderr, "         Using array (x,y) coordinates as the cut-off center";
	else:
		x_halo = int(NAXIS1/2);   x_cut_size = int(NAXIS1/2);
		y_halo = int(NAXIS2/2);   y_cut_size = int(NAXIS2/2);
		print >> sys.stderr, "Warning: No central coordinates were given for snapshot.";
		print >> sys.stderr, "         Using image central pixels as the cut-off center.";

	imagem = pyfits.getdata( input_file );

	x_img_size = NAXIS1;
	y_img_size = NAXIS2;

	y_ini = y_halo - int(y_cut_size/2);
	x_ini = x_halo - int(x_cut_size/2);
	y_fin = y_ini + y_cut_size;
	x_fin = x_ini + x_cut_size;

	if (y_ini < 0):
		y_ini = 0;
	if (x_ini < 0):
		x_ini = 0;
	if (y_fin >= y_img_size):
		y_fin = y_img_size - 1;
	if (x_fin >= x_img_size):
		x_fin = x_img_size - 1;

	y_rest = y_img_size%2;
	x_rest = x_img_size%2;

	imagemnova = imagem[ (y_ini):(y_fin), (x_ini):(x_fin) ]; 

	if (y_fin-y_ini != y_cut_size) or (x_fin-x_ini != x_cut_size):
		print >> sys.stderr, "Warning: Image side sizes are different than given ones.";
		print >> sys.stderr, "         Possible cause: given central position near the edge.";

	if (x_ini >= 0) & (y_ini >= 0) :
		LTV1 = -1*x_ini;
		LTV2 = -1*y_ini;
		CRPIX1 = CRPIX1 + LTV1;
		CRPIX2 = CRPIX2 + LTV2;
	elif (x_ini < 0) :
		LTV2 = -1*y_ini;
		CRPIX2 = CRPIX2 + LTV2;
	elif (y_ini < 0) :
		LTV1 = -1*x_ini;
		CRPIX1 = CRPIX1 + LTV1;

	# Remover CD1_2 e CD2_1 do header
	# Adicionar as seguintes chaves
	#
	WCSDIM = 2
	CDELT1  =  -7.5000002980000E-5
	CDELT2  =  7.50000029800001E-5
	LTM1_1 = 1.0
	LTM2_2 = 1.0
	WAT0_001= 'system=image'
	WAT1_001= 'wtype=tan axtype=ra'
	WAT2_001= 'wtype=tan axtype=dec'

	# Atualizacao do Header para a imagem nova
	#
	if (LTV1 != None) :
		hdr.update('LTV1',LTV1)
		hdr.update('CRPIX1',CRPIX1)
	if (LTV2 != None) :
		hdr.update('LTV2',LTV2)
		hdr.update('CRPIX2',CRPIX2)
	hdr.update('WCSDIM',WCSDIM)
	hdr.update('CDELT1',CDELT1)
	hdr.update('CDELT2',CDELT2)
	hdr.update('LTM1_1',LTM1_1)
	hdr.update('LTM2_2',LTM2_2)
	hdr.update('WAT0_001',WAT0_001)
	hdr.update('WAT1_001',WAT1_001)
	hdr.update('WAT2_001',WAT2_001)

	if (output_file == None):
		return (imagemnova,hdr);

	# Escreve a imagem nova se um nome de arquivo foi dado
	#
	os.system('[ -f %s ] && rm %s' % (output_file,output_file));
	pyfits.writeto( output_file, imagemnova, hdr );

	return (True);
# --------------------

# =========================
if __name__ == "__main__" :

	from optparse import OptionParser;

	parser = OptionParser();

	parser.add_option('-f','--input_fits',
			  dest='infits', default=None,
			  help='Input FITS image');
	parser.add_option('-o','--output_fits',
			  dest='outfits', default='output.fits',
			  help='Output FITS image [output.fits]');
	parser.add_option('--ra',
			  dest='ra', default=None,
			  help='Right Ascension');
	parser.add_option('--dec',
			  dest='dec', default=None,
			  help='Declination');
	parser.add_option('--xo',
			  dest='x', default=None,
			  help='X centre');
	parser.add_option('--yo',
			  dest='y', default=None,
			  help='Y centre');
	parser.add_option('--x_size',
			  dest='dx', default=1000,
			  help='Size of x side on output');
	parser.add_option('--y_size',
			  dest='dy', default=1000,
			  help='Size of y side on output');

	(opts,args) = parser.parse_args();

	infits = opts.infits;
	outfits = opts.outfits;
	ra = opts.ra;
	dec = opts.dec;
	x = opts.x;
	y = opts.y;
	dx = opts.dx;
	dy = opts.dy;

	if ( opts=={} ):
		print"";
		parser.print_help();
		print "";
		sys.exit(2);
	if ( infits==None ):
		print "";
		parser.print_help();
		print "";
		sys.exit(2);

	cutout( input_file=infits, output_file=outfits, RAo=ra, DECo=dec, Xo=x, Yo=y, x_cut_size=dx, y_cut_size=dy );

	sys.exit(0);
# ------------------

