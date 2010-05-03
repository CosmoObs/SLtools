# -*- coding: utf-8 -*-
#!/usr/bin/env python
import sys;

import os
import pyfits
import string
import commands
import numpy as np;


#===================================================
# Recorta pedaco da imagem englobando o cluster
def imcopy( input_file, output_file='output.fits', RAo=None, DECo=None, Xo=None, Yo=None, x_size=1000, y_size=1000 ):


	if ( RAo==None and Xo==None ) or ( DECo==None and Yo==None ):
		print >> sys.stderr, "Choose to input RA,DEC *or* X,Y coordinates";
		return (False);
	elif ( Xo == None or Yo == None ):
		_out = commands.getoutput( 'sky2xy %s %s %s' % (input_file,RAo,DECo) );
		_out = string.split(_out);
		x_halo = int(_out[-2]);
		y_halo = int(_out[-1]);
	else:
		x_halo = int(Xo);
		y_halo = int(Yo);


	imagem, hdr = pyfits.getdata( input_file,header=True );

	CRPIX1 = hdr['CRPIX1'];
	CRPIX2 = hdr['CRPIX2'];

	half_x = int(x_size)/2;
	half_y = int(y_size)/2;

	size_y, size_x = imagem.shape;

	y_ini = y_halo - int(y_size)/2;
	x_ini = x_halo - int(x_size)/2;

	LTV1 = None;
	LTV2 = None;

	if (x_ini >= 0) & (y_ini >= 0) :
		imagemnova = imagem[ (y_ini-1):(y_ini+int(y_size)), (x_ini-1):(x_ini+int(x_size)) ];
		LTV1 = -1*x_ini;
		LTV2 = -1*y_ini;
		CRPIX1 = CRPIX1 + LTV1;
		CRPIX2 = CRPIX2 + LTV2;
	elif (x_ini < 0) :
		imagemnova = imagem[ (y_ini-1):(y_ini+int(y_size)), 0:(x_ini+int(x_size)) ];
		LTV2 = -1*y_ini;
		CRPIX2 = CRPIX2 + LTV2;
	elif (y_ini < 0) :
		imagemnova = imagem[ 0:(y_ini+int(y_size)), (x_ini-1):(x_ini+int(x_size)) ]; 
		LTV1 = -1*x_ini;
		CRPIX1 = CRPIX1 + LTV1;

	if (x_ini < 0) & (y_ini < 0) :
		imagemnova = imagem[ 0:(y_ini+int(y_size)), 0:(x_ini+int(x_size)) ]; 

# Remover CD1_2 e CD2_1 do header
# Adicionar as seguintes chaves
	WCSDIM = 2
	CDELT1  =  -7.5000002980000E-5
	CDELT2  =  7.50000029800001E-5
	LTM1_1 = 1.0
	LTM2_2 = 1.0
	WAT0_001= 'system=image'
	WAT1_001= 'wtype=tan axtype=ra'
	WAT2_001= 'wtype=tan axtype=dec'

# Atualizacao do Header para a imagem nova
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

# Escreve a imagem nova
	os.system('[ -f %s ] && rm %s' % (output_file,output_file));
	pyfits.writeto( output_file, imagemnova, hdr );

	return;

# --

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
	
	if ( infits==None ):
		parser.print_help();
		sys.exit();
	if ( ra==None and x==None ) or ( dec==None and y==None ):
		parser.print_help();
		sys.exit();

	imcopy( input_file=infits, output_file=outfits, RAo=ra, DECo=dec, Xo=x, Yo=y, x_size=dx, y_size=dy );

# -
