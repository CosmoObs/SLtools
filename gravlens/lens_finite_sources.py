from __future__ import division
import os
import logging


# This function generates sersic sources centered at (source_centers) and runs SExtractor on them. Inputs: source parameters (source_centers, eS, rS, nS). Output: fits files with the corresponding images, and the SE outputs (objects fits images and segmentation fits images, the segmentation being ran twice, for 2 different deblending values). Output: a matrix containing the properties of the source for each source


##@package lens_finite_sources
# lenses sources (centered in source_centers) and measures the resulting images, identifing which ones are arcs
#
#@param lens_model 
#@param zs
#@param delta_zs
#@param sextractor_params
#@param arc_selection_params
#@return list with the name of the arcs
def lens_finite_sources(inputlens, a, dimpix, source_centers, ref_magzpt, reference_band, source_model):
	Npix = int( (2*a) / dimpix)
	image_names = []
	f199 = open('sersicinput.txt', 'w')
	f199.write(inputlens)
	for i in range (0,len(source_centers)):
		source_type = source_model[i]['source_type'] # type of the source (can be 'sersic' or 'uniform')

		if source_type == 'sersic':
			f199.write('setsource 1\n')
			# determine 'totalsorceplaneflux', according to the source magnitude
			totalsourceplaneflux = 10**((2/5)*(ref_magzpt - source_model[i][reference_band]))
			#------------------------------------------------------------------------
			f199.write('%s %f %f %f %f %f %f 0 %f macro\n' % (source_type, totalsourceplaneflux, source_centers[i][0], source_centers[i][1], source_model[i]['es'], source_model[i]['thetas'], source_model[i]['rs'], source_model[i]['ns']  ) ) # sersic/uniform F x y e PA halflightr nothing nS macro/micro
			f199.write('0 0 0 0 0 0 0 0\n')
			f199.write('SBmap2 %0.2f %0.2f %d %0.2f %0.2f %d 1 sbmap%05d_%s.fits 3\n' % (-a,a,Npix,-a,a,Npix,i, reference_band) ) # <x lo> <hi> <# steps> <y lo> <hi> <# steps> <Nover> <file> <outtype>	
			image_names.append('sbmap%05d_%s.fits' % (i, reference_band) )

		if source_type == 'uniform':
			f199.write('setsource 1\n')
			# fixed 'totalsorceplaneflux' to 1
			totalsourceplaneflux = 1
			# Nover increase resolution if the source is smaller than the pixel
			nover = source_model[i]['nover']
			while source_model[i]['rs'] < dimpix/nover and nover < 3:
				nover += 1
			source_model[i]['nover'] = nover
			print "nover = ", nover
			#------------------------------------------------------------------------
			f199.write('%s %f %f %f %f %f %f 0 0 macro\n' % (source_type, totalsourceplaneflux, source_centers[i][0], source_centers[i][1], source_model[i]['es'], source_model[i]['thetas'], source_model[i]['rs'] ) ) # sersic/uniform F x y e PA halflightr nothing nS macro/micro
			f199.write('0 0 0 0 0 0 0 0\n')
			f199.write('SBmap2 %0.2f %0.2f %d %0.2f %0.2f %d %d sbmap%05d_%s.fits 3\n' % (-a,a,Npix,-a,a,Npix, nover, i, reference_band) ) # <x lo> <hi> <# steps> <y lo> <hi> <# steps> <Nover> <file> <outtype>	
			image_names.append('sbmap%05d_%s.fits' % (i, reference_band) )

	f199.close()
	print "gravlens is lensing %d finite source(s) (this may take several minutes depending on the resolution)..." % len(source_centers)
	if len(source_centers) > 0:
		status = os.system('gravlens sersicinput.txt > /dev/null')
		logging.debug('Executed gravlens to lens the finite sources and returned status %s' % status)
	return image_names








