#!/usr/bin/python
#
import string
import numpy
import os
import astropy.io.fits as pyfits
import math
import pywcs
import sys



##################################################################################################
def radec2xy(ra,dec,wcs):

	skycrd = numpy.array([[ra,dec]])

	pixcrd = wcs.wcs_sky2pix(skycrd,1)

	x = pixcrd[0][0]

	y = pixcrd[0][1]

	return x,y

##################################################################################################
def modify_star_catalog(sdss_star_catalog,image_name):

#Read and modify the sdss star catalog in order to include the x and y columns 
# Convert sky coordinates (ra,dec) to world coordinates


	sdss_catalog = open(sdss_star_catalog, "r")
	
	input_root, input_extension = os.path.splitext(sdss_star_catalog)

	img_root, img_extension = os.path.splitext(image_name)
	
	band = img_root[-1]

	output_file_name = input_root + "_new2_" + band + input_extension

	output_file = open(output_file_name, "w")

	hdulist = pyfits.open(image_name) # Load the FITS hdulist using pyfits

	wcs = pywcs.WCS(hdulist[0].header) # Parse the WCS keywords in the primary HDU


	for line in sdss_catalog:
		if "#" in line:
			continue
		obj_id = string.split(line,sep=None)[0]
		ra = float(string.split(line,sep=None)[1])
		dec = float(string.split(line,sep=None)[2])
		print ra, dec
		g_mag = string.split(line,sep=None)[3]
		r_mag = string.split(line,sep=None)[4]
		i_mag = string.split(line,sep=None)[5]

		x,y = radec2xy(ra,dec,wcs)

		output_file.write("%s %s %s %s %s %s %s %s \n" %(obj_id,ra,dec,x,y,g_mag,r_mag,i_mag))
	
	return output_file_name

###################################################################################################
def run_sextractor(input_image_name,phot_apert,mag_zpt):
	# Run SExtractor on an image.
	

	hdulist = pyfits.open(input_image_name)
	seeing = 1.2# hdulist[0].header["DIMMSEE"]

	#input_root = string.split(string.split(input_image_name,sep='/')[-1],sep='.')[0] #Separate the image name from the input path and from the extension
	input_root, input_extension = os.path.splitext(input_image_name)



	# Create the SExtractor default configuration file and the parameter file.
	#os.system('sex -dd > default.sex')
	os.system('echo "X_IMAGE" > default.param')
	os.system('echo "Y_IMAGE" >> default.param')
	os.system('echo "MAG_AUTO" >> default.param')
	os.system('echo "MAGERR_AUTO" >> default.param')
	os.system('echo "MAG_BEST" >> default.param')
	os.system('echo "MAGERR_BEST" >> default.param')
	os.system('echo "MAG_APER" >> default.param')
	os.system('echo "MAGERR_APER" >> default.param')
	os.system('echo "MAG_ISO" >> default.param')
	os.system('echo "MAGERR_ISO" >> default.param')
	os.system('echo "MAG_ISOCOR" >> default.param')
	os.system('echo "MAGERR_ISOCOR" >> default.param')
	os.system('echo "CLASS_STAR" >> default.param')

	catalog = input_root + "_aux.cat"
	DEBLEND_MINCONT = 0.005
	DEBLEND_NTHRESH = 20
	DETECT_MINAREA = 10
	THRESH_TYPE = "RELATIVE"
	DETECT_THRESH = 2.0
	ANALYSIS_THRESH = 2.0
	PIXEL_SCALE = 0.154
	

	#Run SExtractor
	os.system('sex %s -c input.sex -CATALOG_NAME %s -SEEING_FWHM %f -PHOT_APERTURES %f -DEBLEND_MINCONT %f -DEBLEND_NTHRESH %d -DETECT_MINAREA %d -THRESH_TYPE %s -DETECT_THRESH %f -ANALYSIS_THRESH %f -PIXEL_SCALE %f -MAG_ZEROPOINT %f -CHECKIMAGE_TYPE SEGMENTATION -CHECKIMAGE_NAME %s' %(input_image_name, catalog, seeing, phot_apert, DEBLEND_MINCONT, DEBLEND_NTHRESH, DETECT_MINAREA, THRESH_TYPE, DETECT_THRESH, ANALYSIS_THRESH, PIXEL_SCALE, mag_zpt,input_root + "_seg.fits")) 
	
	sex_catalog = input_root + ".cat"
	os.system('grep -v "^#" %s > %s' %(catalog, sex_catalog))
	os.system('rm %s' %catalog)

	return sex_catalog

###################################################################################################
def run_sextractor_final_dual(input_image_name,phot_apert,mag_zpt,band,ref_image): 
	# Run SExtractor on an image in dual image mode.
	

	hdulist = pyfits.open(input_image_name)
	seeing = 1.2 #hdulist[0].header["DIMMSEE"]

	#input_root = string.split(string.split(input_image_name,sep='/')[-1],sep='.')[0] #Separate the image name from the input path and from the extension
	input_root, input_extension = os.path.splitext(input_image_name)


	# Create the SExtractor default configuration file and the parameter file.
	#os.system('sex -dd > default.sex')
	os.system('echo "NUMBER" > default.param')
	os.system('echo "X_IMAGE" >> default.param')
	os.system('echo "Y_IMAGE" >> default.param')
	os.system('echo "ALPHA_J2000" >> default.param')
	os.system('echo "DELTA_J2000" >> default.param')
	os.system('echo "MAG_AUTO" >> default.param')
	os.system('echo "MAGERR_AUTO" >> default.param')
	os.system('echo "A_IMAGE" >> default.param')
	os.system('echo "B_IMAGE" >> default.param')
	os.system('echo "ELLIPTICITY" >> default.param')
	os.system('echo "THETA_IMAGE" >> default.param')
	os.system('echo "CLASS_STAR" >> default.param')
	os.system('echo "FLAGS" >> default.param')


	sex_catalog_final = input_root + "_final.cat"
	DEBLEND_MINCONT = 0.005
	DEBLEND_NTHRESH = 20
	DETECT_MINAREA = 10
	THRESH_TYPE = "RELATIVE"
	DETECT_THRESH = 2.0
	ANALYSIS_THRESH = 2.0
	PIXEL_SCALE = 0.154
	

	#Run SExtractor
	os.system('sex %s %s -c input.sex -CATALOG_NAME %s -SEEING_FWHM %f -PHOT_APERTURES %f -DEBLEND_MINCONT %f -DEBLEND_NTHRESH %d -DETECT_MINAREA %d -THRESH_TYPE %s -DETECT_THRESH %f -ANALYSIS_THRESH %f -PIXEL_SCALE %f -MAG_ZEROPOINT %f' %(ref_image,input_image_name, sex_catalog_final, seeing,phot_apert, DEBLEND_MINCONT, DEBLEND_NTHRESH, DETECT_MINAREA, THRESH_TYPE, DETECT_THRESH, ANALYSIS_THRESH, PIXEL_SCALE, mag_zpt)) 
	


	return sex_catalog_final


##########################################################################################################
def run_sextractor_final(input_image_name,phot_apert,mag_zpt,band): 
	# Run SExtractor on an image
	

	hdulist = pyfits.open(input_image_name)
	seeing = 1.2 #hdulist[0].header["DIMMSEE"]

	#input_root = string.split(string.split(input_image_name,sep='/')[-1],sep='.')[0] #Separate the image name from the input path and from the extension
	input_root, input_extension = os.path.splitext(input_image_name)


	# Create the SExtractor default configuration file and the parameter file.
	#os.system('sex -dd > default.sex')
	os.system('echo "NUMBER" > default.param')
	os.system('echo "X_IMAGE" >> default.param')
	os.system('echo "Y_IMAGE" >> default.param')
	os.system('echo "ALPHA_J2000" >> default.param')
	os.system('echo "DELTA_J2000" >> default.param')
	os.system('echo "MAG_AUTO" >> default.param')
	os.system('echo "MAGERR_AUTO" >> default.param')
	os.system('echo "A_IMAGE" >> default.param')
	os.system('echo "B_IMAGE" >> default.param')
	os.system('echo "ELLIPTICITY" >> default.param')
	os.system('echo "THETA_IMAGE" >> default.param')
	os.system('echo "CLASS_STAR" >> default.param')
	os.system('echo "FLAGS" >> default.param')


	sex_catalog_final = input_root + "_final.cat"
	DEBLEND_MINCONT = 0.005
	DEBLEND_NTHRESH = 20
	DETECT_MINAREA = 10
	THRESH_TYPE = "RELATIVE"
	DETECT_THRESH = 2.0
	ANALYSIS_THRESH = 2.0
	PIXEL_SCALE = 0.154
	

	#Run SExtractor
	os.system('sex %s -c input.sex -CATALOG_NAME %s -SEEING_FWHM %f -PHOT_APERTURES %f -DEBLEND_MINCONT %f -DEBLEND_NTHRESH %d -DETECT_MINAREA %d -THRESH_TYPE %s -DETECT_THRESH %f -ANALYSIS_THRESH %f -PIXEL_SCALE %f -MAG_ZEROPOINT %f' %(input_image_name, sex_catalog_final, seeing,phot_apert, DEBLEND_MINCONT, DEBLEND_NTHRESH, DETECT_MINAREA, THRESH_TYPE, DETECT_THRESH, ANALYSIS_THRESH, PIXEL_SCALE, mag_zpt)) 
	


	return sex_catalog_final


##########################################################################################################
def mean_stddev(delta_mag,sig_mag2,delta_mag2):

	# Mean
	
	mean_mag = numpy.sum(delta_mag) / numpy.sum(sig_mag2)		
	
	# Standard deviation
	stddev_mag_sum = 0.

	
	for k in range(len(delta_mag)):
			
			stddev_mag_sum = stddev_mag_sum + ((delta_mag2[k] - mean_mag)**2 ) * sig_mag2[k]


	stddev_mag = stddev_mag_sum / (len(delta_mag) - 1) #((numpy.sum(sig_mag2)) * (len(delta_mag) - 1))

 
	return mean_mag,stddev_mag


###########################################################################################################
def match_stars(star_catalog,sex_catalog,band,mag_zpt,phot_aper,opt,a,b):

	# Matches the stars in the Sextractor catalog with the SDSS stars in the same field by their position (x,y).
	input_root, input_extension = os.path.splitext(sex_catalog)
#	result_file = open(input_root + "_result_" + str(phot_aper) + ".dat","w")
	result_file_final = open(input_root + "_final_result.dat","w")

	sex_obj = open(sex_catalog, "r").readlines()
	sdss_stars = open(star_catalog, "r").readlines()


	index_star_sex = []
	index_star_sdss = []
	sdss_g = []
	sdss_r = []
	sdss_i = []
	

	delta_mag_auto = []
	delta_mag_best = []
	delta_mag_aper = []
	delta_mag_iso = []
	delta_mag_isocor = []
	delta_mag_auto2 = []
	delta_mag_best2 = []
	delta_mag_aper2 = []
	delta_mag_iso2 = []
	delta_mag_isocor2 = []
	sig_mag_auto = []
	sig_mag_best = []
	sig_mag_aper = []
	sig_mag_iso = []
	sig_mag_isocor = []
	sig_mag_auto2 = []
	sig_mag_best2 = []
	sig_mag_aper2 = []
	sig_mag_iso2 = []
	sig_mag_isocor2 = []
	mag_auto = []
	mag_best = []
	mag_aper = []
	mag_iso = []
	mag_isocor = []
	magerr_auto = []
	magerr_best = []
	magerr_aper = []
	magerr_iso = []
	magerr_isocor = []

	mag_type = ["auto","best","aper","iso","isocor"]

	mag = [mag_auto,mag_best,mag_aper,mag_iso,mag_isocor]
	err = [magerr_auto,magerr_best,magerr_aper,magerr_iso,magerr_isocor]
	delta_mag = [delta_mag_auto,delta_mag_best,delta_mag_aper,delta_mag_iso,delta_mag_isocor]
	sig_mag = [sig_mag_auto,sig_mag_best,sig_mag_aper,sig_mag_iso,sig_mag_isocor]
	sig_mag2 = [sig_mag_auto2,sig_mag_best2,sig_mag_aper2,sig_mag_iso2,sig_mag_isocor2]
	delta_mag2 = [delta_mag_auto2,delta_mag_best2,delta_mag_aper2,delta_mag_iso2,delta_mag_isocor2]
	mean_mag = ["mean_mag_auto","mean_mag_best","mean_mag_aper","mean_mag_iso","mean_mag_isocor"]
	stddev_mag = ["stddev_mag_auto","stddev_mag_best","stddev_mag_aper","stddev_mag_iso","stddev_mag_isocor"]

	for i in range( len(sdss_stars)):
		sdss_x = float(sdss_stars[i].split()[3])
		sdss_y = float(sdss_stars[i].split()[4])
		sdss_g.append(float(sdss_stars[i].split()[5]))
		sdss_r.append(float(sdss_stars[i].split()[6]))
		sdss_i.append(float(sdss_stars[i].split()[7]))

		if band == "i":
			sdss_mag = float(sdss_stars[i].split()[7])
		elif band == "g":
			sdss_mag = float(sdss_stars[i].split()[5])
		elif band == "r":		
			sdss_mag = float(sdss_stars[i].split()[6])
		else:
			print "I can't understand the band"
			break

		for j in range(len(sex_obj)):

			sex_x = float(sex_obj[j].split()[0])
			sex_y = float(sex_obj[j].split()[1])
			sex_class_star = float(sex_obj[j].split()[12])

			if  30. <= sex_x <= 2198. and 30. <= sex_y <= 2019. and sex_class_star > 0.85: ##Relacionar com limite superior da imagem

				
				mag_auto.append(float(sex_obj[j].split()[2]))
				magerr_auto.append(float(sex_obj[j].split()[3]))
				mag_best.append(float(sex_obj[j].split()[4]))
				magerr_best.append(float(sex_obj[j].split()[5]))
				mag_aper.append(float(sex_obj[j].split()[6]))
				magerr_aper.append(float(sex_obj[j].split()[7]))
				mag_iso.append(float(sex_obj[j].split()[8]))
				magerr_iso.append(float(sex_obj[j].split()[9]))
				mag_isocor.append(float(sex_obj[j].split()[10]))
				magerr_isocor.append(float(sex_obj[j].split()[11]))

			
				
				if (sdss_x - 5.0) <= sex_x <= (sdss_x + 5.0) and (sdss_y - 5.0) <= sex_y <= (sdss_y + 5.0):
					index_star_sex.append(j)
					index_star_sdss.append(i)

					if band == "i":
					
						if  15.5 <= sdss_mag <= 18.5:

							print sex_x,sex_y,mag[0][-1],sdss_x,sdss_y,sdss_mag
							for k in range(5):
					
								delta_mag[k].append((mag[k][-1] - sdss_mag) )# / (err[k][-1])**2)
	
								sig_mag2[k].append(1.0 ) #/ (err[k][-1])**2)
	
								delta_mag2[k].append(mag[k][-1] - sdss_mag)

					else:
						if  15.5 <= sdss_mag <= 19.:

							print sex_x,sex_y,mag[0][-1],sdss_x,sdss_y,sdss_mag

							for k in range(5):
						
					
								delta_mag[k].append((mag[k][-1] - sdss_mag) )# / (err[k][-1])**2)
	
								sig_mag2[k].append(1.0 ) #/ (err[k][-1])**2)
		
								delta_mag2[k].append(mag[k][-1] - sdss_mag)

	print index_star_sdss


	for i in range(5):
		print i
		flag = 1

		#clipping

		while flag == 1:
	
			mean_mag[i],stddev_mag[i] = mean_stddev(delta_mag[i],sig_mag2[i],delta_mag2[i])
			print mean_mag[i],math.sqrt(stddev_mag[i])
	
			lim1 = 2.5*math.sqrt(stddev_mag[i]) #for negative residuals
			lim2 = 2.5*math.sqrt(stddev_mag[i]) #for positive residuals			
		
			diff = delta_mag2[i] - mean_mag[i]

			index = []

			for kk in range(len(diff)):

				if diff[kk] >= 0 and diff[kk] <= lim2 or diff[kk] < 0 and diff[kk] > -lim1 :
					index.append(kk)
				
					

			#index = numpy.where(diff<lim)[0]
			print index,len(index)
			
			if 0 < len(index) and len(index) < len(delta_mag2[i]):
	
					
				delta_mag_new = []
				sig_mag2_new = []
				delta_mag2_new = []
	
				for k in index:
					
					delta_mag_new.append(delta_mag[i][k])
					sig_mag2_new.append(sig_mag2[i][k])
					delta_mag2_new.append(delta_mag2[i][k])
			
				delta_mag[i] = delta_mag_new
				sig_mag2[i] = sig_mag2_new
				delta_mag2[i] = delta_mag2_new
			else:	
				flag = 0
		
	
	

	#index = stddev_mag.index(numpy.min(stddev_mag)) #Takes the index of the lowest standard deviation
	index = 0  #takes mag_auto
	
	final_mean = mean_mag[index] #Takes the mean magnitude value corresponding to the lowest standard deviation
	final_type = mag_type[index] #Takes the final type magnitude


	#Final magnitude zero-point
	final_mag_zpt = mag_zpt - final_mean
	

	#Preparing plot files

	#color_plot = open(input_root + "_plot.dat","w")
	#color_plot2 = open(input_root + "_plot2.dat","w")
	
	soma = 0.
	soma2 = 0.		

	for i in range(len(index_star_sdss)):
		#SDSS magnitudes for the matched stars
		if band == "i":
			mag_star_sdss = sdss_i[index_star_sdss[i]]
		elif band == "g":
			mag_star_sdss = sdss_g[index_star_sdss[i]]
		elif band == "r":		
			mag_star_sdss = sdss_r[index_star_sdss[i]]
		else:
			print "I can't understand the band"
			break
	

		g_r = sdss_g[index_star_sdss[i]] - sdss_r[index_star_sdss[i]]
		r_i = sdss_r[index_star_sdss[i]] - sdss_i[index_star_sdss[i]]
	
		mag_final = mag[index][index_star_sex[i]] - final_mean #Sextractor magnitude for matched stars corrected by the zero-point magnitude
	
		delta_mag_sdss_sex = mag_final - mag_star_sdss 

		#color_plot.write("%s %s %s \n" %(delta_mag_sdss_sex,g_r,r_i))

		if opt == 1:
		
			mag_final2 = mag_final - a - b * g_r

			delta_mag_sdss_sex2 = mag_final2 - mag_star_sdss

			soma2 = soma2 + delta_mag_sdss_sex2**2
			soma = soma + delta_mag_sdss_sex2
		
			#color_plot2.write("%s %s %s \n" %(delta_mag_sdss_sex2,g_r,r_i))


	sigma = math.sqrt(soma2 / (len(index_star_sdss) - 1))
	result_file_final.write("Star catalog: %s \n" %star_catalog)
	result_file_final.write("SExtractor catalog: %s \n" %sex_catalog)
	result_file_final.write("Phot aperture: %s \n" %phot_aper)
	result_file_final.write("Mean (mag_auto, mag_best, mag_aper, mag_iso, mag_isocor)\n")
	result_file_final.write("%s %s %s %s %s \n"  %(mean_mag[0],mean_mag[1],mean_mag[2],mean_mag[3],mean_mag[4]))
	result_file_final.write("Standard deviation (mag_auto, mag_best, mag_aper, mag_iso, mag_isocor)\n")
	result_file_final.write("%s %s %s %s %s \n"  %(stddev_mag[0],stddev_mag[1],stddev_mag[2],stddev_mag[3],stddev_mag[4]))
	result_file_final.write("Chosen magnitude: %s (%s)\n" %(final_mean,final_type))	

	return final_mag_zpt


##########################################################################################################
def combine_catalogs_dual(sex_catalogs,final_catalog_name):

	#Combine catalogs that were obtained using SExtractor with dual image mode. 

	number = []
	x_image = []
	y_image = []
	alpha_j200 = []
	delta_j200 = []
	mag_auto = []
	magerr_auto = []
	a_image = []
	b_image = []
	ellipticity = []
	theta_image = []
	class_star = []
	flags = []

	for catalog_name in sex_catalogs:
		print catalog_name
		params = numpy.loadtxt(catalog_name, unpack=True)
		number.append(params[0])
		x_image.append(params[1])
		y_image.append(params[2])
		alpha_j200.append(params[3])
		delta_j200.append(params[4])
		mag_auto.append(params[5])
		magerr_auto.append(params[6])
		a_image.append(params[7])
		b_image.append(params[8])
		ellipticity.append(params[9])
		theta_image.append(params[10])
		class_star.append(params[11])
		flags.append(params[12])
		
	#0-> g; 1-> r; 2-> i.

	z = zip(number[0],x_image[0],y_image[0],alpha_j200[0],delta_j200[0],mag_auto[0],magerr_auto[0],mag_auto[1],magerr_auto[1],mag_auto[2],magerr_auto[2],a_image[0],b_image[0],ellipticity[0],theta_image[0],class_star[0],flags[0],flags[1],flags[2])

	final_catalog = open(final_catalog_name,"w")
	final_catalog.write("#   1 NUMBER          Running object number\n")
	final_catalog.write("#   2 X_IMAGE         Object position along x                         [pixel]\n")
	final_catalog.write("#   3 Y_IMAGE         Object position along y                         [pixel]\n")
	final_catalog.write("#   4 ALPHA_J2000     Right ascension of barycenter (J2000)           [deg]\n")
	final_catalog.write("#   5 DELTA_J2000     Declination of barycenter (J2000)               [deg]\n")
	final_catalog.write("#   6 MAG_AUTO_G      Kron-like elliptical aperture magnitude         [mag]\n")
	final_catalog.write("#   7 MAGERR_AUTO_G   RMS error for AUTO magnitude                    [mag]\n")
	final_catalog.write("#   8 MAG_AUTO_R      Kron-like elliptical aperture magnitude         [mag]\n")
	final_catalog.write("#   9 MAGERR_AUTO_R   RMS error for AUTO magnitude                    [mag]\n")
	final_catalog.write("#  10 MAG_AUTO_I      Kron-like elliptical aperture magnitude         [mag]\n")
	final_catalog.write("#  11 MAGERR_AUTO_I   RMS error for AUTO magnitude                    [mag]\n")
	final_catalog.write("#  12 A_IMAGE         Profile RMS along major axis                    [pixel]\n")
	final_catalog.write("#  13 B_IMAGE         Profile RMS along minor axis                    [pixel]\n")
	final_catalog.write("#  14 ELLIPTICITY     1 - B_IMAGE/A_IMAGE\n")
	final_catalog.write("#  15 THETA_IMAGE     Position angle (CCW/x)                          [deg]\n")
	final_catalog.write("#  16 CLASS_STAR      S/G classifier output\n")
	final_catalog.write("#  17 FLAGS_G          Extraction flags\n")
	final_catalog.write("#  18 FLAGS_R          Extraction flags\n")
	final_catalog.write("#  19 FLAGS_I          Extraction flags\n")

	numpy.savetxt(final_catalog, z, fmt="%5i %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %8.3f %8.3f %8.3f %8.2f %6.2f %5i %5i %5i" )


##############################################################################################

def match_catalogs_i_r(r_catalog,i_catalog):

	#Matches i and r band catalogs

	i_cat = numpy.loadtxt(i_catalog)
	r_cat = numpy.loadtxt(r_catalog)

	i_r_catalog = "tmp_ir.cat"
	
	i_r_cat = open(i_r_catalog, "w")
	

	ra_i = i_cat[:,3]
	dec_i = i_cat[:,4]
	ra_r = r_cat[:,3]
	dec_r = r_cat[:,4]

	i_index = []
	r_index = []
	i_match_index = []
	r_match_index = []

	for l in range(len(ra_r)):
		r_index.append(l)


	for k in range(len(ra_i)):
		i_index.append(k)
		
		for j in range(len(ra_r)):


			if ((ra_i[k]- 0.00014) <= ra_r[j] <= (ra_i[k] + 0.00014)) and ((dec_i[k] - 0.00014) <= dec_r[j] <= (dec_i[k] + 0.00014)):
				i_match_index.append(k)
				r_match_index.append(j)
				flag_cris = 0
							

				i_r_cat.write("%d %s  %s  %s  %s  %s  %s  %s  %s  %s  %s  %s  %s  %s  %d %d \n" %(i_cat[k,0],i_cat[k,1],i_cat[k,2],i_cat[k,3],i_cat[k,4],i_cat[k,5],i_cat[k,6],r_cat[j,5],r_cat[j,6],i_cat[k,7],i_cat[k,8],i_cat[k,9],i_cat[k,10],i_cat[k,11],i_cat[k,12],flag_cris))
				

			
	dif_i = numpy.setdiff1d(i_index,i_match_index)

	#Inserting in the i-r catalog the objects in i band that don't matched with objects in r band
	flag_cris = 0
	for j in range(len(dif_i)):
	
		k = dif_i[j]
		i_r_cat.write("%d  %s  %s  %s  %s  %s  %s  99.00  99.00  %s  %s  %s  %s  %s %d %d \n" %(i_cat[k,0],i_cat[k,1],i_cat[k,2],i_cat[k,3],i_cat[k,4],i_cat[k,5],i_cat[k,6],i_cat[k,7],i_cat[k,8],i_cat[k,9],i_cat[k,10],i_cat[k,11],i_cat[k,12],flag_cris))
				


	return i_r_catalog
				
#----------------------------------------------------------------------------------------------------------------------------------------------------	
def combine_catalogs(sex_catalogs,final_catalog):

	
	i_r_catalog = match_catalogs_i_r(sex_catalogs[1],sex_catalogs[2])
	i_r_cat = numpy.loadtxt(i_r_catalog)

	
	g_cat = numpy.loadtxt(sex_catalogs[0])
	
	i_r_g_cat = open(final_catalog, "w")


	ra_i = i_r_cat[:,3]
	dec_i = i_r_cat[:,4]
#	ra_r = i_r_cat[:,6]
#	dec_r = i_r_cat[:,7]
	ra_g = g_cat[:,3]
	dec_g = g_cat[:,4]

	i_r_index = []
	g_index = []
	i_r_match_index = []
	g_match_index = []


	i_r_g_cat.write("#   1 NUMBER          Running object number\n")
	i_r_g_cat.write("#   2 X_IMAGE         Object position along x                         [pixel]\n")
	i_r_g_cat.write("#   3 Y_IMAGE         Object position along y                         [pixel]\n")
	i_r_g_cat.write("#   4 ALPHA_J2000     Right ascension of barycenter (J2000)           [deg]\n")
	i_r_g_cat.write("#   5 DELTA_J2000     Declination of barycenter (J2000)               [deg]\n")
	i_r_g_cat.write("#   6 MAG_AUTO_G      Kron-like elliptical aperture magnitude         [mag]\n")
	i_r_g_cat.write("#   7 MAGERR_AUTO_G   RMS error for AUTO magnitude                    [mag]\n")
	i_r_g_cat.write("#   8 MAG_AUTO_R      Kron-like elliptical aperture magnitude         [mag]\n")
	i_r_g_cat.write("#   9 MAGERR_AUTO_R   RMS error for AUTO magnitude                    [mag]\n")
	i_r_g_cat.write("#  10 MAG_AUTO_I      Kron-like elliptical aperture magnitude         [mag]\n")
	i_r_g_cat.write("#  11 MAGERR_AUTO_I   RMS error for AUTO magnitude                    [mag]\n")
	i_r_g_cat.write("#  12 A_IMAGE         Profile RMS along major axis                    [pixel]\n")
	i_r_g_cat.write("#  13 B_IMAGE         Profile RMS along minor axis                    [pixel]\n")
	i_r_g_cat.write("#  14 ELLIPTICITY     1 - B_IMAGE/A_IMAGE\n")
	i_r_g_cat.write("#  15 THETA_IMAGE     Position angle (CCW/x)                          [deg]\n")
	i_r_g_cat.write("#  16 CLASS_STAR      S/G classifier output\n")
	i_r_g_cat.write("#  17 FLAGS_G          Extraction flags\n")



	for l in range(len(ra_g)):
		g_index.append(l)
		


 	for k in range(len(ra_i)):
 		i_r_index.append(k)

 		
		for j in range(len(ra_g)):
 		

 			if ((ra_g[j]- 0.00014) <= ra_i[k] <= (ra_g[j] + 0.00014)) and ((dec_g[j] - 0.00014) <= dec_i[k] <= (dec_g[j] + 0.00014)):
				
 				i_r_match_index.append(k)
 				g_match_index.append(j)
 				
							

 				i_r_g_cat.write("%5i %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %8.3f %8.3f %8.3f %8.2f %6.2f %5i \n" %(i_r_cat[k,0],i_r_cat[k,1],i_r_cat[k,2],i_r_cat[k,3],i_r_cat[k,4],g_cat[j,5],g_cat[j,6],i_r_cat[k,7],i_r_cat[k,8],i_r_cat[k,5],i_r_cat[k,6],i_r_cat[k,9],i_r_cat[k,10],i_r_cat[k,11],i_r_cat[k,12],i_r_cat[k,13],i_r_cat[k,14]))
				


 			
 	dif_i_r = numpy.setdiff1d(i_r_index,i_r_match_index)

	flag_mag = 99.	
	
	flag_cris = 0
	for j in range(len(dif_i_r)):
	
		k = dif_i_r[j]
		
		i_r_g_cat.write("%5i %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %8.3f %8.3f %8.3f %8.2f %6.2f %5i \n" %(i_r_cat[k,0],i_r_cat[k,1],i_r_cat[k,2],i_r_cat[k,3],i_r_cat[k,4],flag_mag,flag_mag,i_r_cat[k,7],i_r_cat[k,8],i_r_cat[k,5],i_r_cat[k,6],i_r_cat[k,9],i_r_cat[k,10],i_r_cat[k,11],i_r_cat[k,12],i_r_cat[k,13],i_r_cat[k,14]))



##############################################################################################
#Read input file with the list of SOGRAS images.

input_file = open("input.txt","r")

for line in input_file.readlines():
	if "#" in line:
		continue	

	input_image,star_catalog = string.split(line,sep=None)
	
img_root,img_extension = os.path.splitext(input_image)
print img_root

ref_band = "i"	
ref_image = img_root[0:-1] + ref_band + "_register"+  img_extension
print ref_image

bands = ["g","r","i"]
	
sex_catalogs = []

for band in bands:
	print band

	image_name = img_root[0:-1] + band + img_extension
	print image_name

	register_image =  img_root[0:-1] + band + "_register" + img_extension

	print register_image

	#A partir das coordenadas RA,dEC dos catalogo obtenhas as coordenadas x,y para nossas imagens do SOAR.

	star_catalog_new = modify_star_catalog(star_catalog,image_name)
	
	#Photometric calibration: determine the zero point magnitude				

	mag_zpt_ini = 20.0

	phot_aper = 30

	sex_catalog = run_sextractor(image_name,phot_aper,mag_zpt_ini)	

	final_mag_zpt = match_stars(star_catalog_new,sex_catalog,band,mag_zpt_ini,phot_aper,0,0,0)
	os.system ('rm %s' %sex_catalog)

	
	sex_catalog_final = run_sextractor_final(register_image,phot_aper,final_mag_zpt,band)

	sex_catalogs.append(sex_catalog_final)

#Combine catalogs	

final_catalog = img_root[0:-1] + "final_phot.cat"


#combine_catalogs_dual(sex_catalogs,final_catalog)

combine_catalogs(sex_catalogs,final_catalog)

	
