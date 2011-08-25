#!/usr/bin/python
# ===================================================
# Authors:
# Cristina Furlanetto - furlanetto.cristina@gmail.com
# ===================================================
#
# 
'''Module to fit arcs parameters using ArcEllipse analytical expression and Sersic profile.'''

##@module arcfitting
#
#
# This module uses the analytical expression of ArcEllipse combined with a Sersic profile to fit the
# arc parameters from an arc image. The ArcEllipse expression is the result of a deformation of an 
# ellipse such that one of its main axes becomes a circle segment. The Sersic profile is used to 
# describe the radial light distribution of the arc. 
#
# It uses chi2 minimization to find the best fit parameters values. The chi2 minimization is done 
# using Pyminuit (see <A HREF="http://code.google.com/p/pyminuit/">).
#
#
#
from __future__ import division
import os
from numpy import *
import math
import pyfits
import sltools
import string
import minuit
from sltools.image.header_funcs import get_header_parameter
from stools.geometry.get_extrema_pts import *
from sltools.geometry.elementary_geometry import *
from tools.get_bisectrix_pts import *
from tools.circle3points import *


#========================================================================
#Ad-hoc parameters:
sig_frac = 0.95 #Fraction of the total arc signal used for segmentation.
frac_n = 1.E-2
#========================================================================


def f(x,y,x0,y0,rc,ellip,I0,re,n,theta0):

	'''
	Function to compute the value of the surface brightness of the
	(x,y) pixel following the ArcEllipse prescription.

	Input:
	- x: x position (pixel)
	- y: y position (pixel)
	- x0: x position of the point where the arc is centered (pixel)
	- y0: y position of the point where the arc is centered (pixel)
	- r_c: distance of the arc center related to the point where the arc is centered (pixel)
	- ellip: ellipticity of the arc
	- I_0: central surface brightness of the Sersic profile
	- r_e: radius that encloses half of the total luminosity of the Sersic profile (pixel)
	- n: index that describes the slope of the Sersic profile 
	- theta0: orientation of the arc center related to the x-axis (radians)

	Output:
	- the value of the surface brightness of the (x,y) point following the ArcEllipse prescription
	
	'''
	
	x = float(x)
	y = float(y)	

	theta = math.atan2((y - y0),(x - x0))
		
	r = math.sqrt((x - x0)**2 + (y - y0)**2)
	theta_diff = abs(theta - theta0)
	if (theta_diff >= (math.pi)):
		theta_diff = 2. * math.pi - theta_diff

	r_prof = math.sqrt((rc * theta_diff * (1. - ellip))**2 + (r - rc)**2 )
        bn = 1.9992*float(n) - 0.3271
	exponent = (1. / float(n))
	
	exp_arg = - bn * (r_prof/re)**exponent
	
	return I0 * math.exp(exp_arg)	



def segmentation(img_name):

	'''
	Function to segmentate the arc image. 
	
	The arc image segmentation is done by taking all the pixels 
	whose values (counts) are greater than a minimum value. The minimum pixel value or cutoff, Cmin, 
	is the one for which the accumulated signal (sum of the pixel values in decreasing order) is 
	equal to a fraction of the total arc signal.

	Input:
	- img_name: arc image name

	Output:
	- (ximg): list of x pixels of the segmented arc
	- (yimg): list of y pixels of the segmented arc
	- (seg_index): list of the pixel index of the segmented arc in the original image
	- xmax: x position corresponding to the maximum arc intensity
	- ymax: y position corresponding to the maximum arc intensity
	
	'''

	#Get image data
	data = pyfits.getdata(img_name) 
	#Get image parameters
	dim_x,dim_y = get_header_parameter(img_name,'NAXIS1','NAXIS2')

	#Find the maximum arc intensity and the total signal
	vmax = -1
	signal = []
	for i in range(0,dim_x): 
		for j in range(0,dim_y):
			signal.append(data[j][i])
			
			if data[j][i] > vmax:
				vmax = data[j][i]
				xmax =i
				ymax = j

	signal_total = sum(signal)
	
	#Segmentation
	signal_sort = sorted(signal,reverse=True)

	signal_cumul = 0.

	for k in range(len(signal_sort)):
	
		signal_cumul += signal_sort[k]

		if signal_cumul >= sig_frac * signal_total:	
			signal_cut = signal_sort[k]
			break
	
	seg_index = where(data >= signal_cut) # Index of the segmentation image pixels

	ximg = seg_index[1]
	yimg = seg_index[0]

	return (ximg),(yimg),seg_index,xmax,ymax



def initial_guesses(img_name,ximg,yimg,seg_index,xmax,ymax):
	
	'''
	Function to find initial guesses to be used as input for pyminuit for the following 
	parameters:
	- x0: x position of the point where the arc is centered (pixel)
	- y0: y position of the point where the arc is centered (pixel)
	- r_c: distance of the arc center related to the point where the arc is centered (pixel)
	- ellip: ellipticity of the arc
	- I_0: central surface brightness of the Sersic profile
	- r_e: radius that encloses half of the total luminosity of the Sersic profile (pixel)
	- n: index that describes the slope of the Sersic profile 
	- theta0: orientation of the arc center related to the x-axis (radians)

	The initial guesses are obtained using arc measurements.
	First, it finds the extrema points of the arc, called p1 and p2.
	After, it finds the points near the bisectrix of the line joining p1 and p2. 
	The central point of the arc, p3, is obtained as the mean of the bisectrix points weighted by 
	the intensity of the points. The choise of p3 as a point on the bisectrix instead of the point
	of the arc with the maximum intensity was done to guarantee that it is the geometrical centrer
	of the arc. 
	The lenght of the arc is obtained as the sum of the line segments that connect p1 and p3 and 
	p3 and p2.
	The width of the arc was obtained as the variance (weighted by the intensity) of the bisectrix
	points along the bisectrix line. 
	The ellipticity is obtained using the expression: e = 1 - 1/(L/W).
	The x0, y0 and r_c are obtained by finding a circle that fits 3 points of the arc (p1,p2 and p3).
	The thetha_0 is obtained using the expression: thetha_0 = arctan ((yc - y0)/(xc - x0)), 
	where (xc,yc) are the coordinates of the p3 point. 
	The I_0 is obtained as the intensity at p3. 


	Input:
	- img_name: arc image name
	- ximg: list of x pixels of the segmented arc
	- yimg: list of y pixels of the segmented arc
	- seg_index: list of the pixel index of the segmented arc in the original image
	- xmax: x position corresponding to the maximum arc intensity
	- ymax: y position corresponding to the maximum arc intensity

	Output:
	- (x0,y0,rc,ellip,I0,re,theta0): list with the initial guesses values of the parameters 
	
	'''
	
	image_name, image_extension = os.path.splitext(img_name)

	output = open(image_name + "_init_guesses.txt","w")
	
	#Get image data
	data = pyfits.getdata(img_name) 

	#Get extrema points p1 and p2 (to compute length)
	l1,l2=get_extrema_2loops(ximg,yimg,20)

	p1 = (seg_index[1][l1],seg_index[0][l1])
	p2 = (seg_index[1][l2],seg_index[0][l2])
	#output.write("p1 = %s \np2 = %s" %(p1,p2))

	#Get the bisectrix points
	xb,yb,costheta,sentheta,midx,midy = get_bisectrix_pts(ximg, yimg,l1,l2)
	output.write("x_bisectrix = %s \ny_bisectrix= %s\n" %(xb,yb))

	#Get width using the intensity distribution
	#Get the mean point in the bisectrix points
	xb_m = []
	yb_m = []
	intens_m = []
	for k in range(len(xb)):
	
		xb_m.append(data[yb[k]][xb[k]]*xb[k])
		yb_m.append(data[yb[k]][xb[k]]*yb[k])
		intens_m.append(data[yb[k]][xb[k]])

	xb_mean = sum(xb_m)/sum(intens_m)
	yb_mean = sum(yb_m)/sum(intens_m)
	output.write("xbisect_mean = %f, ybisect_mean = %f\n" %(xb_mean,yb_mean))

	#Compute sigma
	yl = []
	for i in range (0,len(xb)):
		yl.append(-sentheta*(xb[i] - midx) + costheta*(yb[i] - midy))

	yl_m = []
	intens_l_m = []
	for k in range(len(yl)):
	
		yl_m.append(data[yb[k]][xb[k]]*yl[k])
		intens_l_m.append(data[yb[k]][xb[k]])

	yl_mean = sum(yl_m)/sum(intens_l_m)

	s_sigmay_l = 0.
	s_intens_l = 0.

	for kk in range(len(yl)):

		s_sigmay_l = s_sigmay_l + (data[yb[kk]][xb[kk]]*yl[kk] -data[yb_mean][xb_mean]*yl_mean)**2
		s_intens_l = s_intens_l + (data[yb[kk]][xb[kk]])**2
		
	sigmay_l = math.sqrt(s_sigmay_l/s_intens_l/(len(yl)-1))
	output.write("yl_mean = %f\n" %(yl_mean))
	output.write("sigmay_l = %f\n" %(sigmay_l))

	#Get point with maximum intensity as the center of the arc
	p3_max = (xmax,ymax)

	#Get lenght using the mean point in the bisectrix points
	p3 = (xb_mean,yb_mean)
	l_1_3 = math.sqrt((p3[0] - p1[0])**2 + (p3[1] - p1[1])**2)
	l_2_3 = math.sqrt((p3[0] - p2[0])**2 + (p3[1] - p2[1])**2)
	length = l_1_3 + l_2_3
	output.write("p1 = %s, p2= %s, p3=%s\n" %(p1,p2,p3))
	
	#Get rc,x0 and y0 by finding a circle that fits the 3 points (p1,p2 and p3)
	x_0,y_0,r_c = three_points_to_circle(p1, p2, p3)

	#Get arc position angle
	theta_0 = math.atan2((yb_mean-y_0),(xb_mean-x_0)) 
	output.write("r_c = %f \nx_0 = %f \ny_0 = %f \ntheta_0 = %f\n" %(r_c,x_0,y_0,theta_0))

	#Get width using arc area
	area = len(ximg)
	width = area / ( math.pi * length)
	output.write("length = %f \narea = %f \nwidth = %f\n" %(length,area,width))

	#Get L/W 
	l_w = length/sigmay_l

	ellip_0 = 1. - 1. / l_w

	I_0 = data[ymax][xmax]

	output.write("l_w = %f  \nellip_0 = %f\n" %(l_w,ellip_0))

	return (x_0,y_0,r_c,ellip_0,I_0,sigmay_l,theta_0)



def chi2(x0,y0,rc,ellip,I0,re,n,theta0):

	'''
	Function to compute the reduced chi2 of the f(x,y,x0,y0,rc,ellip,I0,re,n,theta0) function related 
	to the data. Since we do not have an error associated with the data, we use Poissonic error
	(sigma[i][j] = sqrt(data[i][j]) ).

	Input:
	- x0: x position of the point where the arc is centered (pixel)
	- y0: y position of the point where the arc is centered (pixel)
	- r_c: distance of the arc center related to the point where the arc is centered (pixel)
	- ellip: ellipticity of the arc
	- I_0: central surface brightness of the Sersic profile
	- r_e: radius that encloses half of the total luminosity of the Sersic profile (pixel)
	- n: index that describes the shape of the Sersic profile 
	- theta0: orientation of the arc center related to the x-axis (radians)

	Output:
	- the value of the reduced chi2
	
	'''

	c2 = 0.

	for k in range(len(ximg)):
		i = ximg[k]
		j = yimg[k]

		c2 += ((f(i,j,x0,y0,rc,ellip,I0,re,n,theta0)-data[j][i])**2)/(data[j][i])

	return c2/float(len(ximg) - 1)


	

def find_n_scan(params,n_ini,n_fin,n_n):

	'''
	Function to find the best fit value for n by "brute force", considering other parameter as fixed. 
	It uses the function scan from Pyminuit, which calculates chi2 function values along a lattice 
	in parameter space. 

	Input:
	- params: list with the initial guesses values of the parameters (x0,y0,rc,ellip,I0,re,theta0)
	- n_ini: initial value for n for scan evaluation
	- n_fin: final value for n for scan evaluation
	- n_n: the number of subdivisions to evaluate 

	Output:
	- n_scan: best fit value for n
	'''
	

	print "find n scan"
	m = minuit.Minuit(chi2)

	m.values["x0"] = params[0]
	m.values["y0"] = params[1]
	m.values["rc"] = params[2]
	m.values["ellip"] = params[3]
	m.values["I0"] = params[4]
	m.values["re"] = params[5]
	m.values["theta0"] = params[6]


	m.scan(("n",n_n,n_ini,n_fin))

	n_scan = m.values["n"]
	print n_scan
	return n_scan





def find_n_migrad(params,n_scan):

	'''
	Function to find the best fit value for n. It uses the function migrad from pyminuit.
	It uses best fit value for n obtained with the scan method (n_scan) as initial guess for n.

	Input:
	- params: list with the initial guesses values of the parameters (x0,y0,rc,ellip,I0,re,theta0)
	- n_scan: initial guess value for n

	Output:
	- n_migrad:best fit value for n
	
	'''


	print "find n migrad"
	m = minuit.Minuit(chi2,fix_x0 = True, fix_y0 = True,fix_rc = True, fix_ellip = True,fix_I0 = True, fix_re = True,fix_theta0 = True)

	m.values["x0"] = params[0]
	m.values["y0"] = params[1]
	m.values["rc"] = params[2]
	m.values["ellip"] = params[3]
	m.values["I0"] = params[4]
	m.values["re"] = params[5]
	m.values["theta0"] = params[6]
	m.values["n"] = n_scan
	#m.printMode = 1

	m.limits["n"] = ((0.1,10))
	m.migrad()

	n_migrad = m.values["n"]

	return n_migrad




def find_best_fit(params,nf):

	'''
	Function to find the best fit values for all the parameters, except for n. 
 	The best fit values are found by minimizing the chi2 function, using the function migrad from pyminuit and
	taking n as a fixed parameter.

	Input:
	- params: list with the initial guesses values of the parameters (x0,y0,rc,ellip,I0,re,theta0)
	- n_f: initial guess value for n

	Output:
	- params_new: list with the best fit values of the parameters (x0,y0,rc,ellip,I0,re,theta0)
	- chi2_min: reduced chi2 minimum value
	
	'''

	print "find best fit params"
	m = minuit.Minuit(chi2,fix_n = True)

	m.values["x0"] = params[0]
	m.values["y0"] = params[1]
	m.values["rc"] = params[2]
	m.values["ellip"] = params[3]
	m.values["I0"] = params[4]
	m.values["re"] = params[5]
	m.values["n"] = nf
	m.values["theta0"] = params[6]

	m.migrad() 

	x0 = m.values["x0"]
	y0 = m.values["y0"]
	rc = m.values["rc"]
	ellip = m.values["ellip"]
	I0 = m.values["I0"]
	re = m.values["re"]
	theta0 = m.values["theta0"]

	params_new = (x0,y0,rc,ellip,I0,re,theta0)

	chi2_min = m.fval


	return params_new,chi2_min




def find_best_fit_all(params,nf):

	'''
	Function to find the best fit values for all the parameters. 
 	The best fit values are found by minimizing the chi2 function, using the function migrad from pyminuit.

	Input:
	- params: list with the initial guesses values of the parameters (x0,y0,rc,ellip,I0,re,theta0)
	- n_f: initial guess value for n

	Output:
	- params_new: list with the best fit values of the parameters (x0,y0,rc,ellip,I0,re,theta0)
	- chi2_min: reduced chi2 minimum value
	
	'''

	print "find best fit params all"
	m = minuit.Minuit(chi2)

	m.values["x0"] = params[0]
	m.values["y0"] = params[1]
	m.values["rc"] = params[2]
	m.values["ellip"] = params[3]
	m.values["I0"] = params[4]
	m.values["re"] = params[5]
	m.values["n"] = nf
	m.values["theta0"] = params[6]

	m.migrad() 

	x0 = m.values["x0"]
	y0 = m.values["y0"]
	rc = m.values["rc"]
	ellip = m.values["ellip"]
	I0 = m.values["I0"]
	re = m.values["re"]
	theta0 = m.values["theta0"]
	n = m.values["n"]

	params_new = (x0,y0,rc,ellip,I0,re,theta0,n)
	chi2_min = m.fval

	return params_new,chi2_min

#==========================================================
#Running ArcFitting
#
# It reads the input file containing a list of arc images to be fitted. 
#
# It finds the initial guesses for all the arc parameters, except for n. The initial guess for n is 
# obtained using the function find_n_scan.
#
# The minimization of the reduced chi2 with respect to the parameter n is done separately, but paralell, 
# from the minimization with respect to the other 7 parameters. 
#
# It uses Migrad algorithm, keeping n fixed at its initial values, to obtain the best fit values for 
# the other parameters. With the new values of the parameters, it recalculates the value
# of the parameter n (with find_n_scan) and compares this value with its previous one. If the difference 
# between these values is greater than a given threshold, the process continues and it recalculates 
#the values of the 7 parameters until it finds a convergence in the value of n.
#


images_list = "input_images.txt"
file_in = open(images_list,'r')

for line in file_in.readlines():
	if "#" in line:
		continue
	
	image_name = string.split(line)[0]
	print "\n"
	print image_name

	image_root, image_extension = os.path.splitext(image_name)
	
	output_file = open(image_root + "_fitted_params.txt","w")
	output_file.write("#Image name  x_0(pix)   y_0(pix)   r_c(pix)   ellip   I_0   r_e(pix)   theta_0(rad)   n   chi_2_min \n")

	ximg,yimg,seg_index,xmax,ymax = segmentation(image_name)

	params = initial_guesses(image_name,ximg,yimg,seg_index,xmax,ymax)

	#Using Pyminuit

	data = pyfits.getdata(image_name) #Get image data

	#Find initial n
	n_ini = 0.4
	n_fin = 6.
	n_n = 50

	nf = []
	nf.append(0.)

	nf.append(find_n_scan(params,n_ini,n_fin,n_n))

	params_new,chi2_min = find_best_fit(params,nf[1])

	for i in range(2,1000):
		print i

		n_ini = nf[i-1] - 0.5 * nf[i-1]
	
		if n_ini < 0.4:
			n_ini = 0.4

		n_fin = nf[i-1] + 0.5 * nf[i-1]
		
		nf.append(find_n_scan(params_new,n_ini,n_fin,n_n))
		
		params_new,chi2_min = find_best_fit(params_new,nf[i])
	
		if abs(nf[i] - nf[i-1]) > frac_n and abs(nf[i] - nf[i-2])  > (frac_n / 10.):

			continue
	
		else:
			
			break
	

	print params_new,nf[-1]

	print image_name, params_new[3],nf[-1],params_new[5],params_new[2],params_new[6],params_new[0],params_new[1],params_new[4]
	print "chi2 min: %f" %chi2_min
	output_file.write("%s %s %s %s %s %s %s %s %s %s" %(image_name,params_new[0],params_new[1],params_new[2],params_new[3],params_new[4],params_new[5],params_new[6],nf[-1],chi2_min))


