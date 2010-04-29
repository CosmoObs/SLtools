import os
import numpy

def ordering_catalog(cat_name, column):
	# Order the catalog by the values of the i-th column in increasing order.

	# Separate the root name and the extension of the catalog_name
	root_name, extension = os.path.splitext(cat_name)

	n, x, y, mag, theta, a, elip, cs, fl = numpy.loadtxt(cat_name, unpack=True) #Read the catalog
	
	z = zip(n,x,y,mag,theta,a,elip,cs,fl)

	z.sort(key=lambda x:x[column]) #Sort the catalog by the i-th column

	a = zip(*z)

	b = numpy.transpose(a)
	
	order_catalog = root_name + "_ord.cat"
#	write_array(order_catalog, b, precision = 3)		

	numpy.savetxt(order_catalog, b,fmt="%12.6G")  # Write the ordered catalog.
	
	return order_catalog		
