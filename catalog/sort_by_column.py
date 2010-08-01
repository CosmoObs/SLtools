""" sort_by_column """

##@package sort_by_column

import os
import numpy

def sort_by_column(catalog, column):

    """ 
    Sort the ascii catalog by the i-th column in increasing order 
    
    Inputs:

    catalog : name of the catalog file
    column : number of the column to sort

    Output:

    The sorted catalog

    """

    root_name, extension = os.path.splitext(catalog)

    n, x, y, mag, theta, a, elip, cs, fl = numpy.loadtxt(catalog, unpack=True) 
    
    z = zip(n,x,y,mag,theta,a,elip,cs,fl)
  
    z.sort(key=lambda x:x[column])

    a = zip(*z)

    b = numpy.transpose(a)
	
    sorted_catalog = root_name + "_sorted." + extension

    numpy.savetxt(sorted_catalog, b,fmt="%12.6G")

    return sorted_catalog		
