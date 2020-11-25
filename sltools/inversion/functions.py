# ==================================
# Author: Gabriel Bartosch Caminha - gbcaminha@gmail.com
# ==================================
"""
Package with stand alone functions of gbclib
"""
from __future__ import division

import astropy.io.fits as pyfits

import random

import subprocess

import numpy as np

from math import sqrt

import math

################################################################################
#Alone functions
################################################################################
def strip_list(str_list, str_strip):
    """
    Function to sripe characteres from a list os strings.

    This function expand the method str.strip.
    see http://docs.python.org/library/stdtypes.html?highlight=lstrip#str.strip

    Input:
     - str_list  str : list of strings
     - str_strip str : characteres to be striped
    Output:
     - string list without the characteres
    """
    for i in range(len(str_list)):
        str_list[i] = str_list[i].strip(str_strip)
    return str_list
################################################################################
def extract_parameter(file_in, identifier):
    """
    Function read the values of an ifentifier in a lenstool .par file

    It can be used to read values after the patern 'identifier' in a text file

    Input:
     - file_in         str : file name
     - identifier      str : identifier name
    Output:
     - array      str : identifier parameters
    """
    file_inf = open(file_in, 'r')
    out = []
    identifier = identifier.split()
    size_identifier = len(identifier)
    
    for line in file_inf:
        line_inf = line.split()

        if len(line_inf)>0 and \
            line_inf[:size_identifier] == identifier[:size_identifier]:
            out.append( line_inf[size_identifier:]  )
    file_inf.close()
    return out
################################################################################
def modulus(x_in, y_in):
    """
    Compute the modulus of a vector

    Input:
     - x_in float : x coordenate
     - y_in float : y coordenate
    Output:
     - float sqrt( x_in**2 + y_in**2 )
    """
    return sqrt( x_in**2 + y_in**2 )
################################################################################
def generate_sources(n_src, radius):
    """
    Generate sources within a radius

    Input:
     - n_src    int : number of sources
     - radius float : radius where the sources will be sorted
    Output:
     - [[x_image, y_image],...]
    """
    x_src = 2.0*radius
    y_src = 2.0*radius

    xy_list = []

    for _ in range(n_src):

        while x_src*x_src + y_src*y_src > radius*radius: 
            x_src = random.uniform(-radius, radius)
            y_src = random.uniform(-radius, radius)

        xy_list.append([x_src, y_src])

        # ensure that x*x + y*y > radius*radius true
        x_src = 2.0*radius 
        y_src = 2.0*radius
    return xy_list
################################################################################
def add_err_pos(image_xy, sigma_gauss, mu_gauss = 0.0):
    """
    Add poisson errors in a set of image postitions

    Input:
     - image_xy float [[,],...] : set of image positions
     - sigma              float : gaussian standard deviation
     - mu                 float : gaussian median
    Output:
     - float [[,],...]          : set of image positions with gaussian errors
    """
    image_xy_out = []
    x_add = 0.0
    y_add = 0.0
    for i in range(len(image_xy)):
        x_add = image_xy[i][0] + random.gauss(mu_gauss, sigma_gauss)
        y_add = image_xy[i][1] + random.gauss(mu_gauss, sigma_gauss)
        image_xy_out.append([x_add, y_add])
    return image_xy_out
################################################################################
def generate_substructure(z_lens, n_subs = 5, profil = 1, vdisp_limit = 100):
    """
    Generate substructures. For now, randomly sorts all properties. 

    Input:
     - z_lens float : lens red-shift
     - n_subs   int : number of substructures
     - profil   int : profile type (lenstool classification)
    Output:
     - [{},...]     : vector with dictonaries for each substructure
    """
    if profil != 1:
        print 'Warning: Profile ', profil, ' not implemented yet'
        print 'Warning: creating substructures using profile 1'
        profil = 1
    
    position_limit = 5

    subs = []
    
    for _ in range(n_subs):
        x_sub = random.uniform(-position_limit, position_limit)
        y_sub = random.uniform(-position_limit, position_limit)
        
        vdisp = random.uniform(0.0, vdisp_limit)
        
        ellipticite = random.uniform(0.0, 1.0)
        
        angle_pos = random.uniform(-90.0, 90.0)
        
        dic = {'z_lens': z_lens, 'x_centre': x_sub, 'y_centre': y_sub, \
               'profil': profil, 'angle_pos': angle_pos, \
               'ellipticite': ellipticite, 'v_disp': vdisp}

        subs.append(dic)
    return subs
################################################################################
def extract_second_identifiers(file_in, first_identifier):
    """
    Function to extract all second udentifiers under a firs identifier in a
    lenstool .par file

    Input:
     - file_in    str : file name
     - identifier str : identifier name
    Output:
     - array      str : identifier parameters
    """
    file_inf = open(file_in, 'r')
    out = []
    read_bool = False
    n_frst_id = -1

    for line in file_inf:
        line_inf = line.split()
        if len(line_inf)>0 and line_inf[0] == 'end':
            read_bool = False

        if read_bool:
            out[n_frst_id][line_inf[0]] = line_inf[1:]

        if len(line_inf)>0 and line_inf[0] == first_identifier:
            out.append({})
            read_bool = True
            n_frst_id = n_frst_id + 1

        #print out, '\n\n'
    #print out[1]
    file_inf.close()
    return out
################################################################################
def generate_fits_image(xy_matrx, file_name = 'out.fits'):
    """
    Function to generate a fits image from a matrix

    Input:
     - xy_matrx [[x1y1,x1y2...],[[x2y1,x2y2...]],... ] : file name
    Output:
     - NONE
    """
    hdu = pyfits.PrimaryHDU()
    hdu.data = xy_matrx

    hdu.writeto(file_name, clobber=True)
################################################################################
def add_slash(str_in):
    """
    Add a slash (/) in the end of str_in if it does not have

    Input:
     - str_in str : input string
    Output:
     - str_in + '/' (if str_in does not have a slash in the end)
    """
    if str_in[len(str_in)-1] != '/':
        return str_in + '/'
    return str_in
################################################################################
def is_X_running():
    """
    Function to check if X windows is running, adapted from:
    http://stackoverflow.com/questions/1027894/detect-if-x11-is-available-python

    Input:
     - NONE
    Output:
     - BOOLEAN : True if X window is running, False otherwise
    """
    from subprocess import Popen, PIPE
    p = Popen(["xset", "-q"], stdout=PIPE, stderr=PIPE)
    p.communicate()
    return p.returncode == 0
################################################################################
def randon_function(func, n_pts, x_min, x_max, *args):
    """
    Sort points in a given 'func' distribution

    Input:
     - func function : distribution where the points will be sorted
     - n_pts     int : number of points to be returned
     - x_min    
     - x_max
     - *args   *args : 'func' arguments
    Output:
     - [...]
    """
    x_vec = np.linspace(x_min, x_max, num = 100)
    y_vec = func(x_vec, *args)
    y_min = min(y_vec)
    y_max = max(y_vec)

    x_rand = []

    while len(x_rand) < n_pts:
        xr = np.random.random()
        yr = np.random.random()
        x = x_min + (x_max - x_min) * xr
        y = y_min + (y_max - y_min) * yr
        if y <= func(x, *args):
            x_rand.append(x)
    return x_rand
################################################################################
def rotate_vector(vec_in, theta):
    """
    Function to rotate a vector

    Input:
     - vec_in [] : vector with [x,y]
     - theta  float : angle to rotate the vector in degrees 
    Output:
     - rotated vector
    """
    a = 2.0*math.pi/360.0
    cos = math.cos(a*theta)
    sin = math.sin(a*theta)
    x_out = cos*vec_in[0] - sin*vec_in[1]
    y_out = sin*vec_in[0] + cos*vec_in[1]
    return x_out, y_out
################################################################################
def flatten(x):
    """flatten(sequence) -> list

    Returns a single, flat list which contains all elements retrieved
    from the sequence and all recursively contained sub-sequences
    (iterables).

    Examples:
    >>> [1, 2, [3,4], (5,6)]
    [1, 2, [3, 4], (5, 6)]
    >>> flatten([[[1,2,3], (42,None)], [4,5], [6], 7, MyVector(8,9,10)])
    [1, 2, 3, 42, None, 4, 5, 6, 7, 8, 9, 10]
    """

    result = []
    for el in x:
        #if isinstance(el, (list, tuple)):
        if hasattr(el, "__iter__") and not isinstance(el, basestring):
            result.extend(flatten(el))
        else:
            result.append(el)
    return result
################################################################################
def generate_fits_image(xy_matrx, file_name = 'out.fits'):
    """
    Function to generate a fits image from a matrix

    Input:
     - xy_matrx [[x1y1,x1y2...],[[x2y1,x2y2...]],... ] : file name
    Output:
     - NONE
    """


    hdu = pyfits.PrimaryHDU()
    hdu.data = xy_matrx

    hdu.writeto(file_name, clobber=True)
################################################################################
def fits_diff(image1, image2, file_name):
    file1 = pyfits.open(image1)
    data1 = file1[0].data


    file2 = pyfits.open(image2)
    data2 = file2[0].data
    diff = ( data1 - data2)#*data1[0][0]/data2[0][0])

    generate_fits_image(diff, file_name)

    print data1[0][0], data2[0][0], data1[0][0]/data2[0][0]
    print data1[0][0], data2[0][0]* data1[0][0]/data2[0][0]
################################################################################
def fits_diff_rel(image1, image2, file_name):
    file1 = pyfits.open(image1)
    data1 = file1[0].data


    file2 = pyfits.open(image2)
    data2 = file2[0].data
    diff = ( data1 - data2 ) / data1

    generate_fits_image(diff, file_name)

    print data1[0][0], data2[0][0], data1[0][0]/data2[0][0]
    print data1[0][0], data2[0][0]* data1[0][0]/data2[0][0]
################################################################################
def write_log(string_in, file_name):
    """
    Function to write log files. It opens or creates file_name and add
    string_in to the end of file
    Input
     - string_in str : string to be added to the end of file
     - file_name str : name of the log file
    Output
     - NONE
    #FIXME - finish documentation
    """
    import datetime

    file_out = open(file_name, "a")
    str_tmp = "## " + str(datetime.datetime.now()) + " ##\n"
    file_out.write(str_tmp)
    file_out.write(string_in)

    if ( not string_in[len(string_in) - 1] == "\n" ):
        file_out.write("\n")
    file_out.write("################################\n")
###############################################################################
def select_region(ra_vec, dec_vec, ra_ref, dec_ref, rad_arcmin):
    """
    Functions to select objects in a specifc region
    Input
     -
    Output
     -
    #FIXME - finish documentation
    """
    if (type(ra_vec) != np.array):
        ra_vec = np.array(ra_vec)
    if (type(dec_vec) != np.array):
        dec_vec = np.array(dec_vec)

    dec_min = dec_ref - rad_arcmin/60.0
    ra_min = (ra_ref - rad_arcmin/60.0)#*math.cos(dec_min)

    dec_max = dec_ref + rad_arcmin/60.0
    ra_max = (ra_ref + rad_arcmin/60.0)#*math.cos(dec_max)

    #print dec_min,dec_max
    #print ra_min, ra_max

    msk = (ra_vec >= ra_min) & (ra_vec <= ra_max) & \
          (dec_vec >= dec_min) & (dec_vec <= dec_max)

    ra_vec_msk = ra_vec[msk]
    dec_vec_msk = dec_vec[msk]

#    print ra_vec_msk, dec_vec_msk

#    print "fk5"

    #str_out = ""

    index_out = []

    index_out_np = np.array([])

    for i in range(len(ra_vec_msk)):
        #str_out = "circle(%.10f, %.10f, %.0f\")" % \
        #          (ra_vec_msk[i], dec_vec_msk[i], 3)
        #print str_out
        #print ra_vec_msk[i], dec_vec_msk[i]
        #print np.where( ra_vec == ra_vec_msk[i] )[0], \
        #      np.where( dec_vec == dec_vec_msk[i] )[0]

        index_out.append( np.where(ra_vec == ra_vec_msk[i])[0] )
        index_out_np = np.append(index_out_np, \
                                 np.where(ra_vec == ra_vec_msk[i])[0])

    #print len(index_out_np) - len(np.unique(index_out_np))
    return np.unique(index_out_np)
###############################################################################
def make_circunference(radius, x_0, y_0):
    """
    Function to generate points on a circunference
    Input
     - radius float : radius
     - x_0    float : central 'x' position
     - y_0    float : central 'y' position
    Output
     - x_vec, y_vec
    #FIXME - finish documentation
    """
    npts = 50
    theta = np.linspace(0.0, 2.0*math.pi, num = npts)
    x_out = radius*np.cos(theta) + x_0
    y_out = radius*np.sin(theta) + y_0
    return x_out, y_out
###############################################################################
def s82_ra(vec_ra):
    """
    Convert ra range from 0-360 to -180 to 180
    (useful when ploting S82 regions).
    Input
     - vec_ra : vector with ra coordinates
    Output
     - vec_ra converted
    #FIXME - finish documentation
    """
    #msk = (vec_ra > 180)
    for i in range(len(vec_ra)):
        if vec_ra[i] > 180:
            vec_ra[i] -= 360
    return vec_ra

