# ==================================
# Author: Gabriel Bartosch Caminha - gbcaminha@gmail.com
# ==================================
"""
Package with stand alone functions of gbclib
"""
import pyfits

import random
from math import sqrt


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
            out.append(line_inf[size_identifier:])
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
