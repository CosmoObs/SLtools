# ==================================
# Author: Gabriel Bartosch Caminha - gbcaminha@gmail.com
# ==================================
"""
Package to automate some aplocations of lenstool
- to substitute lenstoolpypeline.py
"""

import os
import numpy as np

import matplotlib.pyplot as plt

from functions import extract_parameter

def generate_par_curves(lenses_vec, file_par = 'generate_curves.par'):
    """
    Generate a lenstool file to compute the tangential caustic and critical 
    curve
    
    Input:
     - lenses_vec [{},...] : list of dictionaries with the lens model
     - plt_show        str : lenstool file name   
    Output:
     - NONE
    """
    file_par_open = open(file_par, 'w')

    file_par_open.write("runmode\n"                  )
    file_par_open.write("        reference   3 0 0\n")
    file_par_open.write("        verbose     0\n"    )
    file_par_open.write("    	source      1 source.mul\n")
    file_par_open.write("        end\n"              )
    file_par_open.write("image\n"                    )
    file_par_open.write("        grid        0\n"    )
    file_par_open.write("        end\n"              )
    file_par_open.write("grille\n"                   )
    file_par_open.write("        nombre      128\n"  )

    string_out = "        nlens       " + str(len(lenses_vec)) + "\n"
    file_par_open.write(string_out)

    file_par_open.write("        nlens_crit  1\n")
    file_par_open.write("        nlens_opt   0\n")
    file_par_open.write("        polaire     1\n")
    file_par_open.write("        end\n"          )

    for i in range(len(lenses_vec)):
        string_out = 'potential ' + str(i) + '\n'
        file_par_open.write(string_out)
        for keys in lenses_vec[i].keys():
            string_out = '        ' + keys + ' ' + str(lenses_vec[i][keys]) + \
                         '\n' 
            file_par_open.write(string_out)
        file_par_open.write('        end\n')

    file_par_open.write("cosmology\n"                )
    file_par_open.write("        H0        70.0\n"   )
    file_par_open.write("        omega     0.3\n"    )
    file_par_open.write("        lambda    0.7\n"    )
    file_par_open.write("        end\n"              )
    file_par_open.write("cline\n"                    )
    file_par_open.write("        nplan 1 2.0\n"      )
    file_par_open.write("        limitHigh 0.2\n"   )
    file_par_open.write("        limitLow 0.1\n"     )
    file_par_open.write("#        algorithm snake\n"  )
    file_par_open.write("        end\n"              )
    file_par_open.write("champ\n"                    )
    file_par_open.write("        xmin -201\n"        )
    file_par_open.write("        xmax 200\n"         )
    file_par_open.write("        ymin -201\n"        )
    file_par_open.write("        ymax 200\n"         )
    file_par_open.write("        dmax 20.\n"         )
    file_par_open.write("        end\n"              )
    file_par_open.write("curve\n"                    )
    file_par_open.write("        nfile 2\n"          )
    file_par_open.write("        namein 1 0 ce.dat\n")
    file_par_open.write("        column 1 5\n"       )
    file_par_open.write("               2 3\n"       )
    file_par_open.write("        namein 2 0 ci.dat\n")
    file_par_open.write("        column 2 3\n"       )
    file_par_open.write("               2 3\n"       )
    file_par_open.write("        end\n"              )
    file_par_open.write("fini\n"                     )
    file_par_open.close()
################################################################################
def generate_par_images(lenses_vec, file_par = 'generate_images.par'):
    """
    Generate a lenstool file to compute the images from a source position
        
    Input:
     - lenses_vec [{},...] : list of dictionaries with the lens model
     - plt_show        str : lenstool file name   
    Output:
     - NONE
    """
    file_par_open = open(file_par, 'w')

    file_par_open.write("runmode\n"                        )
    file_par_open.write("        reference   3 0 0\n"      )
    file_par_open.write("        verbose     0\n"          )
    file_par_open.write("    	source      1 source.mul\n")
    file_par_open.write("        end\n"                    )
    
    file_par_open.write("image\n"                  )
    file_par_open.write("        grid        0\n"  )
    file_par_open.write("        n_source    100\n")
    file_par_open.write("        end\n"            )

    file_par_open.write("source\n"                 )
    file_par_open.write("        grid        500\n")
    file_par_open.write("        end\n"            )
    
    file_par_open.write("grille\n"                 )
    file_par_open.write("        nombre      128\n")
    
    string_out = "        nlens       " + str(len(lenses_vec)) + "\n"
    file_par_open.write(string_out)

    file_par_open.write("        nlens_crit  1\n"  )
    file_par_open.write("        nlens_opt   0\n"  )
    file_par_open.write("        polaire     1\n"  )
    file_par_open.write("        end\n"            )
    
    for i in range(len(lenses_vec)):
        string_out = 'potential ' + str(i) + '\n'
        file_par_open.write(string_out)
        for keys in lenses_vec[i].keys():
            string_out = '        ' + keys + ' ' + str(lenses_vec[i][keys]) + \
                         '\n' 
            file_par_open.write(string_out)
        file_par_open.write('        end\n')


    file_par_open.write("cosmology\n"             )
    file_par_open.write("        H0        70.0\n")
    file_par_open.write("        omega     0.3\n" )
    file_par_open.write("        lambda    0.7\n" )
    file_par_open.write("        end\n"           )
    
    file_par_open.write("champ\n"            )
    file_par_open.write("        xmin -101\n")
    file_par_open.write("        xmax 100\n" )
    file_par_open.write("        ymin -101\n")
    file_par_open.write("        ymax 100\n" )
    file_par_open.write("	     dmax 10.\n" )
    file_par_open.write("        end\n"      )

    file_par_open.write("fini\n")
    file_par_open.close()
################################################################################
def compute_curves_lenstool(lenses_vec, file_par = 'generate_curves.par'):
    """
    Compute the tangential caustics and critical curve using lenstool
    
    Input:
     - lenses_vec [{},...] : list of dictionaries with the lens model
     - plt_show        str : lenstool file name used to create the curves   
    Output:
     - NONE
    """
    generate_par_curves(lenses_vec, file_par)
    string_out = "lenstool " + file_par + " -n"
    os.system(string_out)
################################################################################
def plot_curves_lenstool(label = ' ', courves_file = 'ce.dat', marker = 'm.', \
                         plt_show = False):
    """
    Function to plot caustics and critical curves only from a lenstool 
    output file
    
    Input:
     - courves_file  str : file that contains the curves (lenstool output)
     - plt_show     bool : if true show the plot in a window
    Output:
     - NONE
    """
    x_ca, y_ca = np.loadtxt( courves_file, usecols=(3, 4), unpack=True )
    x_cc, y_cc = np.loadtxt( courves_file, usecols=(1, 2), unpack=True )

    plt.figure(1, figsize=(16, 8))

    plt.subplot(1, 2, 1).set_aspect(1)
    plt.plot(x_ca, y_ca, marker, linewidth=3, label = label)
    axis_limit = max( 1.075*max(np.absolute(x_ca)), \
    1.075*max(np.absolute(y_ca)) )
    plt.axis([-axis_limit, axis_limit, -axis_limit, axis_limit])

    plt.subplot (1, 2, 2).set_aspect(1)
    plt.plot(x_cc, y_cc, marker, linewidth=3, label = label)
    axis_limit = max( 1.075*max(np.absolute(x_cc)), \
    1.075*max(np.absolute(y_cc)) )
    plt.axis([-axis_limit, axis_limit, -axis_limit, axis_limit])
    plt.legend()

    if plt_show:
        plt.show()
        plt.close()
################################################################################
def plot_points(source = [], images = [], plt_show = False, \
                color = 'darkgreen', marker = 'x', label = ''): #FIXME - this function should me moved to a more generic module (plot_utils.py for instance)
    """
    Function to plot caustics and critical curves only from a lenstool 
    output file
    
    Input:
     - source   [[x,y],...] : vector with pairs to be plotet on the source plane
     - images   [[x,y],...] : vector with pairs to be plotet on the lens plane
     - plt_show        bool : if true show the plot in a window
    Output:
     - NONE
    """
    print source, images
    if len(source) == 0 and len(images) == 0:
        return 0

    
    plt.figure(1, figsize=(16, 8))

    if len(source) > 0:
        x_src = [row[0] for row in source]
        y_src = [row[1] for row in source]
        plt.subplot(1, 2, 1).set_aspect(1)
        plt.plot(x_src, y_src, markersize = 7, color = color, marker = marker, \
                    label = label, linewidth=0)
        plt.legend()

    if len(images) > 0:
        x_img = [row[0] for row in images]
        y_img = [row[1] for row in images]
        plt.subplot (1, 2, 2).set_aspect(1)
        plt.plot(x_img, y_img, markersize = 7, color = color, marker = marker, \
                    label = label, linewidth=0)
        plt.legend()

    if plt_show:
        plt.show()
        plt.close()
################################################################################
def generate_arcs_lenstool(source, lenses_vec, z_source = 2.0):
    """
    Function to compute the images from a list of point sources
    
    Input:
     - source     [[x,y],...] : list with source positions
     - lenses_vec    [{},...] : list of dictionaries with the lens model
     - z_source         float : source redshift
    Output:
     - NONE
    """
    print 'generate_arcs_lenstool in implementation'

    x_list = [row[0] for row in source]
    y_list = [row[1] for row in source]
    
    file_sources = open('source.mul', 'w') #FIXME - give the possibilite to change the file name
    file_sources.write('#REFERENCE 3 0.0000000 0.0000000 \n')

    for i in range(len(x_list)):
        str_out = str(i) + ' ' + str(x_list[i]) + ' ' + str(y_list[i]) + \
                  ' 0.000000 0.000000 0.00000 ' + str(z_source) + ' 0.00' + '\n'
        file_sources.write(str_out)
    file_sources.close()

    generate_par_images(lenses_vec, file_par = 'generate_images.par')#FIXME - give the possibilite to change the file name

    str_sys = 'lenstool ' + 'generate_images.par' + ' -n' #FIXME- see just abovegive the possibilite to change the file name
    os.system (str_sys)

    file_sources.close()
    return read_images_lenstool()
################################################################################
def read_images_lenstool(image_file = 'image.all'):
    """
    Function to split the multiples images from a lenstool image file

    Input:
     - image_file   str : file with images
    Output:
     - np.array   float : array with splited multiple image systems
    """
    #FIXME - is it a good idea to keep this one line function?
    x_img = np.loadtxt(image_file, usecols = (1, 2))
    return x_img
################################################################################
def make_inversion_lenstool(image_positons, z_source, \
                                inversion_par_file = 'make_inversion_siep.par'):
    """
    Function to perform le lens inversion for point images using lenstool for 
    siep model

    Input:
     - image_positons [[x,y],...] : a list with multiple image positions
     - source_z              float : source (image) redshift
     - inversion_par_file     str : par file to perform the inversion
    Output:
     - dictionarie with parameters: ellipticite, angle_pos, v_disp
    """
    write_multfile(image_positons, z_source, file_name = 'multfile.in')
    
    str_sys = 'lenstool ' + inversion_par_file + ' -n > lenstool.out'
    os.system(str_sys)
    
    ell = float( extract_parameter('best.par', 'ellipticite')[0][0] )
    angle = float( extract_parameter('best.par', 'angle_pos')[0][0] )
    vdisp = float( extract_parameter('best.par', 'v_disp')[0][0] )

    return {'ellipticite' : ell, 'angle_pos' : angle, 'v_disp' : vdisp}
################################################################################
def write_multfile(image_coords, source_z, file_name = 'multfile.in'):
    """
    Function to write a image file for lenstool use as imput for the lens
    inversion

    Input:
     - image_coords [[x,y],...] : array with image coordinates
     - source_z           float : source (image) redshift
     - file_name            str : file name to write the images
    Output:
     - NONE
    """
    file_in = open(file_name, 'w')
    file_in.write('#REFERENCE 3  0.0 0.0\n')

    for i in range(len(image_coords)):
        image_id = 'A' + str(i+1) + ' '
        data     = str(image_coords[i][0]) + ' ' + str(image_coords[i][1]) \
                   + str(' 0.2 0.2 0 ') + str(source_z) + ' 0'
        final    = image_id + data + '\n'
        file_in.write(final)
    file_in.close()
################################################################################
def write_src_multfile(src_coords, source_z, file_name = 'source.mul'):
    """
    Function to write a image file for lenstool use as imput for the lens
    inversion

    Input:
     - image_coords [[x,y],...] : array with image coordinates
     - source_z           float : source (image) redshift
     - file_name            str : file name to write the images
    Output:
     - NONE
    """
    file_in = open(file_name, 'w')
    file_in.write('#REFERENCE 3  0.0 0.0\n')

    for i in range(len(src_coords)):
        image_id = str(i+1) + ' '
        data     = str(src_coords[i][0]) + ' ' + str(src_coords[i][1]) \
                   + str(' 0.2 0.2 0 ') + str(source_z) + ' 0'
        final    = image_id + data + '\n'
        file_in.write(final)
    file_in.close()
################################################################################
def read_amplification(amp_file = 'dist.dat'):
    """
    Function to split the multiples images from a lenstool image file

    Input:
     - image_file   str : file with images
    Output:
     - np.array   float : array with amplification for each image
    """
    n_img, amp_img = np.loadtxt(amp_file, usecols=(0, 6), unpack=True)

    amp = []

    amp_tmp = []

    count = 1

    for i in range(len(n_img)):
        if count == n_img[i]:
            amp_tmp.append( amp_img[i] )
        else:
            amp.append(amp_tmp)

            amp_tmp = []

            amp_tmp.append( amp_img[i] )

            count = count + 1
    amp.append(amp_tmp)

    return amp