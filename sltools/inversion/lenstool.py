# ==================================
# Author: Gabriel Bartosch Caminha - gbcaminha@gmail.com
# ==================================
"""
Package to automate some aplocations of lenstool
- to substitute lenstoolpypeline.py
"""

import os

import shutil

import numpy as np

import functions as fc

import astropy.io.fits as pyfits

import math

#change the matplotlib display to allow to run in a terminal with no X window
#FIXME - not tested yet
if not fc.is_X_running():
    import matplotlib
    matplolib.use('agg')

import matplotlib.pyplot as plt

from functions import extract_parameter #FIXME use: import functions as fc



def generate_par_curves(lenses_vec, args, file_par = 'generate_curves.par'):
    """
    Generate a lenstool file to compute the tangential caustic and critical 
    curve
    
    Input:
     - lenses_vec [{},...] : list of dictionaries with the lens model
     - args      namespace : name space with float 'source_final_z'
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
    file_par_open.write("        nplan 1 " + str(args.source_final_z) + "\n")
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
    
    file_par_open.write("champ\n"              )
    file_par_open.write("        xmin   -201\n")
    file_par_open.write("        xmax   200\n" )
    file_par_open.write("        ymin   -201\n")
    file_par_open.write("        ymax   200\n" )
    file_par_open.write("	     dmax   50.\n" )
    file_par_open.write("        end\n"        )

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
def plot_cc_lenstool(label = "", courves_file = "ce.dat", marker = "m.", \
                     plt_show = False):
    """
    Function to plot critical curves from a lenstool output file
    Input
     - label         str : points label
     - courves_file  str : file that contains the curves (lenstool output)
     - marker        str : marker for the plot
     - plt_show     bool : if true show the plot in a window

    Output
     - NONE
    """
#    plt.figure(1, figsize=(8, 8))
    x_cc, y_cc = np.loadtxt( courves_file, usecols=(1, 2), unpack=True )
    plt.plot(x_cc, y_cc, marker, linewidth=3, label = label)
    axis_limit = max( 1.075*max(np.absolute(x_cc)), \
    1.075*max(np.absolute(y_cc)) )
    plt.axis([-axis_limit, axis_limit, -axis_limit, axis_limit])
    plt.legend()
    plt.axes().set_aspect("equal")

    if plt_show:
        plt.show()
        plt.close()


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
#    x_ca, y_ca = np.loadtxt( courves_file, usecols=(3, 4), unpack=True )
#    x_cc, y_cc = np.loadtxt( courves_file, usecols=(1, 2), unpack=True )

    x_cc, y_cc, x_ca, y_ca = \
       np.loadtxt( courves_file, usecols=(1, 2, 3, 4), unpack=True )

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
    #print source, images
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

    str_sys = 'lenstool ' + 'generate_images.par' + ' -n' #FIXME- see just above give the possibilite to change the file name
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
    Function to write a image file for lenstool use as input for the lens
    inversion

    Input:
     - image_coords [[x,y],...] : array with image coordinates
     - source_z           float : source (image) redshift
     - file_name            str : file name to write the images
    Output:
     - NONE
    """
    print 'write_multfile'
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
def write_src_multfile(src_coords, args, file_name = 'source.mul',):
    """
    Function to write a image file for lenstool use as input for the lens
    inversion

    Input:
     - image_coords [[x,y],...] : array with image coordinates
     - source_z           float : source (image) redshift
     - file_name            str : file name to write the images
    Output:
     - NONE
    """
    print 'write_src_multfile'
    file_in = open(file_name, 'w')
    file_in.write('#REFERENCE 3  0.0 0.0\n')

    for i in range(len(src_coords)):
        image_id = str(i+1) + ' '
        data     = str(src_coords[i][0]) + ' ' + str(src_coords[i][1]) \
                   + str(' 0.2 0.2 0 ') + str(args.source_real_z) + ' 0'
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
################################################################################
def read_inversion_info(file_dic):
    """
    Function to read the information from the lens inversion

    Input:
     - 
    Output:
     - 
    """
    #print_file_test = open('file_test.txt','w')

    if not ( check_inversion_files(file_dic) ):
        print 'error(read_inversion_info): problem with lenstool file names'
        return 0
    
    file_generate_arcs = file_dic['file_generate_arcs']
    info_input_lens = fc.extract_second_identifiers( file_generate_arcs, \
                                                     'potential' )
#-------------------------------------------------------------------------------

    file_source = file_dic['file_source']
    info_src = np.loadtxt(file_source, unpack=False)
    if len(info_src) == 8 and np.isscalar(info_src[0]):
    #FIXME - check if the second condition is all we need
        info_src = [info_src]
#-------------------------------------------------------------------------------

    file_make_inversion = file_dic['file_make_inversion']
    info_fited_param = fc.extract_second_identifiers( file_make_inversion, \
                                                      'limit' )
    info_forme = fc.extract_parameter(file_make_inversion, 'forme')[0][0]

#-------------------------------------------------------------------------------

    file_best_fit = file_dic['file_best_fit']
    info_best_fit = fc.extract_second_identifiers( file_best_fit, \
                                                   'potentiel' )

    info_xi2 = fc.extract_parameter(file_best_fit, '#Chi2pos:')

#-------------------------------------------------------------------------------
    file_chires = file_dic['file_chires']

    info_chires = extract_parameter(file_chires, '0')
    rmss_mean = [0.0, 0.0]
    rmsi_mean = [0.0, 0.0]
    for i in info_chires:
        if i[0] != 'A':
            rmss_mean[0] = rmss_mean[0] + float(i[7])
            rmss_mean[1] = rmss_mean[1] + 1.0
            
            rmsi_mean[0] = rmsi_mean[0] + float(i[8])
            rmsi_mean[1] = rmsi_mean[1] + 1.0

    rmss_mean = rmss_mean[0]/rmss_mean[1]
    rmsi_mean = rmsi_mean[0]/rmsi_mean[1]
#-------------------------------------------------------------------------------
    out_dict = { 'xi2' : float(info_xi2[0][0]), \
                 'best_fit_lens' : info_best_fit, \
                 'rmsi_mean' : rmsi_mean, \
                 'rmss_mean' : rmss_mean, \
                 'fited_parameters' : info_fited_param[0].keys(), \
                 'input_lens' : info_input_lens[len(info_input_lens) - 1], \
                 'forme' : info_forme \
               }
    #for i in out_dict.keys():
    #    print i, out_dict[i]
    return out_dict
################################################################################
def copy_inversion_files(count, file_dic):
    """
    Function to copy files created by lenstool when doing the lens inversion
    
    Input:
     - count    str : srting to be added to the end of all copied files
     - file_dic dic : dictionarie with file names created by lenstool inversion
    Output:
     - NONE
    """
    if not ( check_inversion_files(file_dic) ):
        print 'error(copy_inversion_files): problem with lenstool file names'
        return 0
    results_directory = 'invertion_stat_results'
    if results_directory[len(results_directory)-1] != '/':
        results_directory = results_directory + '/'
    #print results_directory

    file_generate_arcs = file_dic['file_generate_arcs']
    file_source = file_dic['file_source']
    file_make_inversion = file_dic['file_make_inversion']
    file_best_fit = file_dic['file_best_fit']
    file_chires = file_dic['file_chires']

    file_generate_arcs_out = results_directory + str(count) + file_generate_arcs
    file_source_out = results_directory + str(count) + file_source
    file_make_inversion_out = results_directory + str(count) + \
                              file_make_inversion
    file_best_fit_out = results_directory + str(count) + file_best_fit
    file_chires_out = results_directory + str(count) + file_chires

    if not os.path.exists(results_directory):
        os.makedirs(results_directory)
    else:
        print 'Directory', results_directory, 'already exists'

    shutil.copyfile(file_generate_arcs, file_generate_arcs_out)
    shutil.copyfile(file_source, file_source_out)
    shutil.copyfile(file_make_inversion, file_make_inversion_out)
    shutil.copyfile(file_best_fit, file_best_fit_out)
    shutil.copyfile(file_chires, file_chires_out)
################################################################################
def check_inversion_files(file_dic):
    """
    Check if 'file_dic' has all file names created by lenstool inversion
    
    Input:
     - file_dic dic : dictionarie with file names created by lenstool inversion
    Output:
     - bool : True if everything is ok, False otherwise 
    """
    patern_dic = {'file_generate_arcs' : "", \
                  'file_source' : "", \
                  'file_make_inversion' : "", \
                  'file_best_fit' : "", \
                  'file_chires' : ""}

    for i in patern_dic.keys():
        os.path.isfile(file_dic[i])
        if not ( i in file_dic.keys() ):
            print 'error(check_inversion_files): check_inversion_files: item', \
                  i, 'is not defined in the input dictionary'
            return False
        if not os.path.isfile(file_dic[i]):
            print 'error(check_inversion_files): check_inversion_files: file', \
                  file_dic[i], 'does not exists'
            return False
    return True
################################################################################
def compute_kappa_map(lens_vec, size, size_map):
    """
    docstring for compute_kappa_map
    Input
     -
    Output
     -
    #FIXME - finish documentation
    """

    par_file_name = "kappa_map.par"
    fit_file_name = "kappa_map.fits"
    z_source = 2.0
    size_map = size_map * 1.05

    file_map = open(par_file_name, 'w')

    conv_lens_vec(lens_vec)

    file_map.write("runmode\n"                  )
    file_map.write("        reference   3 0 0\n")
    file_map.write("        verbose     0\n"    )
    file_map.write("        mass        3 " + str(size) + " " + \
                   str(lens_vec[0]["z_lens"]) + " " + fit_file_name + "\n")
    file_map.write("        end\n")
    file_map.write("source\n")
    file_map.write("    z_source     " + str(z_source) + "\n")
    file_map.write("    end\n")
    file_map.write("grille\n")
    file_map.write("        nombre      128\n")
    file_map.write("        nlens       4\n")
    file_map.write("        nlens_crit  1\n")
    file_map.write("        nlens_opt   0\n")
    file_map.write("        polaire     1\n")
    file_map.write("        end\n")


    for i in range(len(lens_vec)):
        string_out = 'potential ' + str(i) + '\n'
        file_map.write(string_out)
        #print string_out,
        for keys in lens_vec[i].keys():
            string_out = '        ' + keys + ' ' + str(lens_vec[i][keys]) + \
                     '\n'
            #print string_out,
            file_map.write(string_out)
        file_map.write('        end\n')

    file_map.write("cosmology\n")
    file_map.write("        H0        70.0\n")
    file_map.write("        omega     0.3\n")
    file_map.write("        lambda    0.7\n")
    file_map.write("        end\n")
    file_map.write("champ\n")
    file_map.write("        xmin -101\n")
    file_map.write("        xmax 100\n")
    file_map.write("        ymin -101\n")
    file_map.write("        ymax 100\n")
    file_map.write("        dmax " + str(size_map) + "\n")
    file_map.write("        end\n")
    file_map.write("fini\n")

    file_map.close()
################################################################################
def conv_lens_vec(lens_vec):
    """
    Function to substitute srting for float in lens models read using
    extract_second_identifiers
    Input
     - lens_vec   vec : lens vector
    Output
     -
    #FIXME - finish documentation
    """
    for i_lens in range(len(lens_vec)):
        #print i_lens
        for i_prop in lens_vec[i_lens]:
            lens_vec[i_lens][i_prop] = float(lens_vec[i_lens][i_prop][0])
    #print lens_vec
################################################################################
def compute_projmass(args):
    """
    Main function
    Input
     -
    Output
     -
    #FIXME - finish documentation
    """
    radius = args.radius/3600.0

    k_map = pyfits.open(args.kappa_map)
    k_data = k_map[0].data
    k_data_tmp = k_data

    pix_dim = math.fabs(k_map[0].header["CDELT1"])
    pix_unit = k_map[0].header["CUNIT1"]
    shape = k_map[0].data.shape

    x_axis = np.linspace(-(shape[0] - 1.0)/2.0*pix_dim , \
                          (shape[0] - 1.0)/2.0*pix_dim, shape[0])
    y_axis = np.linspace(-(shape[1] - 1.0)/2.0*pix_dim , \
                          (shape[1] - 1.0)/2.0*pix_dim, shape[1])

    if pix_unit != "deg":
        print "Error, pixel unit not in deg"
    if (x_axis.max() - x_axis.min())/2.0 < radius:
        print "Error, the radius is larger than the image limits"


    proj_mass = 0.0
    for i_x in range(shape[0]):
        for i_y in range(shape[1]):
            if x_axis[i_x]**2.0 + y_axis[i_y]**2.0 <= radius**2.0:
                #k_data_tmp[i_x][i_y] = 0.0
                proj_mass += k_data_tmp[i_x][i_y]

    print "%e M_sol" % (proj_mass*1E12)

    if args.plot_cont:
        circ = fc.make_circunference(radius*3600, 0, 0)
        plt.plot(circ[0], circ[1], "k--", linewidth = 2)
        plt.contour(x_axis*3600.0, y_axis*3600.0, k_data)
        plt.show()

    return proj_mass


################################################################################