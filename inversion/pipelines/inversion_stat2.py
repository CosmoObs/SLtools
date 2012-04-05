#!/usr/bin/env python
# ==================================
# Author: Gabriel Bartosch Caminha - gbcaminha@gmail.com
# ==================================
"""Pipeline to make several lens inversion using lenstool"""
import functions as fc

import lenstool as lt

import numpy as np

import plot_utils as pu

import plot_curves as pc

import parser as ps

import os

def add_subs_parfile(lens_file, n_subs):
    """
    Function to add substructures in a lenstool par file

    Input:
     - lens_file      str: par file containing the main lenses
     - lens_out_file  str: par file where the substructures are added
    Output:
     - NONE
    """
    print 'Adding substructures to:', lens_file
    subs_lens = fc.generate_substructure( z_lens = 1.0, n_subs = n_subs, \
                                          vdisp_limit = 150 )

    main_lens = fc.extract_second_identifiers(lens_file, 'potential')

    print '    The file', lens_file, 'has', len(main_lens), 'main lenses'
    print '    Adding', len(subs_lens), 'more lenses'
    
    all_lens = []
    all_lens.append(subs_lens)
    for i in range(len(main_lens)):
        for j in main_lens[i].keys():
            main_lens[i][j] = float(main_lens[i][j][0])
        all_lens[0].append(main_lens[i])
    
    #lt.generate_par_curves(all_lens[0], file_par = lens_out_file)
    
    #print '    Final system added at', lens_out_file
    return all_lens[0]
################################################################################
def make_single_inversion(args, plt_show = False, inversion_id = 1):
    """
    Function to simulate and make the inversion for a single strong lensing
    system

    Input:
     - lens_file     struct : arguments from the parser
     - plt_show        bool : flag to show the critical curves and caustics
     - inversion_id     int : id for the inversion
    Output:
     - NONE
    """
    if args.substructure:
    #create a par file with the main lens + substructures
        all_systems = add_subs_parfile(lens_file = args.lens_file, n_subs = 3)
        lt.generate_par_images(all_systems, file_par = 'im_' + args.lens_file)
        lt.generate_par_curves(all_systems, file_par = 'cc_' + args.lens_file)
    else:
    #create a par file only with the main lens
        all_systems = add_subs_parfile(lens_file = args.lens_file, n_subs = 0)
        lt.generate_par_images(all_systems, file_par = 'im_' + args.lens_file)
        lt.generate_par_curves(all_systems, file_par = 'cc_' + args.lens_file)

    #run lenstool first to generate de caustics and critical curves
    str_run_lt = 'lenstool ' + 'cc_' + args.lens_file + ' -n'
    os.system(str_run_lt)

    # read the caustics and define the radius within the sources will be sorted
    courves_file = 'ce.dat'
    xy_ca = np.loadtxt( courves_file, usecols=(3, 4), unpack=True )
    ca_radial_max = np.sqrt(max(xy_ca[0]**2 + xy_ca[1]**2))

    lenstool_img_pos = []
    #sort the sources until find multiplicitie > 3
    while len(lenstool_img_pos) < 4:
        source = fc.generate_sources(n_src = 1, radius = ca_radial_max)
        lt.write_src_multfile( src_coords = source, source_z = args.src_z, \
                               file_name = 'source.mul' )
        str_run_lt = 'lenstool ' + 'im_' + args.lens_file + ' -n > lenstool.out'
        os.system(str_run_lt)


        lenstool_img_pos = list(lt.read_images_lenstool(\
                                                      image_file = 'image.all'))

    #removing the image with lower magnification (usualy < 1) to disconsider
    #the central image
    if len(lenstool_img_pos) > 4:
        lenstool_img_mag = lt.read_amplification()
        min_index = lenstool_img_mag[0].index(min(lenstool_img_mag[0]))

        lenstool_img_mag[0].pop( min_index )
        lenstool_img_pos.pop( min_index )


    #adding error to the image positions
    if args.position_error:
        lenstool_img_pos_shape = np.shape(lenstool_img_pos)
        im_err = float(args.position_error_value) * \
                 np.random.randn(lenstool_img_pos_shape[0], \
                                 lenstool_img_pos_shape[1])
        lenstool_img_pos = lenstool_img_pos + im_err

    #FIXME find a way to plot cc from real and optimized system
    if plt_show:
        pc.main_plot_curves_py(plot_show=False)
        lt.plot_points(images=lenstool_img_pos, plt_show=True)

    opt_result = lt.make_inversion_lenstool(lenstool_img_pos, args.src_z, \
                                 inversion_par_file = 'make_inversion_siep.par')




    #Saving files:
    file_dic = {
                'file_generate_arcs' : 'im_' + args.lens_file, \
                'file_source' : 'source.mul', \
                'file_make_inversion' : 'make_inversion_siep.par', \
                'file_best_fit' : 'best.par', \
                'file_chires' : 'chires.dat', \
               }
    lt.copy_inversion_files(inversion_id, file_dic)
    return opt_result
    #str_os = 'rm ce.dat ci.dat'
    #os.system(str_os)
    


################################################################################
def main_inversion_stat_py():
    """
    Main function

    Input:
     - NONE 
    Output:
     - NONE
    """
    import argparse

    parser = argparse.ArgumentParser(\
        description='Make lens inversion using lenstool.', \
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument('--lens-file', default='generate_images.par', \
                      help='lenstool file used to generate the multiple images')

    parser.add_argument('--src-file', default='source.mul', \
                        help='lenstool file with sources information')

    parser.add_argument('--input-file', default='gbclib_in.cfg', \
                        help='configure file')

    parser.add_argument('--src-z', default=2.0, \
                        help='source red shift - NOT WORKING!')

    args = parser.parse_args()
    ps.parser_file(args)
    
    print '------------------------------------------------------------'
    print 'Module to make lens inversion for several simulated systems '
    print '------------------------------------------------------------'
    optimization_result = []
    optm_ell = []
    optm_vdisp = []
    optm_angle = []
    f_optm_out = open('optm_out.txt', 'w')
    for i in range(int(args.number_systems)):
        optimization_result.append( make_single_inversion(args, \
                                    plt_show = args.plot_curves, \
                                    inversion_id = i) )
        optm_ell.append( optimization_result[i]['ellipticite'] )
        optm_angle.append( optimization_result[i]['angle_pos'] )
        optm_vdisp.append( optimization_result[i]['v_disp'] )
        str_out = str(optm_ell[i]) + ' ' + str(optm_vdisp) + ' ' + \
                  str(optm_angle)
        f_optm_out.write(str_out)
    f_optm_out.close()
    #print optm_ell
    #print optm_angle
    #print optm_vdisp
    #pu.plot_scatter_hist(optm_ell, optm_vdisp)
    plt_fname = 'pos' + str(args.position_error_value) + '_e' + \
                str(args.number_systems)
    plt_fname = plt_fname + '.png'
    pu.plot_histogram(optm_ell, xlabel = r'$\varepsilon$', \
                      filename = plt_fname, \
                      plot=False, plt_range = (0, 0.6))

    plt_fname = 'pos' + str(args.position_error_value) + '_v' + \
                str(args.number_systems)
    plt_fname = plt_fname + '.png'
    pu.plot_histogram(optm_vdisp, xlabel = r'\Tex $\sigma$', \
                      filename = plt_fname, \
                      plot=False, plt_range = (450, 550))

    plt_fname = 'pos' + str(args.position_error_value) + '_t' + \
                str(args.number_systems)
    plt_fname = plt_fname + '.png' 
    pu.plot_histogram(optm_angle, xlabel = r'\Tex $\theta$', \
                      filename = plt_fname, \
                      plot=False, plt_range = (0, 80))
################################################################################
# MAIN #########################################################################
################################################################################
if __name__ == '__main__':
    main_inversion_stat_py()
