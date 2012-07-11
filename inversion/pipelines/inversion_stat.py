#!/usr/bin/env python
# ==================================
# Author: Gabriel Bartosch Caminha - gbcaminha@gmail.com
# ==================================
"""Pipeline to make several lens inversion using lenstool"""
import numpy as np

import functions as fc

import lenstool as lt

import plot_curves as pc

import parser as ps

import os

import shutil

import glob

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
                                          vdisp_limit = 200 )

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
    str_log = ""
    if args.substructure:
    #create a par file with the main lens + substructures
        all_systems = add_subs_parfile(lens_file = args.lens_file, n_subs = 3)
    else:
    #create a par file with the main lens only
        all_systems = add_subs_parfile(lens_file = args.lens_file, n_subs = 0)

#    if args.lens_position_error:
#        for i in all_systems:
#            lens_err_x = np.random.randn() \
#                         * float(args.lens_position_error_value)
#            lens_err_y = np.random.randn() \
#                         * float(args.lens_position_error_value)

#            str_log = str_log + "lens_xy_real  " + str(i['x_centre']) + "  " + \
#                       str(i['y_centre']) + "\n"
#            i['x_centre'] = i['x_centre'] + lens_err_x
#            i['y_centre'] = i['y_centre'] + lens_err_y
#            str_log = str_log + "lens_xy_final  " + str(i['x_centre']) + \
#                      "  " + str(i['y_centre']) + "\n"

    str_log = str_log + "lens_xy_final  " + \
              str(all_systems[len(all_systems) - 1]["x_centre"]) + "  " + \
              str(all_systems[len(all_systems) - 1]["y_centre"]) + "\n"

    if args.lens_position_error:
        lens_err_x = np.random.randn() \
                     * float(args.lens_position_error_value)
        lens_err_y = np.random.randn() \
                     * float(args.lens_position_error_value)
        all_systems[len(all_systems) - 1]["x_centre"] = lens_err_x + \
                               all_systems[len(all_systems) - 1]["x_centre"]

        all_systems[len(all_systems) - 1]["y_centre"] = lens_err_y + \
                               all_systems[len(all_systems) - 1]["y_centre"]

    str_log = str_log + "lens_xy_final  " + \
              str(all_systems[len(all_systems) - 1]["x_centre"]) + "  " + \
              str(all_systems[len(all_systems) - 1]["y_centre"]) + "\n"



    #Adding error on the lens red-shift
    if args.lens_z_error:
        for i in all_systems:
            lens_err_z = np.random.randn() * float(args.lens_z_error_value)
            i['z_lens'] = i['z_lens'] + lens_err_z

    #Adding error on the source red-shift
    if args.source_z_error:
        source_err_z = np.random.randn() * float(args.source_z_error_value)
        args.source_final_z = args.source_real_z + source_err_z
        #z_source can be < z_lens after adding the z_source error, this while
        #avoid this
        while args.source_final_z < args.lens_real_z:
            source_err_z = np.random.randn() * float(args.source_z_error_value)
            args.source_final_z = args.source_real_z + source_err_z
    else:
        args.source_final_z = args.source_real_z

    #print all_systems
    lt.generate_par_images(all_systems, file_par = 'im_' + args.lens_file)
    lt.generate_par_curves(all_systems, args, file_par = 'cc_' + args.lens_file)
    #run lenstool to generate de caustics and critical curves
    str_run_lt = 'lenstool ' + 'cc_' + args.lens_file + ' -n'
    os.system(str_run_lt)

    #read the caustics and define the radius within the sources will be sorted
    courves_file = 'ce.dat'
    xy_ca = np.loadtxt( courves_file, usecols=(3, 4), unpack=True )
    ca_radial_max = np.sqrt(max(xy_ca[0]**2 + xy_ca[1]**2))

    lenstool_img_pos = []
    count = 0
    #n_inter_max = 200
    #sort the sources until find multiplicitie > 3
    while len(lenstool_img_pos) < 4: #and count < n_inter_max:
        source = fc.generate_sources(n_src = 1, radius = ca_radial_max)
        lt.write_src_multfile( src_coords = source, args = args, \
                               file_name = 'source.mul' )
        str_run_lt = 'lenstool ' + 'im_' + args.lens_file + ' -n > lenstool.out'
        os.system(str_run_lt)

        lenstool_img_pos = list(lt.read_images_lenstool(\
                                                      image_file = 'image.all'))
        count = count + 1
    #removing the image with lower magnification (usualy < 1) to disconsider
    #the central image
    if len(lenstool_img_pos) > 4:
        print "Found more than 4 multiple images"
        lenstool_img_mag = lt.read_amplification()
        min_index = lenstool_img_mag[0].index(min(lenstool_img_mag[0]))

        lenstool_img_mag[0].pop( min_index )
        lenstool_img_pos.pop( min_index )

    #saving multiple image positions
    for i in lenstool_img_pos:
        str_log = str_log + "im_xy_real  " + str(i[0]) + "  " + str(i[1]) + "\n"

    #adding error to the image positions
    if args.position_error:
        #lenstool_img_pos_shape = np.shape(lenstool_img_pos)
        im_err = float(args.position_error_value) * \
                 np.random.randn(np.shape(lenstool_img_pos)[0], \
                                 np.shape(lenstool_img_pos)[1])

        lenstool_img_pos = lenstool_img_pos + im_err

    for i in lenstool_img_pos:
        str_log = str_log + "im_xy_final  " + str(i[0]) + "  " + str(i[1]) + \
                  "\n"
    fc.write_log( str(inversion_id), "test.txt")
    fc.write_log( str_log, str(inversion_id) + "extra.txt" )


    #FIXME find a way to plot cc from real and optimized system
    if plt_show:
        pc.main_plot_curves_py(plot_show=False)
        lt.plot_points(images=lenstool_img_pos, plt_show=True)

    opt_result = lt.make_inversion_lenstool(lenstool_img_pos, \
            args.source_final_z, inversion_par_file = args.inv_file)


    #Saving files:
    file_dic = {
                'file_generate_arcs' : 'im_' + args.lens_file, \
                'file_source' : 'source.mul', \
                'file_make_inversion' : args.inv_file, \
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

    parser.add_argument("--inv-file", default = "make_inversion_siep.par", \
                        help = "Lentool file to make the inversion")

    args = parser.parse_args()
    ps.parser_file(args)

    if not os.path.isfile(args.input_file):
        print "error:", args.input_file, "does not exists"
        return 0

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

    print optimization_result

    #Renaming the folder with outputs
    args.ellipticite = lt.extract_parameter(args.lens_file, 'ellipticite')[0][0]
    args.v_disp = lt.extract_parameter(args.lens_file, 'v_disp')[0][0]
    run_id = 1
    dst_dir = str(run_id) + 'err' + str(args.position_error_value) + '_ell' + \
           str(args.ellipticite) + '_vdisp' + str(args.v_disp)

    while os.path.isdir(dst_dir):
        run_id = run_id + 1
        dst_dir = str(run_id) + 'err' + str(args.position_error_value) + '_ell' + \
               str(args.ellipticite) + '_vdisp' + str(args.v_disp)

    #print os.path.isdir(dst_dir)
    original_dir = 'invertion_stat_results'
    os.rename(original_dir, dst_dir)

    shutil.copyfile(args.input_file, dst_dir + '/' + args.input_file)
    os.makedirs(dst_dir + "_extra/")
    for i in glob.glob("*extra.txt*"):
        shutil.move(i, dst_dir + "_extra/")
    print 'Result saved at', dst_dir
################################################################################
# MAIN #########################################################################
################################################################################
if __name__ == '__main__':
    main_inversion_stat_py()
