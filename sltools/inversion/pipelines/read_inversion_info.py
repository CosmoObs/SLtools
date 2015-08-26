#!/usr/bin/env python
# ==================================
# Author: Gabriel Bartosch Caminha - gbcaminha@gmail.com
# ==================================
"""
Package to test the function read_inversion_info
"""

import lenstool as lt

import numpy as np

import functions as fc

import parser as ps

import composite_model as cm

import os

import matplotlib.pyplot as plt

from matplotlib import rc

import re

rc('text', usetex = True)
FONT = {'family' : 'serif', \
        'size' : 15}

rc('font', **FONT)

def main_test_read_inversion_info():
    import argparse

    bool_sigma_eff = False

    parser = argparse.ArgumentParser(\
        description='Make plots from several lenstool inversions - UD!.', \
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('directory', default='invertion_stat_results', \
                      help='directory where lenstool inversion files are')

    args = parser.parse_args()

    args.directory = fc.add_slash(args.directory)

    number_files_single_inversion = 5

    config_file_name = 'gbclib_in.cfg'

    listing = np.array(os.listdir(args.directory))
    listing = np.delete( listing, np.where(listing == config_file_name) )
    listing.sort()


    #check if the number of files is consistent

    if len(listing) % number_files_single_inversion:
        print 'error: Number of files for a single inverson is !=', \
               number_files_single_inversion
        return 0

    listing = listing.reshape( len(listing)/number_files_single_inversion, \
                               number_files_single_inversion )


    dictionari_vec = []
    for i in listing:
        dictionari_vec.append( {
                'file_best_fit' : args.directory + i[0], \
                'file_chires' : args.directory + i[1], \
                'file_generate_arcs' : args.directory + i[2], \
                'file_make_inversion' : args.directory + i[3], \
                'file_source' : args.directory + i[4], \
               } )

    xi2_vec = []
    rmsi_vec = []
    rmss_vec = []
    forme_vec = []
    sigmav_eff = []


    tmp = lt.read_inversion_info( dictionari_vec[0] )
    opt_parameters = tmp['fited_parameters']

    dic_best_fit = {}
    dic_imput_vec = {}
    for i in opt_parameters:
        dic_best_fit[i] = []
        dic_imput_vec[i] = []

    args.input_file = args.directory + config_file_name
    ps.parser_file(args, verbose = False)

    #read the data from all files
    for i in dictionari_vec:
        tmp = lt.read_inversion_info(i)
        xi2_vec.append(tmp['xi2'])

        input_parameters = tmp['input_lens']
        for j in opt_parameters:
            dic_best_fit[j].append(float(tmp['best_fit_lens'][0][j][0]))
            dic_imput_vec[j].append(float(input_parameters[j][0]) )
        rmsi_vec.append(tmp['rmsi_mean'])
        rmss_vec.append(tmp['rmss_mean'])
        forme_vec.append(tmp['forme'])

        radius = cm.SisComp(i["file_best_fit"], args.source_real_z).sis_par[0]
        print "-system:", i["file_best_fit"][22:], "---------------------------"
        if bool_sigma_eff:
            lens_in = cm.SisComp(i["file_generate_arcs"], args.source_real_z)
            sigmav_eff_tmp = \
                     lens_in.proj2sigmav( lens_in.int_proj_kappa(radius), 1, 2 )
            sigmav_eff.append( sigmav_eff_tmp)
        #print float(tmp['best_fit_lens'][0][j][0]) - sigmav_eff_tmp
        #print sigmav_eff_tmp

    #print np.array(sigmav_eff) - dic_best_fit["v_disp"]
    #print dic_best_fit["v_disp"]
    #plt.hist( dic_best_fit["v_disp"] - np.array(sigmav_eff) + 500, range = (460, 540) )
    #plt.show()
    #plotting the results from several optimizations
    subplot_lines = 3
    subplot_num = 1

    if bool_sigma_eff:
        dic_best_fit["v_disp"] = \
                             dic_best_fit["v_disp"] - np.array(sigmav_eff) + 500

    plt.figure(1, figsize=(4*len(opt_parameters), \
                           len(opt_parameters)*subplot_lines))


    #Verify if the value for the image error are not 0.0 to put the correct \
    #indication on the plot 
    if args.position_error_value < 1E-5:
        args.position_error = False

    str_out = r'$\rm Substrucure : ' + str(args.substructure) + \
              ',\;Im Pos Err : ' + str(args.position_error) + \
              ',\;Lens Pos Err : ' + str(args.lens_position_error) + '$'
    plt.figtext( x = 0.005, y = 0.98, s = str_out)
    str_out = r'$\rm forme:' + forme_vec[0] + ',\; n\_systems = %i' % \
              (len(forme_vec)) + '$'
    plt.figtext( x = 0.005, y = 0.95, s = str_out)

    #plotting histograms
    for i in opt_parameters:
        plt.subplot(subplot_lines, len(opt_parameters), subplot_num)
        plt.hist( dic_best_fit[i] )
        plt.axvline( dic_imput_vec[i][0], color = 'red', linewidth = 1.25 )
        subplot_num = subplot_num + 1
        #plt.set_yticklabels([])
        mean_tmp = '%.2f' % np.mean(dic_best_fit[i])
        std_tmp = '%.2f' % np.std(dic_best_fit[i])
        #mean_tmp = str(np.mean(dic_best_fit[i]))
        #std_tmp =  str(np.std(dic_best_fit[i]))
        plt.title( r'\Tex $' + mean_tmp + '\pm' + std_tmp + '$')

    #plotting xi2 scatter
    plt.subplot(subplot_lines, len(opt_parameters), subplot_num)
    plt.ylabel(r'\Tex $\chi^2$')
    for i in opt_parameters:
        plt.subplot(subplot_lines, len(opt_parameters), subplot_num)
        plt.plot(dic_best_fit[i], xi2_vec, '.')
        plt.axvline( dic_imput_vec[i][0], color = 'red', linewidth = 1.5 )
        subplot_num = subplot_num + 1

    #plot the rms on the lens plane or on the source plane, depending on where
    #the optimizations was performed
    if forme_vec[0] in ['-11', '-1']:
        #plotting rmsi scatter
        plt.subplot(subplot_lines, len(opt_parameters), subplot_num)
        plt.ylabel(r'\Tex $rms_{lens}$')
        for i in opt_parameters:
            plt.subplot(subplot_lines, len(opt_parameters), subplot_num)
            plt.plot(dic_best_fit[i], rmsi_vec, '.')
            plt.xlabel(r'\Tex $\rm ' + re.sub('_', '\_', i) + '$')
            plt.axvline( dic_imput_vec[i][0], color = 'red', linewidth = 1.5 )
            subplot_num = subplot_num + 1
    else:
        #plotting rmss scatter
        plt.subplot(subplot_lines, len(opt_parameters), subplot_num)
        plt.ylabel(r'\Tex $rms_{source}$')
        for i in opt_parameters:
            plt.subplot(subplot_lines, len(opt_parameters), subplot_num)
            plt.plot(dic_best_fit[i], rmss_vec, '.')
            plt.xlabel(r'\Tex $\rm ' + re.sub('_', '\_', i) + '$')
            plt.axvline( dic_imput_vec[i][0], color = 'red', linewidth = 1.5 )
            subplot_num = subplot_num + 1

    filename = args.directory[:len(args.directory)-1] + '.pdf'
    plt.savefig(filename)
    plt.show()
    print 'File saved :', filename

if __name__ == '__main__':
    main_test_read_inversion_info()
