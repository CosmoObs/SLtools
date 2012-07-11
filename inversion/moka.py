#!/usr/bin/env python
# ==================================
# Author: Gabriel Bartosch Caminha - gbcaminha@gmail.com
# ==================================

import os

import numpy as np

import functions as fc

import random as rd

import shutil

import math

#change the matplotlib display to allow to run in a terminal with no X window
#FIXME - not tested yet
if not fc.is_X_running():
    import matplotlib
    matplotlib.use('agg')

import matplotlib.pyplot as plt

from matplotlib import rc
rc('text', usetex=True)
font = {'family' : 'serif', \
        'size' : 17}
rc('font', **font)


MOKA_INPUT = { \
              "!...Omega0" : 0.3, \
              "!...OmegaL" : 0.7, \
              "!...h0" : 0.7, \
              "!...w_dark_energy" : -1, \
              "!...lens_redshift" : 0.25, \
              "!...halo_mass" : 1.e15, \
              "!...source_redshift" : 2., \
              "!...field_of_view" : 0, \
              "!...JAFFE_HERNQUIST_or_DM_for_the_BCG_or_not" : "not", \
              "!...number_of_haloes" : 1, \
              "!...npixels" : 300, \
              "!...write_fits_for_what_component" : "YES", \
              "!...write_SkyLens_file_FORTRAN_or_CPP" : "FORTRAN", \
              "!...input_file_mass_concentration" : "NO", \
              "!...scatter_in_concentration_log_normal" : 0.25, \
              "!...mass_resolution_subhalo_mass_function" : 1.0e10, \
              "!...SISc_or_NFWc_for_satellitess" : "SIS", \
              "!...satellite_distribution_model_NFW_or_GAO" : "GAO", \
              "!...mstar_zl" : 0*5.21206e+12, \
              "!..elliptical_or_spherical_halo" : 1, \
              "!...n_pix_smooth_fourier" : 0, \
              "!..bcg_velocity_disp_prof_and_lensing_comps_prof" : "NO" \
             }

def write_input(moka_dic):
    """
    Function to write MOKA input file
    
    Input:
     - moka_dic dictionarie : dictionarie with all information for INPUT
    Output:
     - NONE
    """

    file_out = open('INPUT', 'w')

#Since moke INPUT file demands the parameters in a specific order, it is 
#more suitable make  each print hard-coded
    i = "!...Omega0"
    file_out.write(i + "\n" + str(moka_dic[i]) + "\n")

    i = "!...OmegaL"
    file_out.write(i + "\n" + str(moka_dic[i]) + "\n")

    i = "!...h0"
    file_out.write(i + "\n" + str(moka_dic[i]) + "\n")

    i = "!...w_dark_energy"
    file_out.write(i + "\n" + str(moka_dic[i]) + "\n")

    i = "!...lens_redshift"
    file_out.write(i + "\n" + str(moka_dic[i]) + "\n")

    i = "!...halo_mass"
    file_out.write(i + "\n" + str(moka_dic[i]) + "\n")

    i = "!...source_redshift"
    file_out.write(i + "\n" + str(moka_dic[i]) + "\n")

    i = "!...field_of_view"
    file_out.write(i + "\n" + str(moka_dic[i]) + "\n")

    i = "!...JAFFE_HERNQUIST_or_DM_for_the_BCG_or_not"
    file_out.write(i + "\n" + str(moka_dic[i]) + "\n")

    i = "!...number_of_haloes"
    file_out.write(i + "\n" + str(moka_dic[i]) + "\n")

    i = "!...npixels"
    file_out.write(i + "\n" + str(moka_dic[i]) + "\n")

    i = "!...write_fits_for_what_component"
    file_out.write(i + "\n" + str(moka_dic[i]) + "\n")

    i = "!...write_SkyLens_file_FORTRAN_or_CPP"
    file_out.write(i + "\n" + str(moka_dic[i]) + "\n")

    i = "!...input_file_mass_concentration"
    file_out.write(i + "\n" + str(moka_dic[i]) + "\n")

    i = "!...scatter_in_concentration_log_normal"
    file_out.write(i + "\n" + str(moka_dic[i]) + "\n")

    i = "!...mass_resolution_subhalo_mass_function"
    file_out.write(i + "\n" + str(moka_dic[i]) + "\n")

    i = "!...SISc_or_NFWc_for_satellitess"
    file_out.write(i + "\n" + str(moka_dic[i]) + "\n")

    i = "!...satellite_distribution_model_NFW_or_GAO"
    file_out.write(i + "\n" + str(moka_dic[i]) + "\n")

    i = "!...mstar_zl"
    file_out.write(i + "\n" + str(moka_dic[i]) + "\n")

    i = "!..elliptical_or_spherical_halo"
    file_out.write(i + "\n" + str(moka_dic[i]) + "\n")

    i = "!...n_pix_smooth_fourier"
    file_out.write(i + "\n" + str(moka_dic[i]) + "\n")

    i = "!..bcg_velocity_disp_prof_and_lensing_comps_prof"
    file_out.write(i + "\n" + str(moka_dic[i]) + "\n")

    file_out.close()


def run_moka():
    #file_seed_ell = open('rand1.d', 'w')
    file_seed_sub = open('rand3.d', 'w')

    if os.path.exists("info_haloes.dat"):
        os.remove("info_haloes.dat")
    if os.path.exists("nfw_info.dat"):
        os.remove("nfw_info.dat")
    if os.path.exists("sis_info.dat"):
        os.remove("sis_info.dat")
    if os.path.exists("fits"):
        shutil.rmtree("fits")
    if os.path.exists("fits_xy"):
        shutil.rmtree("fits_xy")
    if os.path.exists("fits_yz"):
        shutil.rmtree("fits_yz")
    if os.path.exists("fits_zx"):
        shutil.rmtree("fits_zx")

    os.system("time MOKA.3.0.0b")

    #file_seed_ell.close()
    file_seed_sub.close()

    #os.remove("rand1.d")
    os.remove("rand3.d")

def load_moka_data(plt_show = False):
    """
    Function to load informations from a single moka run
    
    Input:
     - lenses_vec [{},...] : list of dictionaries with the lens model
     - args      namespace : name space with float 'source_final_z'
     - plt_show        str : lenstool file name
    Output:
     - NONE
    """


    if os.path.exists("axes.d") and os.path.exists("projection.d"):
        print ".d files exist"

    nfw_info_path = "nfw_info.dat"
    nfw_data = { "kappa_s" : \
                         fc.extract_parameter(nfw_info_path, "kappa_s")[0][0], \
                 "cl" : fc.extract_parameter(nfw_info_path, "cl")[0][0] \
               }

    sis_info_path = "sis_info.dat"
    if os.path.exists(sis_info_path):
        sis_data_tmp = fc.flatten( fc.extract_parameter(sis_info_path, "cl") )
    else:
        sis_data_tmp = []

    sis_data = []
    i = 0
    while i < len(sis_data_tmp):
        sis_data.append(float(sis_data_tmp[i]))
        i = i + 2


    if os.path.exists(sis_info_path):
        stat_info_path = "satellites/satinfo.0.sub"
        stat_data = np.loadtxt(stat_info_path, unpack = True)
    else:
        stat_data = [[]]

    main_info_path = "info_haloes.dat"
    main_data = np.loadtxt(main_info_path, unpack = True)

    rvir0 = fc.extract_parameter( "rvir0.dat", "Rvir0")[0][0]

    if plt_show:
        plt.figure(1, figsize=(8, 8))
        plt.hist( np.log10(stat_data[0]), log = True, bins = 15 )
        plt.xlabel(r"$\rm \log_{10}(Mass)$")
        plt.xlim( min(np.log10(stat_data[0])), max(np.log10(stat_data[0])) )
        plt.ylabel(r"$\rm counts$")
        plt.title(r"$\rm Substructure \; Mass \; Distribution$")
        plt.show()

        plt.subplot(1, 1, 1).set_aspect(1)
        size_1 = stat_data[0]/max(stat_data[0])*50 + 1
        axis_limit = max( 1.075*max(stat_data[2]), 1.075*max(stat_data[3]) )

        plt.scatter(stat_data[2], stat_data[3], s = size_1 , color = "b", \
                    edgecolors = 'none')
        plt.axis([-axis_limit, axis_limit, -axis_limit, axis_limit])
        plt.xlabel("$x$")
        plt.ylabel("$y$")
        plt.show()

        #plt.hist(np.log10(sis_data), log = True)
        #plt.show()

    ell = 0.0

    file_gravlens = open('gravlens.gravlens', 'w')

    file_gravlens.write("startup " + str(len(stat_data[0]) + 1) + " 1\n")

    file_gravlens.write("    nfw " + \
                        str(float(nfw_data["kappa_s"])/(1.0 - ell)) + \
                        " 0.0 0.0 " + str(ell) + " 0.0 0.0 0.0 " + \
                         str(float(nfw_data["cl"])/(1.0 - ell)) + " 0.0 0.0\n")

    for i in range(len(sis_data)):
        file_gravlens.write( "    alphapot " + str(sis_data[i]) + " " + \
                   str(stat_data[2][i]) + " " + str(stat_data[3][i]*1.4) + \
                             " 0.0 0.0 0.0 0.0 0.0 0.0 1.0\n")

    for i in range(len(sis_data) + 1):
        file_gravlens.write( "    0 0 0 0 0 0 0 0 0 0\n")

    file_gravlens.write( "plotdef1 pot_gl.txt -" + rvir0 + " " + rvir0 + \
                         " 256 -" + rvir0 + " " + rvir0 + " 256\n" )

    #file_gravlens.write( "plotkappa kappa_gl.fits 3 -" + rvir0 + " " + rvir0 + \
    #                     " 512 -" + rvir0 + " " + rvir0 + " 512\n" )

    file_gravlens.write( "plotcrit crit.dat" )

    file_gravlens.close()


























