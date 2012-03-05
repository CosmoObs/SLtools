#!/usr/bin/env python
# ====================================================
# Authors:
# Bruno Moraes - bruno.a.l.moraes@gmail.com - 01/Nov/2011
# ====================================================


"""Module to calculate the limiting magnitude for a given FITS catalogue"""

##@package mag_pipeline
#
# This package calculates the limiting magnitude for a given FITS catalogue.
# It has four mandatory arguments: the (input) FITS catalogue path and name,
# the magnitude field name to be used for the calculations (i.e. MAG_AUTO,
# MAG_APER, etc...), the corresponding error field name (i.e. MAGERR_AUTO,
# MAGERR_APER, etc...) and the inferior magnitude value for which a fit of
# of the object number count per magnitude bin should be performed. 
# For a detailed description of the procedure,
# check http://twiki.linea.gov.br/blabla".
# 
# The code can run on an arbitrary list of input FITS catalogues. The outputs
# are a check plot of the number count fitting procedure and an ASCII file,
# with each line corresponding to a FITS catalogue and containing information 
# such as the center coordinates of the image, the mean seeing, the total
# number of objects, the number of galaxies and stars per arcmin^2, the
# limiting magnitude, etc. This output can be used as input for a QA analysis
# of the collection of catalogues.
#
# Executable package: YES
#
# To see the package help message, type:
#
# > python mag_pipeline.py --help
#
# To run the code:
#
# > python mag_pipeline.py 'input_file' 'mag_field_name' 'mag_error_field_name' mag_inf_value plot_x_lim_inf plot_x_lim_sup



from __future__ import division
from ast import literal_eval as lit
import sys
import os
import errno
import numpy as np
import math
import cmath
from scipy import optimize
import collections
import fit_mag as fmg
import glob
import pyfits


# Hard-coded use of sltools functions. Need to make a build of the library to get rid of this path.

sys.path.append("/home/brunomor/Documents/Trabalho/Repositorio") ### Append Bruno
from sltools.catalog import fits_data as fd
from sltools.coordinate import wcs_conversion as wcscv

# =================================================================================================================

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST:
            pass
        else: raise


def mag_pipeline(input_file, field_names, folder_path, mag_inf, bin_size, gal_cut, mag99, S2N_cut, stell_idx, plot_inf, plot_sup):

    """
    Calculate the limiting magnitude for a given FITS catalogue.

    This function performs specified cuts to the catalogue, proceeds to bin 
    the data and then fit it to an analytical (number count) x (magnitude)
    function. A difference of more than 10% between the data and the fit
    defines the limiting magnitude. This value and several other catalogue
    properties are returned as a list of strings. Additionally, a check plot
    of the fitting procedure is saved to folder_path.


    Input:
     - input_file          str : Input file name
     - field_names  [str, str] : Magnitude and magnitude error field names
     - folder_path         str : Folder path
     - mag_inf           float : Inferior magnitude to perform fitting
     - bin_size          float : (Uniform) binning size
     - gal_cut            bool : Galaxy or Star fit (according to stell_idx cut)
     - mag99              bool : Inclusion of unidentified magnitudes
     - S2N_cut           float : Signal-to-noise minimum cut
     - stell_idx         float : CLASS_STAR value (defines what is a Galaxy/Star)
     - plot_inf          float : Check plot inferior x-axis value
     - plot_sup          float : Check plot superior x-axis value

    Output:
     - [str,...,str] : a list of the following properties of the catalog in
                       string format: catalogue name, original image center RA,
                       original image center DEC, original image seeing, fitting 
                       bin size, S/N cut value, CLASS_STAR value, constant 
                       coefficient of the linear fit, linear coefficient of the
                       linear fit, inferior magnitude of the performed fit,
                       superior magnitude of the performed fit, limiting
                       magnitude, total number of objects, total number of
                       galaxies (stars), total number of identified-magnitude
                       galaxies (stars), total number of identified magnitude
                       galaxies (stars) above S/N cut, percentage of FLAGS = 0
                       objects over total number of objects, number of galaxies 
                       per arcmin^2, number of stars per arcmin^2.

     - check plot    : plot showing the binned number count data and the
                       corresponding analytical fit, saved to folder_path.
    
    ---
    """

# Getting data and info from the header

    hdulist = pyfits.open(input_file,ignore_missing_end=True,memmap=True)

    if fmg.check_if_ldac(hdulist):

        hdudata = hdulist[2].data
    
        hdulist_string = hdulist[1].data[0][0]

        tile_ra = fmg.get_entry_ldac_header(hdulist_string,'CRVAL1')
        tile_dec = fmg.get_entry_ldac_header(hdulist_string,'CRVAL2')
        tile_seeing = fmg.get_entry_ldac_header(hdulist_string,'SEEING')
    
        tile_name = fmg.get_tile_name(input_file)
        tile_area = (21000*0.187/60)**2 # in arcmin^2, hardcoded for Megacam pixel scale, approximative value
    
    else:
        
        # Improve this part later for regular fits files

        hdudata = hdulist[1].data

        tile_ra = 1
        tile_dec = 1
        tile_seeing = 1
        tile_name = 'full'
        tile_area = 177*(21000*0.187/60)**2


# Performing cuts in different columns

    magdata = hdudata.field(field_names[0]).astype(np.float64)
    
    magerrdata = hdudata.field(field_names[1]).astype(np.float64)
 
    class_star = hdudata.field('CLASS_STAR').astype(np.float64)

    flags = hdudata.field('FLAGS').astype(np.float64)

     
    if gal_cut == True:
        stell_mask = (class_star < stell_idx)
 
    else:
        stell_mask = (class_star > stell_idx)

 
    if mag99 == False:
        mag99_mask = (magdata != 99)
 
    else:
        mag99_mask = True

    flag_mask = (flags == 0)

    flag_mask = stell_mask

    S2N_mask = magerrdata < 1.086/S2N_cut
     
    objs = []
 
    objs.append(len(magdata))
    objs.append(len(magdata[stell_mask]))
    objs.append(len(magdata[~stell_mask]))
    objs.append(len(magdata[stell_mask & mag99_mask]))
 
 
    data = magdata[stell_mask & mag99_mask & S2N_mask]

    objs.append(len(data))

    objs.append(len(magdata[flag_mask]))

 
    # Bin and cut the data
 
    binned_data = fmg.make_histogram(data, bin_size)

    binned_data = fmg.cut_unpopulated_bins(binned_data, 10) # Attention! Hard-coded value for the minimum number of counts for keeping bins 

    mag_sup = fmg.find_mag_sup(binned_data)

    cut_data = fmg.cut_mag_data(binned_data, mag_inf, mag_sup)


    # Perform fit

    a0 = 0; b0 = 0.3;  # Attention! Hard-coded initial parameter values for the fit

    fit_params = fmg.fit_mag_data(cut_data, a0, b0)


    # Get limit magnitude

    cut_inf_data = fmg.cut_mag_data(binned_data, mag_sup, 99)

    mag_lim = fmg.find_mag_lim(cut_inf_data, fit_params, bin_size)

    
    # Make the plots and save them to a given folder
    
    mkdir_p(folder_path + '/Plots/')

    fmg.make_plot_mag_pipe(binned_data, tile_name, folder_path, fit_params, mag_inf, mag_sup, mag_lim, bin_size, gal_cut, S2N_cut, stell_idx, plot_inf, plot_sup)

    
    # Return results of the fit
    
    return [tile_name, str(tile_ra), str(tile_dec), str(tile_seeing), str(bin_size), str(S2N_cut), str(stell_idx), str(fit_params[0]), str(fit_params[1]), str(mag_inf), str(mag_sup), str(mag_lim), str(objs[0]), str(objs[1]), str(objs[3]), str(objs[4]),str(objs[5]/objs[0]),str(objs[1]/tile_area),str(objs[2]/tile_area)]


# \cond
# ============================================================================

# Main function: Include ifs (or try/excepts to prevent crashing), and other minor changes


if __name__ == "__main__" :

    
    from optparse import OptionParser;

    usage="\n\n python %prog [flags] [options] 'input_file' 'mag_field_name' 'mag_error_field_name' \n\n This pipeline calculates the limiting magnitude for a given FITS catalogue. It has four mandatory arguments: the (input) FITS catalogue path and name, the magnitude field name to be used for the calculations (i.e. MAG_AUTO, MAG_APER, etc...) and the corresponding error field name (i.e. MAGERR_AUTO, MAGERR_APER, etc...). For a description of the procedure, check http://twiki.linea.gov.br/blabla\n\n # WARNING!: In its current state, this pipeline will only work with FITS_LDAC files generated by SExtractor (tested with v2.13.2)."
    
    parser = OptionParser(usage=usage);

    parser.add_option('-b','--bin_size', dest='bin_size', default=0.2,
                      help='magnitude bin size [default=0.2]');

    parser.add_option('-g','--galaxy_cut', dest='gal_cut', default='True',
                      help='perform analysis on galaxies (CLASS_STAR < stell_idx) if True, stars  (CLASS_STAR > stell_idx) if False [default=True]');

    parser.add_option('-i','--include_unid_mags', dest='mag99', default='False', 
                      help='include objects with unidentified magnitudes (MAG = 99) [default=False]');

    parser.add_option('-m','--mag_inf', dest='mag_inf', default=21, 
                      help='inferior magnitude value for which a fit of the object number count per magnitude bin should be performed [default=21]');

    parser.add_option('-n','--signal2noise_cut', dest='S2N_cut', default=5,
                      help='S/N cut value [default=5]');

    parser.add_option('--plot_x_inf', dest='plot_inf', default=16,
                      help='x-axis inferior limit of the plots [default=16]');

    parser.add_option('--plot_x_sup', dest='plot_sup', default=25,
                      help='x-axis superior limit of the plots [default=25]');

    parser.add_option('-s','--stellarity_index', dest='stell_idx', default=0.95,
                      help='stellarity index value [default=0.95]');


    (opts,args) = parser.parse_args();

    bin_size = float(opts.bin_size);

    gal_cut = (opts.gal_cut == 'True'); # Ugly, improve this comparison later
    mag99 = (opts.mag99 == 'True'); # Ugly, improve this comparison later
    mag_inf = float(opts.mag_inf);
    plot_inf = float(opts.plot_inf);
    plot_sup = float(opts.plot_sup);
    S2N_cut = float(opts.S2N_cut);
    stell_idx = float(opts.stell_idx);

    # How do I write a condition like this?

    '''
    if (len(sys.argv)!=4):
        print "\n\nThe number of arguments is wrong. Check documentation for usage instructions:\n\n";
        parser.print_help();
        print "";
        sys.exit(1);
    '''

# Get all the files in the given folder

    filenames = glob.glob(sys.argv[1])

    field_names = [sys.argv[2],sys.argv[3]]

    folder_path = os.path.dirname(filenames[0])


# Run the program for each file and send the results to the folder in question. The sending of plots is inside the pipeline, the sending of the final results for all files is here.

    results = []

    for i in range(len(filenames)):

        print i
       
        results.append(mag_pipeline(filenames[i], field_names, folder_path, mag_inf, bin_size, gal_cut, mag99, S2N_cut, stell_idx, plot_inf, plot_sup))

    results = np.array(results)


    np.savetxt(folder_path + '/fit_results_bins_'+ str(bin_size) + '_S2N_'+ str(int(S2N_cut))+ '_stell_'+str(stell_idx)+'.txt', results,fmt='%7s       %s        %s        %.4s      %s         %s       %.4s     %.8s       %.5s       %s        %s        %s           %s      %6s      %6s      %6s      %s      %s      %s')

    # Numpy 1.5 doesn't have header options for savetxt yet (to come in Numpy 2.0). Add a header using file commands
    
     
    file = open(folder_path + '/fit_results_bins_'+ str(bin_size) + '_S2N_'+ str(int(S2N_cut))+ '_stell_'+str(stell_idx)+'.txt', "r+")
     
    old = file.read() # read everything in the file
     
    file.seek(0) # rewind

    if gal_cut:
        file.write("ID_test       ra         dec         seeing    bin_size    S/N_cut    class_star     A_fit         b_fit     mag_inf     mag_sup    mag_comp    N_tot        N_gal      N_gal_99    N_gal_99_S/N      N_flags/N_tot    N_gal/arcmin^2       N_star/arcmin^2\n" + old)

    else:
        file.write("ID_test       ra         dec         seeing    bin_size    S/N_cut    class_star     A_fit         b_fit     mag_inf     mag_sup    mag_comp    N_tot        N_star      N_star_99    N_star_99_S/N     N_flags/N_tot    N_star/arcmin^2      N_star/arcmin^2\n" + old)
     
    file.close()
    

    sys.exit(0);



# ------------------
# \endcond
