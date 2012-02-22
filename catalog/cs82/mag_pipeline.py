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

sys.path.append("/home/brunomor/Documents/Trabalho/Repositorio") ### Append Bruno
from sltools.catalog import fits_data as fd
from sltools.coordinate import wcs_conversion as wcscv


def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST:
            pass
        else: raise


def mag_pipeline(input_file, field_names, folder_path, mag_inf, bin_size, gal_cut, mag99, S2N_cut, stell_idx, plot_inf, plot_sup):


#######################################################

    hdudata = pyfits.open(input_file,ignore_missing_end=True)[2].data
    hdulist_string = pyfits.open(input_file,ignore_missing_end=True)[1].data[0][0]

    tile_ra = fmg.get_entry_tile_cs82_ldac(hdulist_string,'CRVAL1')
    tile_dec = fmg.get_entry_tile_cs82_ldac(hdulist_string,'CRVAL2')
    tile_seeing = fmg.get_entry_tile_cs82_ldac(hdulist_string,'SEEING')
    
    tile_name = fmg.get_tile_name(input_file)
    tile_area = (21000*0.187/60)**2 # in arcmin^2
 
    magdata = hdudata.field(field_names[0])
    magdata = magdata.astype(np.float64)
    
    magerrdata = hdudata.field(field_names[1])
    magerrdata = magerrdata.astype(np.float64)
 
    class_star = hdudata.field('CLASS_STAR')
    class_star = class_star.astype(np.float64)

    flags = hdudata.field('FLAGS')

     
    if gal_cut == True:
        stell_mask = (class_star < stell_idx)
 
    else:
        stell_mask = (class_star > stell_idx)

 
    if mag99 == False:
        mag99_mask = (magdata != 99)
 
    else:
        mag99_mask = True

    flag_mask = (flags == 0) 

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

 
    binned_data = fmg.bin_mag_data(data, bin_size)

#######################################################

# Versao lenta demais com o sample_entries

    '''
    # Get the data from the file and the tile name


    tile_name = fmg.get_tile_name(input_file)


    hdu = pyfits.open(input_file,ignore_missing_end=True)[2]

    # Count objects in different cuts
    
    objs = []

    objs.append(len(hdu.data))

   
    # Perform the cuts on the data


    # Stellarity cut
    
    if gal_cut == True:
         hdu = fd.sample_entries(hdu,CLASS_STAR=(0,stell_idx))

    else:
         hdu = fd.sample_entries(hdu,CLASS_STAR=(stell_idx,1))

    objs.append(len(hdu.data))


    # Mag != 99 cut

    if mag99 == False:
        hdu = fd.sample_entries(hdu,MAG_AUTO=(0,98))

    objs.append(len(hdu.data))

    
    # S/N cut

    hdu = fd.sample_entries(hdu, MAGERR_AUTO=(0,1.086/S2N_cut))

    objs.append(len(hdu.data))
    

    # Bin the data and get rid of undesired bins
    
    magdata = hdu.data.field(field_names[0])
    magdata = magdata.astype(np.float64)

    binned_data = fmg.bin_mag_data(magdata, bin_size)
    '''
####################################################################

    binned_data = fmg.cut_unpopulated_bins(binned_data, 10) # Attention! Hard-coded value for the minimum number of counts for keeping bins 

    mag_sup = fmg.find_mag_sup(binned_data)

    cut_data = fmg.cut_mag_data(binned_data, mag_inf, mag_sup)



    # Perform fit

    a0 = 0;  # Attention! Hard-coded initial parameter values for the fit

    b0 = 0.3;

    p = fmg.fit_mag_data(cut_data, a0, b0)

    fit_params = list(p)


    # Get limit magnitude

    cut_inf_data = fmg.cut_mag_data(binned_data, mag_sup, 99)

    mag_lim = fmg.find_mag_lim(cut_inf_data, fit_params, bin_size)

    
    # Make the plots and save them to a given folder
    
    mkdir_p(folder_path + '/Plots/')

    fmg.make_plot_mag_pipe(binned_data, tile_name, folder_path, fit_params, mag_inf, mag_sup, mag_lim, bin_size, gal_cut, S2N_cut, stell_idx, plot_inf, plot_sup)
    
    # Return results of the fit
    
    hdulist_string = 0
    hdudata = 0
    magdata = 0
    magerrdata = 0
    classstar = 0
    data = 0
    cut_inf_data = 0
    cut_data = 0
    binned_data = 0
    cut_inf_data = 0
    stell_mask = 0
    mag99_mask = 0
    S2N_mask = 0
    flags = 0
    flag_mask = 0
    
    del hdudata
    del hdulist_string


    return [tile_name, str(tile_ra), str(tile_dec), str(tile_seeing), str(bin_size), str(S2N_cut), str(stell_idx), str(fit_params[0]), str(fit_params[1]), str(mag_inf), str(mag_sup), str(mag_lim), str(objs[0]), str(objs[1]), str(objs[3]), str(objs[4]),str(objs[5]/objs[0]),str(objs[1]/tile_area),str(objs[2]/tile_area)]


# ==========================

# Main function: Improve documentation of usage and optional arguments, include ifs (or try/excepts to prevent crashing), and other minor changes


if __name__ == "__main__" :

    
    from optparse import OptionParser;

    usage="\n  %prog [flags] [options] ??? ??? "
    
    parser = OptionParser(usage=usage);

    parser.add_option('-b','--bin_size', dest='bin_size', default=0.5, help='bla');

    # parser.add_option('-f','--flag_cut', dest='flag_cut', default=False, help='bla'); # Add these cuts later

    parser.add_option('-g','--galaxy_cut', dest='gal_cut', default='True', help='bla');

    parser.add_option('-i','--include_unidentified_magnitudes', dest='mag99', default='False', help='bla');

    parser.add_option('-n','--signal2noise_cut', dest='S2N_cut', default=5, help='bla');

    parser.add_option('-s','--stellarity_index', dest='stell_idx', default=0.8, help='bla');


    (opts,args) = parser.parse_args();

    bin_size = float(opts.bin_size);

    gal_cut = (opts.gal_cut == 'True'); # Ugly, improve this comparison later
    mag99 = (opts.mag99 == 'True'); # Ugly, improve this comparison later
    S2N_cut = float(opts.S2N_cut);
    stell_idx = float(opts.stell_idx);
    

# Get all the files in the given folder

    filenames = glob.glob(sys.argv[1])

    field_names = [sys.argv[2],sys.argv[3]]

    '''
    if gal_cut:
        folder_path = filenames[0][0:filenames[0].find('S82')] + 'Results/' + field_names[0] + '_GAL/'

    else:
        folder_path = filenames[0][0:filenames[0].find('S82')] + 'Results/' + field_names[0] + '_STAR/'
    '''

    folder_path = os.path.dirname(filenames[0])

    mag_inf = float(sys.argv[4])

    plot_inf = float(sys.argv[5])
    plot_sup = float(sys.argv[6])


# Run the program for each file and send the results to the folder in question. The sending of plots is inside the pipeline, the sending of the final results for all files is here.

    results = []

    for i in range(len(filenames)):

        print i
       
        results.append(mag_pipeline(filenames[i], field_names, folder_path, mag_inf, bin_size, gal_cut, mag99, S2N_cut, stell_idx, plot_inf, plot_sup))

    results = np.array(results)


    np.savetxt(folder_path + '/fit_results_bins_'+ str(bin_size) + '_S2N_'+ str(int(S2N_cut))+ '_stell_'+str(stell_idx)+'.txt', results,fmt='%7s       %s        %s        %.4s      %s         %s       %.4s     %.8s       %.5s       %s        %s        %s           %s      %6s      %6s      %6s      %s      %s      %s')

    # np.savetxt(folder_path + '/fit_results_bins_'+ str(bin_size) + '_S2N_'+ str(int(S2N_cut))+ '_stell_'+str(stell_idx)+'.txt',results,fmt='%7s       %s       %s        %.4s      %s         %s       %.4s     %.8s       %.5s       %s        %s        %s           %s      %6s      %6s      %6s      %s      %s      %s')

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
