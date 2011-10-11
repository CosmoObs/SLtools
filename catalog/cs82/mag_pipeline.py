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


def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST:
            pass
        else: raise


def mag_pipeline(input_file, field_names, folder_path, mag_inf, bin_size, gal_cut, mag99, S2N_cut, stell_idx, plot_inf, plot_sup):

    # Get the data from the file and the tile name


    tile_name = fmg.get_tile_name(input_file)


    hdu = pyfits.open(input_file,ignore_missing_end=True)[2]

   
    # Perform the cuts on the data


    # Stellarity cut
    
    if gal_cut == True:
         hdu_stell_cut = fd.sample_entries(hdu,CLASS_STAR=(0,stell_idx))

    else:
         hdu_stell_cut = fd.sample_entries(hdu,CLASS_STAR=(stell_idx,1))


    # Mag != 99 cut

    if mag99 == False:
        hdu_stell_mag99_cut = fd.sample_entries(hdu_stell_cut,MAG_AUTO=(0,98))

    else:
        hdu_stell_mag99_cut = hdu_stell_cut

    
    # S/N cut

    hdu_all_cuts = fd.sample_entries(hdu_stell_mag99_cut, MAGERR_AUTO=(0,1.086/S2N_cut))
    

    # Count objects in different cuts
    
    objs = [len(hdu.data),len(hdu_stell_cut.data),len(hdu_stell_mag99_cut.data),len(hdu_all_cuts.data)]



    # Bin the data and get rid of undesired bins
    
    magdata = hdu_all_cuts.data.field(field_names[0])
    magdata = magdata.astype(np.float64)

    binned_data = fmg.bin_mag_data(magdata, bin_size)

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

    # Return results of the fit


    return [tile_name, str(bin_size), str(S2N_cut), str(fit_params[0]), str(fit_params[1]), str(mag_inf), str(mag_sup), str(mag_lim), str(objs[0]), str(objs[1]), str(objs[2]), str(objs[3])]


    # Make the plots and save them to a given folder
'''
    fmg.make_plot(binned_data, field_names, folder_path, fit_params, mag_inf, mag_sup, mag_lim, bin_size, gal_cut, S2N_cut, image_number, plot_inf, plot_sup)

'''

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
       
        results.append(mag_pipeline(filenames[i], field_names, folder_path, mag_inf, bin_size, gal_cut, mag99, S2N_cut, stell_idx, plot_inf, plot_sup))

    results = np.array(results)


    np.savetxt(folder_path + '/fit_results_bins_'+ str(bin_size) + '_S2N_'+ str(int(S2N_cut))+'.txt', results,fmt='%s       %s         %s       %.8s       %.5s       %s        %s        %s      %s      %s      %s      %s')

    # np.savetxt(folder_path + '/fit_results_bins_'+ str(bin_size) + '_S2N_'+ str(int(S2N_cut))+'.txt',results,fmt='%c          %-3.2f         %.0f       %f    %.4f     %.2f      %.2f      %.2f      %-6.0f      %-6.0f        %-6.0f        %-6.0f')

    # Numpy 1.5 doesn't have header options for savetxt yet (to come in Numpy 2.0). Add a header using file commands
    
    S''' 
    file = open(folder_path + '/fit_results_bins_'+ str(bin_size) + '_S2N_'+ str(int(S2N_cut))+'.txt', "r+")
     
    old = file.read() # read everything in the file
     
    file.seek(0) # rewind

    if gal_cut:
        file.write("ID_test    bin_size    S/N_cut        A_fit         b_fit     mag_inf     mag_sup    mag_comp    N_tot        N_gal      N_gal_99    N_gal_99_S/N \n" + old)

    else:
        file.write("ID_test    bin_size    S/N_cut        A_fit         b_fit     mag_inf     mag_sup    mag_comp    N_tot        N_star      N_star_99    N_star_99_S/N \n" + old)
     
    file.close()
    '''

    sys.exit(0);
