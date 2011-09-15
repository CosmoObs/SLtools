from __future__ import division
import sys
import os
import numpy as np
import math
import cmath
from scipy import optimize
import collections
import fit_mag as fmg
import glob
import pyfits

def mag_pipeline(input_file, field_names, mag_inf, bin_size, gal_cut, mag99, S2N_cut, stell_idx):

    # Get the data from the file

    hdu = pyfits.open(input_file,ignore_missing_end=True)

    magdata = hdu[2].data.field(field_names[0])
    magerrdata = hdu[2].data.field(field_names[1])

    class_star = hdu[2].data.field('CLASS_STAR')

    # Get the file number

    image_name = input_file.replace('/home/brunomor/Documents/Trabalho/Repositorio/sltools/catalog/cs82/S82p28m_y.V2.3A.swarp.cut_', '')
    image_number = float(image_name.replace('.fit', ''))
    
    

    # Perform the cuts on the data
    
    if gal_cut == True:
         stell_mask = (class_star < stell_idx)

    else:
         stell_mask = (class_star > stell_idx)


    if mag99 == False:
        mag99_mask = (magdata != 99)

    else:
        mag99_mask = (magdata != -1000)

    # S2N_mask = magerrdata < 1.086/S2N_cut
    

    data = magdata[stell_mask & mag99_mask] # & S2N_mask]


    # Bin and cut the data

    binned_data = fmg.bin_mag_data(data, bin_size)

    mag_sup = fmg.find_mag_sup(binned_data)

    cut_data = fmg.cut_mag_data(binned_data, mag_inf, mag_sup)

    # Perform fit and get limit magnitude 

    a0 = 0;

    b0 = 0.3;

    p = fmg.fit_mag_data(cut_data, a0, b0)

    fit_params = list(p)

    cut_inf_data = fmg.cut_mag_data(binned_data,mag_inf,99)

    mag_lim = fmg.find_mag_lim(cut_inf_data,fit_params)

    # Make the plots and save them to a given folder. (still to do)

    # Return results of the fit

    return [image_number, fit_params[0], fit_params[1], mag_inf, mag_sup, mag_lim]


# ==========================

if __name__ == "__main__" :

    '''
    from optparse import OptionParser;

    usage="\n  %prog [flags] [options] ??? ??? "
    
    parser = OptionParser(usage=usage);

    parser.add_option('-b','--bin_size', dest='bin_size', default=0.5, help='bla');

    # parser.add_option('-f','--flag_cut', dest='flag_cut', default=False, help='bla');

    parser.add_option('-g','--galaxy_cut', dest='gal_cut', default=True, help='bla');

    parser.add_option('-i','--include_unidentified_magnitudes', dest='mag99', default=False, help='bla'); #Check this definition later

    parser.add_option('-n','--signal2noise_cut', dest='S2N_cut', default=5, help='bla');

    parser.add_option('-s','--stellarity_index', dest='stell_idx', default=0.8, help='bla');


    (opts,args) = parser.parse_args();

    bin_size = opts.bin_size;
    gal_cut = opts.galaxy_cut;
    mag99 = opts.mag99;
    S2N_cut = opts.S2N_cut;
    stell_idx = opts.stell_idx;

    # if m_inf is not given, take the minimum for which N = 10.

    if ( opts=={} ):
        print"";
        parser.print_help();
        print "";
        sys.exit(1);

    if ( input_file==None or ra_0 == None or dec_0 == None ):
        print "";
        parser.print_help();
        print "";
        sys.exit(1);

# Include here a series of tests for the optional parameters to avoid inconsistencies.

# Get all the files in the given folder.
'''
    bin_size = 0.5;
    gal_cut = True;
    mag99 = False;
    S2N_cut = 5;
    stell_idx = 0.8;

    filenames = glob.glob(sys.argv[1])

    field_names = [sys.argv[2],sys.argv[3]]

    mag_inf = float(sys.argv[4])

# Run the program for each file and send the results to the folder in question. The sending of plots is inside the pipeline, the sending of the final results for all files is here.

    results = []

    for i in range(len(filenames)):
    
        print mag_pipeline(filenames[i], field_names, mag_inf, bin_size, gal_cut, mag99, S2N_cut, stell_idx)
    '''    
        results.append(mag_pipeline(filenames[i], field_names, mag_inf, bin_size, gal_cut, mag99, S2N_cut, stell_idx))
 
    results = np.array(results)

    print results
    '''
    # The savetxt command must be here.


    sys.exit(0);
