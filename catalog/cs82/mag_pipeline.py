# A first draft of the complete pipeline.

# Import necessary packages.

def mag_pipeline(image_name(with_path), bin_size, gal_cut, mag99, mag_inf, S2N_cut, stell_idx):

    # Get the data from the file

    NUMBER, ALPHA_J2000, DELTA_J2000, X_IMAGE, Y_IMAGE, MAG_ISO, MAGERR_ISO, MAG_ISOCOR, MAGERR_ISOCOR, MAG_APER, MAGERR_APER, MAG_AUTO, MAGERR_AUTO, MAG_BEST, MAGERR_BEST, FLUX_RADIUS, CLASS_STAR, MU_MAX, FLAGS = np.loadtxt(filename, usecols=(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14 ,15, 16, 17, 18), unpack=True)

    data = np.array(data)

    #cut_mag99
   
    #cut_stellarity

    #cut_signal2noise

    binned_data = bin_mag_data(data, bin_size, m_inf)

    mag_sup = find_mag_sup(binned_data)

    cut_data = cut_mag_data(binned_data, m_inf, mag_sup)

    p = fit_mag_data(cut_data, a0, b0)

    fit_params = list(p)

    cut_inf_data = cut_mag_data(binned_data,m_inf,99)

    find_mag_lim(cut_inf_data,fit_params)

# Make the plots and save them to a given folder.


# Create the main function: include the necessary options (different cuts, inputs, etc with parser blablabla).

if __name__ == "__main__" :

    from optparse import OptionParser;

    parser = OptionParser();

    parser.add_option('-b','--bin_size', dest='bin_size', default=0.5, help='bla');

    # parser.add_option('-f','--flag_cut', dest='flag_cut', default=False, help='bla');

    parser.add_option('-g','--galaxy_cut', dest='gal_cut', default=True, help='bla');

    parser.add_option('-i','--include_unidentified_magnitudes', dest='mag99', default=False, help='bla'); #Check this definition later

    parser.add_option('-m','--mag_inf_cut', dest='m_inf', default=19, help='bla'); #Input obrigatório

    parser.add_option('-n','--signal2noise_cut', dest='S2N_cut', default=5, help='bla');

    parser.add_option('-s','--stellarity_index', dest='stell_idx', default=0.8, help='bla');


    (opts,args) = parser.parse_args();

    bin_size = opts.bin_size;
    gal_cut = opts.galaxy_cut;
    mag99 = opts.mag99;
    m_inf = opts.m_inf;
    S2N_cut = opts.S2N_cut;
    stell_idx = opts.stell_idx;

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

    filenames = glob.glob('*.cat')

# Run the program for each file and send the results to the folder in question. The sending of plots is inside the pipeline, the sending of the final results for all files is here.

    for i in range(len(filenames)):
        
        result = mag_analysis(path_name_filenames[i],...)

        "savetxt(result) bla bla to the necessary folder"

        return


    sys.exit(0);
