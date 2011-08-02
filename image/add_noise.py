""" add_noise """

##@package add_noise

"""
Add poisson noise to a list of FITS images
"""

import pyfits
import string
import numpy as np
from numpy.random import poisson

def add_poisson_noise(list_fitsFiles):

    """ Add poissonian noise to given list of FITS files.

    Input:
    -> list of fits filenames
	
    Output:
    <- list of fits filenames with poisson noise added

    """

    list_noisyFiles = []
    for Arc_filename in list_fitsFiles:

        # Unfortunately the convolution process outputs to a file..
        # Lets open it to an array..
        #
        ArcImg_array,ArcImg_header = pyfits.getdata( Arc_filename, header=True );

        ###########################################################################
        # Here I am treating a strange "behavior" (negative counts)
        #
        Arc_tmp = np.zeros( ArcImg_array.shape, dtype=float );
        non_null = np.where( ArcImg_array > 0.0 );
        Arc_tmp[ non_null ] = ArcImg_array[ non_null ];
        ArcImg_array = Arc_tmp.copy();
        del Arc_tmp;
        ###########################################################################


        # Re-sort pixel intensities of arc image in a poissonian fashion
        #    >>> arco = poisson(arco)
        #
        ArcImg_array = poisson( ArcImg_array );


        # Outputs Image to fits files
        #
        Filename_OUT = ''
        for i in range( len(Arc_filename) - 6): # the last 6 characters corresponds to 'color.fits'
            Filename_OUT = Filename_OUT + Arc_filename[i]
        Filename_OUT = Filename_OUT + 'ns_' + Arc_filename[-6] + '.fits' # Arc_filename[-6] is the 6th str from the end (the colour)


#        _part = string.split( Arc_filename, sep="_");
#        root = string.join( _part[:-1], sep="_");
#        ext = _part[-1];
#        Filename_OUT = root+'_ns_'+ext;


        pyfits.writeto( Filename_OUT, ArcImg_array, ArcImg_header );
        
        list_noisyFiles.append( Filename_OUT );

    return (list_noisyFiles);
