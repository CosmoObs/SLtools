'''@file 
	get_fits_data module
'''


'''@package get_fits_data 
	Get a list of variables from a FITS catalog 
'''

import pyfits 
def get_data( hdulist, *args ): 

    ''' Get a list of variables from a FITS catalog

    @param hdulist of the catalog
    @param args is a comma separated list of variables to be read 

    @return a dictionary with the variable names and values
    '''
	
    if ( len(args) == 0 ): 
        return (None); 

    tbdata = hdulist[1].data; 

    dic = {}; 
    for arg in args: 
        try:
            dic[arg] = tbdata.field(arg); 
        except KeyError:
            print "Variable %s does not exist in this catalog." % arg
    
    return (dic); 

