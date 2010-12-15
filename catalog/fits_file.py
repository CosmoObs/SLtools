""" Module to deal with FITS catalog read/access """

##@package fits_file

import sys;

import pyfits;
import string;

# ---
def get_entries(tbdata, ID, *indx):
    """ Read entries (lines) from given HDU (catalog)
    
    Input:
     - tbdata : Numpy Record-Array (hdulist[#].data)
     - ID   :   str
     - *indx   :   indices (#ID) to read
         
    Output:
     -> <Numpy Record-Array> : selected lines from 'tbdata'
        
    """
    
    _inds = [ tbdata.field(ID).tolist().index(i) for i in indx ];
    tabind = tbdata.take(_inds);

    return tabind;

# ---
def get_columns( tbdata, *args ): 
    """ Get a list of variables from a FITS catalog

    Input:
     - tbdata : Numpy Record-Array (hdulist[1].data)
        Header instance (pyfits.open()) of a FITS catalog
     - args : [str,]
        Comma separated list of variables to be read from 'hdulist'

    Output:
     -> {*args} : a dictionary with given *args as keys and their values

    """
	
    if ( len(args) == 0 ): 
        return (None); 

    dic = {}; 
    for arg in args: 
        try:
            dic[arg] = tbdata.field(arg); 
        except KeyError:
            print "Variable %s does not exist in this catalog." % arg
    
    return (dic); 

def get_fits_data(hdulist, *args):
    """ Get a list of variables from a FITS catalog

    Input:
     - hdulist : pyfits-header instance
        Header instance (pyfits.open()) of a FITS catalog
     - args : [str,]
        Comma separated list of variables to be read from 'hdulist'

    Output:
     -> {*args} : a dictionary with given *args as keys and their values

    """
    tbdata = hdulist[1].data;
    
    return get_columns(tbdata, *args);

# ---
def open_catalog( catalog_file ):
    """ Open given FITS catalog

    Input:
     - catalog_file : str
        FITS catalog filename

    Output:
     -> hdulist : header instance of given catalog

    """
     
    if ( string.find(catalog_file,'.fit') == -1 ):
        print >>sys.stderr,"Error: Does not seem to be a .fit(s) filename."
        return (False); 

    # Open fits table
    hdulist = pyfits.open(catalog_file)

    return (hdulist); 

open_fits_catalog=open_catalog;
# -
