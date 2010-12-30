""" Module to deal with FITS catalog read/access """

##@package fits_file

import sys;

import pyfits;
import string;
import numpy as np;
import re;

# ---

def select_rows(tbhdu, *indices):
    """ Read rows(indexes) from given HDU (catalog)
    
    Input:
     - tbhdu : pyfits.open('data.fit')[?]
        Table HDU, often "?" equals to 1
     - *indices : int,
        List of indexes to read from tbhdu
    
    Output:
     -> (new) BinTableHDU : sub-selection (rows) of tbhdu
    
    """
    
    _inds = [ i for i in indices ];
    cols = tbhdu.data.take(_inds);
    names = tbhdu.columns.names;
    formats = tbhdu.columns.formats;
    c = [];
    for i in range(len(names)):
        c.append(pyfits.Column(name=names[i],format=formats[i],array=list(cols[i])));
    
    return pyfits.new_table(pyfits.ColDefs(c));

# ---

def select_entries(tbhdu, field, *values):
    """ Read entries (lines) from given HDU (catalog)
    
    Input:
     - tbhdu : pyfits.open('data.fit')[?]
        Table HDU, often "?" equals to 1
     - field : str
        Field (column) name that 'values' should match
     - *values : <tbdata.field(field).dtype>
        List of values to match in 'field' entries
    
    Output:
     -> (new) BinTableHDU : sub-selection (rows) of tbhdu
    
    """
    
    _inds = [ tbhdu.data.field(field).tolist().index(_v) for _v in values ];

    return select_rows(_inds);

# ---

def select_columns(tbhdu, *fields): 
    """ Get a list of variables from a FITS catalog

    Input:
     - tbhdu : pyfits.open('data.fit')[?]
        Table HDU, often "?" equals to 1
     - cols : str,
        Comma separated list of variables to be read from 'hdulist'

    Output:
     -> (new) BinTableHDU, with just the selected fields
    
    """
	
    coldefs = tbhdu.columns;
    tbdata = tbhdu.data;
    inds = [ tbdata.names.index(id.upper()) for id in fields ];
    cols = [];
    for i in inds:
        cols.append( pyfits.Column(name=coldefs[i].name, format=coldefs[i].format, array=tbdata.field(i)) );
    coldefs = pyfits.ColDefs(cols);
    
    return pyfits.new_table(coldefs);

# ---
'''
def read_SE_data(segimg,tbdata,centroids):
    """ Read table data for requested objects

    Objects referenced on segmented image 'segimg'
    by 'centroids' positions are searched for
    in 'tbdata' (fits hdulist.data).
    
    Input:
     - segimg : <ndarray>
     - tbdata : <numpy record-array> (fits cat)
     - centroids : list of tuples [(x,y),]
    
    Output:
     -> {'*tbdata.names'},{'x','y'}
     1st dict: cross-checked objects entries/features from 'tbdata'
     2nd dict: non-checked object centroids
     
    """
    
    # Get the object IDs of selected objects (ds9)..
    #
    objIDs = segobjs.centroid2ID(segimg,centroids=centroids);
    
    # and see how many objects were not deteted.
    # Store these guys in somewhere, and remove them from objIDs list.
    Dnon = {};
    if (objIDs.count(0)):
        non_detected_centroids = [ centroids[i] for i in range(len(centroids)) if objIDs[i]==0 ];
        non_detected_indexes = [i for i in range(len(centroids)) if objIDs[i]==0 ];
    Dnon['x'],Dnon['y'] = zip(*non_detected_centroids);
    Dnon['ind'] = non_detected_indexes;
    
    objIDs = filter(lambda x: x, objIDs);
    
    # Read SE's catalog entries corresponding the ones in DS9's catalog..
    #
    params = tbdata.names;
    tbdata = fts.get_entries(tbdata, 'NUMBER', *objIDs);
    Ddata = fts.get_columns(tbdata, *params);
    
    return (Ddata,Dnon)
'''

# ---

def dict_to_tbHDU(dic, tbname='', *args):
    """ Return a table HDU with 'dic' keys|lists
    
    This function is designed to handle lists of number/string
    values (not arrays). The 'keys' used to label 'dic' values
    (list) are used to label the output data table and each
    column format (data type) will be estimated from data.

    It can be used one-dimensional numpy arrays in 'dic' entries.
    Valid data-types: 'int', 'int32', 'int64', 'float', 'float32',
    'float64', 'complex', 'long/object' and string(character)
    Maximum size of string values is 10 characters.

    Input:
     - dic : {'key':[],}
        Dictionary with a collection of key:column data
    
    Output:
     - <table HDU>
      A HDU associated to a table with 'dic' contents
    
    """
    
    if not len(args):
        args = dic.keys();

    c = [];
    for _arg in args:

        if type(dic[_arg]) == type(str()):
            continue;
        
        try:
            tipo = np.max(np.abs(dic[_arg])).dtype
            ndim = np.asarray(dic[_arg]).ndim
            if tipo == 'int64':
                tipo = str(ndim)+'K'
            elif tipo == 'int32' or tipo == 'int':
                tipo = str(ndim)+'J'
            elif tipo == 'float64':
                tipo = str(ndim)+'D'
            elif tipo == 'float32' or tipo == 'float':
                tipo = str(ndim)+'E'
            elif tipo == 'complex':
                tipo = str(ndim)+'C'
            elif tipo == 'bool':
                tipo = str(ndim)+'L'
            elif tipo == 'object':
                tipo = str(ndim)+'K'
            else:
                tipo = '10A'
        except:
            tipo = '20A'
            
        c.append(pyfits.Column(name=str(_arg).upper(),format=tipo,array=np.asarray(dic[_arg])));
    
    newtable = pyfits.new_table(pyfits.ColDefs(c));
    newtable.name = tbname;

    return newtable;
    
# ---

def dict_from_tbHDU(tbhdu, *fields):
    """ Get a list of variables from a FITS table HDU

    Input:
     - tbhdu : pyfits.open('data.fit')[1]
        Table HDU from a FITS catalog
     - fields : str,
        Comma separated list of variables to be read from 'hdulist'

    Output:
     -> {*fields} : a dictionary with given *args as keys and their values

    """
    
    if not len(fields):
        fields = tbhdu.data.names;

    dic = {}; 
    for arg in fields: 
        try:
            dic[arg] = tbdata.field(arg); 
        except KeyError:
            print "Variable %s does not exist in this catalog." % arg
    
    return (dic); 

# ---

def open_catalog_table( catalog_file ):
    """ Open given FITS catalog

    Input:
     - catalog_file : str
        FITS catalog filename

    Output:
     -> tbHDU : table HDU (1st FITS extension)
    """
     
    if ( string.find(catalog_file,'.fit') == -1 ):
        print >>sys.stderr,"Error: Does not seem to be a .fit(s) file."
        return (False); 

    return (pyfits.open(catalog_file)[1]);

# ---

# 
# Functions below are preserved for compatibility reasons:
#
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
    
    if ( len(args) == 0 ): 
        return (None); 

    dic = {}; 
    for arg in args: 
        try:
            dic[arg] = tbdata.field(arg); 
        except KeyError:
            print "Variable %s does not exist in this catalog." % arg
    
    return (dic); 


def open_fits_catalog( catalog_file ):
    """ Open given FITS catalog

    Input:
     - catalog_file : str
        FITS catalog filename

    Output:
     -> hdulist : FITS catalog "header" instance
    """
     
    if ( string.find(catalog_file,'.fit') == -1 ):
        print >>sys.stderr,"Error: Does not seem to be a .fit(s) filename."
        return (False); 

    return (pyfits.open(catalog_file)); 
