""" Module to deal with FITS catalog read/access """

##@package fits_data

import sys;

import pyfits;
import string;
import numpy as np;
import re;

# ---

def merge_tbHDU(tbhdu_A, tbhdu_B):
    """
    Merge two tables (HDU) columns
    
    A new table (HDU) containing table 'A' and 'B'
    columns is output. It is expected that both tables
    lines are the same (same order). If one of them 
    has fewer lines, 0/NULL value is given to complete
    the lines/columns.
    
    Input:
     - tbhdu_A : pyfits.BinTableHDU
     - tbhdu_B : pyfits.BinTableHDU
    
    Output:
     - tbhdu : pyfits.BinTableHDU
        Result from input tables merge
    
    ---
    """
    
    new_tb = tbhdu_A.columns + tbhdu_B.columns;
    
    return pyfits.new_table(new_tb);

# ---

def extend_tbHDU(tbhdu_A, tbhdu_B):
    """
    Extend first tbHDU with second entries

    The output is a table (HDU) with column 'B' lines
    extending 'A' entries. Column names on both input
    table need to be the same.
    
    Input:
     - tbhdu_A : FITS binary table HDU
     - tbhdu_B : FITS binary table HDU
    
    Output:
     - tbhdu : FITS binary table HDU
        Result from extention, A+B
    
    ---
    """
    
    Nrows_A = tbhdu_A.header['NAXIS2'];
    Nrows_B = tbhdu_B.header['NAXIS2'];
    Nrows = Nrows_A + Nrows_B;
    
    new_tb = pyfits.new_table(tbhdu_A.columns,nrows=Nrows);
    
    for name in tbhdu_A.columns.names:
        new_tb.data.field(name)[Nrows_A:] = tbhdu_B.data.field(name);
    
    return new_tb;
    
# ---

def select_columns(tbhdu, *fieldnames): 
    """
    Select particular columns from given table
    
    A new table with only the asked columns ('fieldnames')
    is output.

    Input:
     - tbhdu : pyfits.open('data.fit')[?]
        Table HDU, often "?" equals 1
     - cols : str,
        Comma separated list of variables to be read from 'hdulist'

    Output:
     -> (new) BinTableHDU, with just the selected fields
    
    ---
    """
	
    coldefs = tbhdu.columns;
    tbdata = tbhdu.data;
    inds = [ tbdata.names.index(id.upper()) for id in fieldnames ];
    cols = [];
    for i in inds:
        cols.append( pyfits.Column(name=coldefs[i].name, format=coldefs[i].format, array=tbdata.field(i)) );
    coldefs = pyfits.ColDefs(cols);
    
    return pyfits.new_table(coldefs);

# ---

def select_rows(tbhdu, *indices):
    """
    Read rows(indexes) from given HDU (catalog)
    
    A new table with only the asked table indexes,
    'indices', is output.
    
    Input:
     - tbhdu     [BinTableHDU] : FITS table HDU
     - *indices         [int,] : List of indexes to read from tbhdu
    
    Output:
     -> (new) BinTableHDU : sub-selection (rows) of tbhdu

    ---
    """
    
    _inds = [ i for i in indices ];
    data = tbhdu.data.take(_inds);

    return pyfits.BinTableHDU(data);
    
# ---

def select_entries(tbhdu, fieldname, *values):
    """
    Read entries (lines) from given HDU (catalog)
    
    'values' matching entries in column 'fieldname'
    is used to select rows from given 'tbhdu'. A new
    table (HDU) is generated from selected lines.
    
    Input:
     - tbhdu     [BinTableHDU] : Table HDU, often "?" equals to 1
     - field            [str]  : Field (column) name that 'values' should match
     - *values  [field.dtype]  : List of values to match in 'field' entries
    
    Output:
     -> (new) BinTableHDU : sub-selection (rows) of tbhdu
    
    ---
    """
    
    _inds = [ tbhdu.data.field(fieldname).tolist().index(_v) for _v in values ];

    return select_rows(tbhdu,*_inds);

# ---

def sample_entries(tbhdu, **kwargs):
    """
    Retrieve a slice of given table HDU
    
    'kwargs' are keyworded arguments to pass threshold
    parameter values for the filtering. kwargs' keys are
    expected to be fieldnames in tbhdu, and the values
    can be just a single number - used as minimum threshold - 
    or a tuple - for (min,max) thresholds.
    
    Input:
     - tbhdu : pyfits.BinTableHDU
        FITS table HDU
     - kwargs : {'key1'=value1,}
        Key=Value, a dictionary

    Output:
     - tbhdu : pyfits.BinTableHDU
        Binary table HDU with selected lines
    
    ---
    """

    tbsamp = tbhdu.copy();

    for _key in kwargs.keys():
    
        try:
            _min,_max = kwargs[_key];
            _in = np.where(tbsamp.data.field(_key) >= _min)[0].tolist();
            _ax = np.where(tbsamp.data.field(_key) <= _max)[0].tolist();
            inds = list(set.intersection(set(_in),set(_ax)));
            tbsamp = select_rows(tbsamp,*inds);
            
        except:
            _min = kwargs[_key];
            inds = np.where(tbsamp.data.field(_key) >= _min)[0].tolist();
            tbsamp = select_rows(tbsamp,*inds);
            
    return tbsamp;

# ---

def dict_to_tbHDU(dic, tbname='', *fieldnames):
    """ Return a table HDU with 'dic' keys|lists
    
    This function is designed to handle lists of number/string
    values (not arrays). The 'keys' used to label 'dic' values
    (list) are used to label the output data table and each
    column format (data type) will be estimated from data.

    It can be used one-dimensional numpy arrays in 'dic' entries.
    Valid data-types: 'int', 'int32', 'int64', 'float', 'float32',
    'float64', 'complex', 'long/object' and string(character)
    Maximum size of string values is 10 characters.
    
    If 'fieldnames' are given, only them (from 'dic') are read
    and written to output (new) table.

    Input:
     - dic : {'key':[],}
        Dictionary with a collection of key:column data
     - tbname : str
        Output tbale header name
     - fieldnames : str,
        [optional] keys to read from 'dic'
    
    Output:
     - tbhdu : pyfits.BinTableHDU
      A HDU associated to a table with 'dic' contents
    
    ---
    """
    
    args = fieldnames;
    
    if not len(args):
        args = dic.keys();

    # Function verify datatypes for convertion
    #
    def is_numeric(lit):
        try:
            return int(lit);
        except ValueError:
            pass;
        try:
            return float(lit);
        except ValueError:
            pass;
        try:
            return complex(lit);
        except:
            pass;
        return lit;
        
    c = [];
    to_header = {};
    for _arg in args:

        if type(dic[_arg]) == type(str()):
            to_header[_arg] = dic[_arg];
            continue;
        
        valores = [ is_numeric(s_i) for s_i in dic[_arg] ];
        
        try:
            tipo = np.max(np.abs(valores)).dtype
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
                tipo = '20A'
        except:
            tipo = '40A'
            
        c.append(pyfits.Column(name=str(_arg).upper(),format=tipo,array=np.asarray(dic[_arg])));

    newtable = pyfits.new_table(pyfits.ColDefs(c));
    newtable.name = tbname;
    
    for _k in to_header.keys():
        newtable.header.update('hierarch'+str(_k), to_header[_k], 'Entry automatically added(not a list)');

    return newtable;
    
# ---

def dict_from_tbHDU(tbhdu, *fieldnames):
    """ Get a list of variables from a FITS table HDU

    Input:
     - tbhdu : pyfits.open('data.fit')[1]
        Table HDU from a FITS catalog
     - fieldnames : str,
        Comma separated list of variables to be read from 'hdulist'

    Output:
     -> {*fieldnames} : a dictionary with given *args as keys and their values

    """
    
    if not len(fieldnames):
        fieldnames = tbhdu.data.names;

    dic = {}; 
    for arg in fieldnames: 
        try:
            dic[arg] = tbhdu.data.field(arg); 
        except KeyError:
            print "Variable %s does not exist in this catalog." % arg
    
    return (dic); 



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
