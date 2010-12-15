import sys

""" Module to deal with ascii (CSV) catalog files """

##@package text_file

"""
The following functions are meant to help in reading and
writing text/csv catalogs given fieldnames and necessary structure.
"""

# ---
def dict2csv(dictionary, fieldnames, filename, mode='w', delimiter=',', quotechar='"'):
    """ Write a CSV catalog from dictionary contents
    
    Input:
     - filename : str
        Name of csv catalog to write
     - filednames : [str,]
        Fieldnames to read from 'dictionary'
     - dictionary : {str,}
        Contents to be write in csv catalog
     - mode : str
        Write a new catalog, 'w', or append to an existing one, 'a'.
    
    Output:
     <bool>
    
    """
    import csv;
    
    if not (list(set(fieldnames)-set(dictionary.keys())) == []):
        print "Error: Given dictionary does not contain every requested fieldnames.";
        return (False);

    # Initialize csv catalog
    catFile = open(filename,mode);
    catObj = csv.writer(catFile, delimiter, quotechar);
    catObj.writerow(fieldnames);
    
    LL = [ dictionary[_k] for _k in fieldnames ];
    for _row in zip(*LL):
        catObj.writerow(_row);
    catFile.close();

    return (True);

# ---
def dict_from_csv(filename, fields, remove_n_lines=1, delimiter=',', quotechar='"'):
    """ Read CSV catalog and return a dictionary with the contents
    
    It is assumed that the first line to be read is the catalog header,
    where the column names are defined('fieldnames'), that's why
    'remove_n_lines=1' by default.
    If necessary to remove a certain number of lines from the beginning
    of file ('filename') till the data line, use 'remove_n_lines' to
    give the number of "invalid" lines from the head of the catalog.
    
    Input:
     - filename : str
        Name of csv catalog to read
     - fieldnames : [str,]
        Fieldnames to be read from catalog
     - remove_n_lines : int
        Number of lines to remove from the head of 'filename'
     - delimiter : str
        Delimiter to use between columns in 'filename'
     - quotechar : str
        Char used for limiting fields in 'filename'
    
    Output:
     -> {*fieldnames}
     
    """
    import csv;

    # Initialize output dictionary
    Dout = {};
    for k in fieldnames:
        Dout[k] = [];
    
    #Initialize csv reader
    catFile = open(filename,'r');

    lixo_head = [ catFile.next() for i in range(remove_n_lines) ];

    catObj = csv.DictReader(catFile,fieldnames,delimiter,quotechar);
    for row in catObj:
        for k in fieldnames:
            Dout[k].append(row[k]);

    return Dout;
