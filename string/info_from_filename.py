#!/usr/bin/env python
# ====================================================================
# Authors:
# Maria Elidaiana - mariaeli@cbpf.br
# Current Fev/2014
# ====================================================================


import re
import os

def get_catalog_name(pathtofile):
    """
    Function to get the catalog name from a full path. 

    This function receives a string with the entire path to a catalog 
    and return only the name of catalog. If there is no catalog or the
    path do not set for a valid catalog (.fits, .fit or .cat) an error
    message is displayed and the code stop.     


    Input:
     - pathtofile     str : The path or name to the catalog 

    Output:
    - A string with the catalog name

    ---
    """


    filename = os.path.basename(pathtofile)

    if re.search(".fits$",filename.lower()) or re.search(".fit$",filename.lower()):
        pass  
    elif re.search(".cat$",filename.lower()):    
        pass  
    else:
        print "Catalog not found! Try again."
        exit()

    return filename

def get_survey_info(filename):
    """
    Function to find out the survey infos through the name of a catalog. 

    This function receives a string with the name of a catalog and find
    out on which survey it belongs. After that, it extracts tile number,
    band name and, if it exists, the profile fitted. These infos are 
    returned in a list in this order: [tile, band, model]. 

    Input:
     - filename   str : Catalog name 

    Output:
    - list with the name of tile, band and model, respectively

    ---
    """


    cs82pattern = re.compile(r'S82(\w\d*\w)_(\w*).(\w*.\w*).swarp.cut.(\w*).(\w*)')
    obj1 = re.search(cs82pattern, filename)
    if obj1:
        surveyname = 'CS82'
    elif obj1==None:
        pass

    vics82pattern = re.compile(r'S82(\w\d*\w)_(\w*)_(\w*).(\w*)$')
    obj2 = re.search(vics82pattern, filename)
    if obj2:
        surveyname = 'VICS82'
    elif obj2==None:
        pass


    if surveyname.lower()=='cs82': 
        tilenumber = cs82pattern.search(filename).groups()[0]
        bandname = cs82pattern.search(filename).groups()[1] 
        modelname = cs82pattern.search(filename).groups()[3]
    elif surveyname.lower()=='vics82':
        tilenumber = vics82pattern.search(filename).groups()[0]
        bandname = vics82pattern.search(filename).groups()[1] 
        modelname = ''
    return [tilenumber, bandname, modelname]









