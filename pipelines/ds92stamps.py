#!/usr/bin/python
import sys;

""" Module for automate the process of object images segmentation """

##@package ds92stamps

""" This module executes a set of functions (some of them defined here)
with the goal to have a set of features and poststamps of pre-selected
objects from a given image.

"""

import os;
import pyfits;
import numpy as np;
from numpy import asarray as ndaa;
import re;

from sltools.io import check_dependencies as ckdep;
from sltools.catalog import ds9region as ds9reg;
from sltools.image import sextractor as se;
from sltools.image import segobjs;
from sltools.image import imcp;
from sltools.catalog import fits_file as fts;


# GLOBAL variables:
#
INSTRUMENT='HST';
PARAMS=['NUMBER','X_IMAGE','Y_IMAGE','X2_IMAGE','Y2_IMAGE','XY_IMAGE','FWHM_IMAGE','ELONGATION','ELLIPTICITY','ALPHA_J2000','DELTA_J2000'];


# ---

def read_SEcat_data(segimg,tbdata,centroids):
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
    Dnon['x'],Dnon['y'] = zip(*non_detected_centroids);
    
    objIDs = filter(lambda x: x, objIDs);
    
    # Read SE's catalog entries corresponding the ones in DS9's catalog..
    #
    params = tbdata.names;
    tbdata = fts.get_entries(tbdata, 'NUMBER', *objIDs);
    Ddata = fts.get_columns(tbdata, *params);
    
    return (Ddata,Dnon)

# ---

def write_csvcat(filename, fieldnames, dictionary, mode='w'):
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
        Default is 'w'.
    
    Output:
     <bool>
    
    """
    import csv;
    
    if not (list(set(fieldnames)-set(dictionary.keys())) == []):
        print "Error: Given dictionary does not contain every requested fieldnames.";
        return (False);

    # Initialize csv catalog
    catFile = open(filename,mode);
    catObj = csv.writer(catFile, delimiter=',', quotechar='\"');
    catObj.writerow(fieldnames);
    
    LL = [ dictionary[_k] for _k in fieldnames ];
    for _row in zip(*LL):
        catObj.writerow(_row);
    catFile.close();

# ---

def fits2array(imagename,PARAMS):
    
    # Run Sextractor for specific params..
    #
    Dsex = se.run_segobj(imagename, params=PARAMS, preset=INSTRUMENT);
    if not Dsex:
        print >> sys.stderr, "Error: SE's run raised an error. Finishing this thing...";
        sys.exit(1);
    
    # Open the output from Sextractor run, images and catalog..
    #
    header = pyfits.getheader(imagename);
    objimg = pyfits.getdata(Dsex['OBJECTS']);
    segimg = pyfits.getdata(Dsex['SEGMENTATION']);
    tbdata = fts.open_catalog(Dsex['CATALOG'])[1].data;
    
    return (header,objimg,segimg,tbdata);

# ---

def run(regionfile, imagename, outfits='.pstamp_', catrootname='objectsin_'):
    
    # Run sextractor over 'imagename' and take/open their outputs..
    #
    hdr,obj,seg,tab = fits2array(imagename,PARAMS);


    # Read DS9 regionfile..
    #
    Dds9 = ds9reg.read_cat(regionfile);

    centroids = zip(ndaa(Dds9['x'],dtype=int),ndaa(Dds9['y'],dtype=int));

    imagename_ds9 = Dds9['image'].rstrip('\n');
    if (imagename != imagename_ds9):
        print >> sys.stderr,"Image filenames do not match: %s %s" % (imagename,imagename_ds9);
        print >> sys.stderr,"Continuing with image %s " % (imagename);


    # Detect objects corresponding to given centroids..
    #
    DSE,non_detected = read_SE_data(seg,tab,centroids);


    # Concatenate objects properties from ds9 & SE..
    #
    Dout = {};
    outfits = re.sub(".fits","",imagename)+outfits;
    objIDs = DSE['NUMBER'];
    Dout['poststamp'] = [ outfits+"%s.fits" % (objIDs[i]) for i in range(len(objIDs)) ];
    Dout.update(Dds9);
    Dout.update(DSE);

    PARAMSout = ['x','y','tag','poststamp'];
    PARAMSout.extend(PARAMS)
    catname = re.sub(".fits","",imagename)+"_cat.filtered.csv";
    write_csvcat(catname,PARAMSout,Dout);
    

    # Now lets take the detected objects snapshots..
    #
    objs,hdrs = imcp.segstamp(obj, seg, header=hdr, increase=2, relative_increase=True, objIDs=objIDs);
    
    # OK. End of story. Print/write output..
    #
    outfits = re.sub(".fits","",imagename)+outfits;
    for i in range(len(objIDs)):
        if not objIDs[i]:
            print "No objects detected at position ",centroids[i];
            continue;
        outname = outfits+"%s.fits" % (objIDs[i]);
        os.system('[ -f %s ] && rm -f %s' % (outname,outname));
        pyfits.writeto( outname, objs[i], hdrs[i] );
        print "Object: %d, File: %s, Centroid: " % (objIDs[i],outname),centroids[i];



# ==========================

if __name__ == '__main__' :

    # Initializing code, command-line options..
    #
    import optparse;
    print "";
    
    usage="Usage:  %prog [options] <ds9regionfile.reg> <image.fits>"
    parser = optparse.OptionParser(usage=usage);

    parser.add_option('-o',
                        dest='outPstamp', default='.pstamp_',
                        help='Output (image) file rootname ["pstamp_"]');
    parser.add_option('-c',
                        dest='outCatalog', default='objectsin_',
                        help='Name for output identified objects (csv) catalog');

    (opts,args) = parser.parse_args();

    outfitsname = opts.outPstamp;
    catrootname = opts.outCatalog;
    
    if len(args) < 2 :
        parser.print_help();
        sys.exit(1);
        
    regionfile = args[0];
    imagename = args[1];

    run(regionfile, imagename, outfitsname, catrootname):

    print >> sys.stdout, "Done.";
    sys.exit(0);
