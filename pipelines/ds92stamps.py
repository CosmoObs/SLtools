#!/usr/bin/env python
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

#from sltools.io import check_dependencies as ckdep;
from sltools.catalog import ds9region as ds9reg;
from sltools.image import sextractor as se;
from sltools.image import segobjs;
from sltools.image import imcp;
from sltools.catalog import fits_file as fts;
from sltools.catalog import ascii_file as asc;


# GLOBAL variables:
#
PARAMS=['NUMBER','X_IMAGE','Y_IMAGE','X2_IMAGE','Y2_IMAGE','XY_IMAGE','FWHM_IMAGE','ELONGATION','ELLIPTICITY','ALPHA_J2000','DELTA_J2000'];


# ---

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

# ---

# ---

def run(regionfile, imagename, outname_stamps='pstamp_', outname_cat='_cat.filtered.csv', sex_preset=''):
    """ Function to take snapshots from objects listed in a DS9's regions file
    
    Sextractor is run over given 'imagename' and take the snapshots from objects listed
    in DS9's 'regionfile'.
    
    Input:
     - regionfile : str
        DS9' regions file
     - imagename : str
        FITS image filename
     - outname_stamps : str
        String to add to output name of the poststamp images.
     - outname_cat : str
     
    Output:
     -> <bool>
    
    """

    # Read DS9 regionfile..
    #
    D_in = ds9reg.read_cat(regionfile);

#    centroids = zip(ndaa(Dds9_in['x'],dtype=float),ndaa(Dds9_in['y'],dtype=float));

    imagename_ds9 = D_in['image'].rstrip('\n');
    if (imagename != imagename_ds9):
        print >> sys.stderr,"Warning: Image filenames do not match: %s %s" % (imagename,imagename_ds9);
        print >> sys.stderr,"Warning: Continuing with image %s " % (imagename);


    # Run sextractor over 'imagename' and take/open their outputs..
    #
    Dsex = se.run_segobj(imagename, params=PARAMS, preset=preset);
    if not Dsex:
        print >> sys.stderr, "Error: SE's run raised an error. Finishing this thing...";
        sys.exit(1);
    
    # Open the output from Sextractor run, images and catalog..
    #
    header = pyfits.getheader(imagename);
    objimg = pyfits.getdata(Dsex['OBJECTS']);
    segimg = pyfits.getdata(Dsex['SEGMENTATION']);
    tbdata = fts.open_catalog(Dsex['CATALOG'])[1].data;

    # Detect objects corresponding to given centroids..
    #
    x = D_in['x'];
    y = D_in['y'];
    objIDs = segobjs.centroid2ID(segimg,zip(x,y));

    # Now lets take the detected objects snapshots..
    #
    outname_1 = string.replace(".fits","",imagename)+"_sextracted_";
    outname_0 = string.replace(".fits","",imagename)+"_noxtracted_";
    idXY = zip(objIDs,x,y);
    for i in range(len(idXY)):
        if idXY[i][0]:
            _otmp,_htmp = imcp.segstamp(objimg,segimg,hdr,objIDs=[idXY[i][0]],increase=2);
            pyfits.writeto(outname_1+str(i)+str(idXY[i][0])+".fits",_otmp[0],_htmp[0]);
        else:
            _otmp,_htmp = imcp.snapshot(objimg,hdr,(idXY[i][1],idXY[i][2]),shape=(100,100))
            pyfits.writeto(outname_0+str(i)+"0"+".fits",_otmp,_htmp);


    # Detect objects corresponding to given centroids..
    #
#    DSE,non_detected = read_SE_data(seg,tab,centroids);
    # ..update ds9 dictionary/cat:
#    DS9 = {};
#    valid_inds = list(set(range(len(centroids)))-set(non_detected['ind']));
#    for _key in Dds9_in.keys():
#        DS9[_key] = [ Dds9_in[_key].pop(i) for i in indices ]

#    del Dds9_in;
    DS9_detected = {}
    for i in range(len(idXY)):
        if idXY[i][0]:
            DS9_detected['x']

<<<<<<< HEAD
    imagename_ds9 = Dds9['image'];
    if (imagename != imagename_ds9):
        print >> sys.stderr,"Image filenames do not match: %s %s" % (imagename,imagename_ds9);
        print >> sys.stderr,"Continuing with image %s " % (imagename);
=======
    
    
    # Now lets take the detected objects snapshots..
    #
    objIDs = DSE['NUMBER'];
    objs,hdrs = imcp.segstamp(obj, seg, header=hdr, increase=2, relative_increase=True, objIDs=objIDs);
>>>>>>> 9dae1419c5627c775db5da7455e87343c75a0b67


    # OK. End of story. Lets write down the outputs..

    Dout = {};

    # Now, concatenate objects properties from ds9 & SE...
    #
    Dout.update(DS9);
    Dout.update(DSE);

    PARAMSout = ['x','y','color','poststamp'];
    PARAMSout.extend(PARAMS)

    outname_cat = re.sub(".fits","",imagename)+outname_cat;

    # ...and write the csv catalog:
    asc.dict_to_csv(Dout,PARAMSout,outname_cat);
    
    return True;


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
<<<<<<< HEAD
                        help='Output (image) file name appended to image file rootname [".pstamp_"]');
    parser.add_option('-c',
                        dest='outCatalog', default='objectsin_',
                        help='Name for output identified objects (csv) catalog ["objectsin_"]');
=======
                        help='Output (image) file name to append [".pstamp_"]');
    parser.add_option('-c',
                        dest='outCatalog', default='_cat.filtered.csv',
                        help='Name for output identified objects (csv) catalog');
    parser.add_option('--sex_preset',
                    dest='preset', default='',
                    help='Preset Sextractor configuration. See sltools.image.sextractor for more');
>>>>>>> 9dae1419c5627c775db5da7455e87343c75a0b67

    (opts,args) = parser.parse_args();

    outfitsname = opts.outPstamp;
    catrootname = opts.outCatalog;
    sex_preset = opts.preset;
    
    if len(args) < 2 :
        parser.print_help();
        sys.exit(1);
        
    regionfile = args[0];
    imagename = args[1];

<<<<<<< HEAD
    run(regionfile, imagename, outfitsname, catrootname)
=======
    run(regionfile, imagename, outfitsname, catrootname,sex_preset):
>>>>>>> 9dae1419c5627c775db5da7455e87343c75a0b67

    print >> sys.stdout, "Done.";
    sys.exit(0);
