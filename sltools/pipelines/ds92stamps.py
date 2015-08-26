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

import logging;

#from sltools.io import check_dependencies as ckdep;
from sltools.image import sextractor as se;
from sltools.image import segobjs;
from sltools.image import imcp;
from sltools.catalog import fits_data as fts;
from sltools.catalog import ascii_data as asc;
from sltools.io import log;


# GLOBAL variables:
#
#PARAMS=['NUMBER','X_IMAGE','Y_IMAGE','X2_IMAGE','Y2_IMAGE','XY_IMAGE','FWHM_IMAGE','ELONGATION','ELLIPTICITY','ALPHA_J2000','DELTA_J2000'];

def run(D_in, objimg, segimg=None, shape=(100,100), objsfile='pstamp_', hdr=None):
    """
    Function to take snapshots from objects listed in a DS9's regions file
    
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

    from numpy import asarray;
    

    imagename_ds9 = D_in['filename'];

    # Detect objects corresponding to given centroids..
    #
    XY = zip(asarray(D_in['x']),asarray(D_in['y']));
    
    logging.debug("Input ['x'],['y'] points: %s" % (XY));
    
    if (segimg!=None):
        logging.info("Segmentation image given, looking for objects for each (x,y) on it.");
        objIDs = [ segobjs.centroid2ID(segimg,_xy) for _xy in XY ];
    else:
        logging.info("No segmentation image given.");
        objIDs = [0]*len(XY);
    
    logging.debug("ObjIDs: %s " % (objIDs));

    D_out = D_in;

    D_out['id'] = objIDs;
    D_out['segmented'] = [ _id != 0  for _id in objIDs ];
    D_out['objfile'] = [];
    del D_out['size'];

    # Now lets take the detected objects snapshots..
    #
    for i in range(len(objIDs)):

        outname = objsfile+str(i)+"_id_"+str(objIDs[i])+".fits";
        D_out['objfile'].append(outname);
        logging.info("Output objects/files (#,ID,x,y,filename): %s",(i,objIDs[i],XY[i],outname));

        _otmp,_htmp = imcp.snapshot( objimg, hdr, centroid=(XY[i][0],XY[i][1]), shape=shape );
        try:
            pyfits.writeto(outname,_otmp,_htmp);
        except IOError:
            os.remove(outname);
            pyfits.writeto(outname,_otmp,_htmp);

    
    return D_out;


# ==========================

if __name__ == '__main__' :

    import optparse;
    
    usage="\n  %prog [options] <ds9regionfile.reg> <image.fits>"
    parser = optparse.OptionParser(usage=usage);

    parser.add_option('-s', '--seg_image',
                    dest='segfile', default=None,
                    help="Segmentation image to verify the points in 'ds9regionfile'");
    parser.add_option('--shape',
                    dest='xy_sizes', default='100,100',
                    help="Output images shape [100,100]");
    parser.add_option('-o',
                    dest='objsname', default='pstamp_',
                    help="Output object filenames prefix [pstamp_]");
    parser.add_option('-c',
                    dest='catname', default='regtable_',
                    help="Name for output objects catalog [regtable_]");
    parser.add_option('-t',
                    dest='cattype', default='fits',
                    help="Type of output objects catalog (fits|csv) [fits]");

    (opts,args) = parser.parse_args();

    cattype = opts.cattype;
    catfile = opts.catname;
    objsfile = opts.objsname;
    shape = opts.xy_sizes.split(",");
    segimg = opts.segfile;
    
    if len(args) < 2 :
        parser.print_help();
        sys.exit(0);
        
    regionfile = args[0];
    imagename = args[1];

    shape = (int(shape[0]),int(shape[1]));
    

    # Logging handler:
    logfile = regionfile+'.log'
    logging = log.init(logfile,debug=True,verbose=True);


    # Read DS9 regionfile..
    #
    D_in = asc.read_ds9cat(regionfile);
    objimg = pyfits.getdata(imagename,header=False);
    if segimg:
        segimg = pyfits.getdata(segimg);
        
    D_in['filename'] = imagename;
    
    D_out = run(D_in, objimg, segimg, shape=shape, objsfile=objsfile, hdr=None);

    if cattype.upper() == 'FITS':
        outname = catfile+re.sub(".reg",".fit",regionfile);
        tbhdu = fts.dict_to_tbHDU(D_out, D_out['filename']+"_stamps", 'filename','x','y','color','segmented','id','objfile');
        try:
            tbhdu.writeto(outname);
        except IOError:
            os.remove(outname);
            tbhdu.writeto(outname);
    else:
        outname = catfile+re.sub(".reg",".csv",regionfile);
        asc.dict_to_csv(D_out,filename=outname,fieldnames=['filename','x','y','color','segmented','id','objfile']);
    
    print >> sys.stdout, "Done.";
    sys.exit(0);
