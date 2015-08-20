#!/usr/bin/env python
"""
Segment, cutout and filter a list of object positions

Function(s) defined here are meant to read informations from a DS9 region file (.reg)
and proceed with segmentation, postage stamps and catalogs creation, etc.

DS9 SAOImage imaging and data visualization : http://hea-www.harvard.edu/RD/ds9/
"""

##@package pipelines
##@file eyescan

import sys
import pyfits
from numpy import asarray
import logging
import os
import re

from sltools.io import log
from sltools.catalog import ascii_data as asc
from sltools.catalog import fits_data as fts
from sltools.image import imcp
from sltools.image import segmentation_ids as segobjs
from sltools.image import sextractor as SE
from sltools.catalog import table_matching as tbmtch

# ---

# Sextractor Params:
PARAMS=['NUMBER','ELONGATION','ELLIPTICITY','ISOAREA_IMAGE','A_IMAGE','B_IMAGE','THETA_IMAGE','X_IMAGE','Y_IMAGE','X2_IMAGE','Y2_IMAGE','XY_IMAGE']

def run(regionfile, imagefile, args={}, preset='', shape=(100,100), tablefile="objectsTb",stampsfile="pstamp_"):
    """
    Runs the pipeline for arc postage stamps creation and segmentation (in this order)
    
    For each point inside (DS9) 'regionfile' take a window of size 'shape' around it.

    Sextractor config parameters can be passed using 'args'. So far, catalog output 
    parameters are declared through the global variable 'PARAMS'.
    Use 'preset' for optimized settings. See sltools.image.sextractor for more info.
    
    Input:
     - regionfile      str : DS9 region file
     - imagefile       str : FITS filename
     - args  {'key',value} : Sextractor config params
     - preset          str : See sltools.image.sextractor
     - shape     (int,int) : Object cutouts shape (dx,dy)
     - tablefile       str : FITS table output filename
     - stampsfile      str : Output file rootnames
    
    Output:
     Write to 'tablefile' FITS table the segmented objects properties ('PARAMS')
    
    ---
    """
    
    rootname = imagefile[:-5]
    
    # Object's centroid will be placed over postamp's central point
    x_i = shape[0]/2
    y_i = shape[1]/2
    x_size,y_size = shape

    tbList = []
    flList = []
    
    # Read DS9 regionfile and input imagefile..
    #
    D_in = asc.read_ds9cat(regionfile)

    X = asarray(D_in['x'])
    Y = asarray(D_in['y'])
    centroids = XY = zip(X,Y)
    R = asarray(D_in['size'])
    logging.debug("Centroids: %s",XY)
    if len(XY) == 0:
        logging.warning("No objects in given Centroids list. Finishing.")
        return

    img,hdr = pyfits.getdata(imagefile,header=True)

    # For each "centroid":
    #
    for i in xrange(len(XY)):

        logging.info("Process %s, point: %s",i,XY[i])
        
        # Take a snapshot
        #
        obj,hdr = imcp.cutout(img,hdr,xo=X[i],yo=Y[i],x_size=x_size,y_size=y_size)
        file_i = rootname+"_pstamp_"+str(X[i])+"_"+str(Y[i])+".fits"
        pyfits.writeto(file_i,obj,hdr,clobber=True)

        logging.info("Snapshot %s created",file_i)
        flList.append(file_i)
        
        # Run Sextractor over newly created sanpshot..
        #
        Dsex = SE.run_segobj(file_i, PARAMS,preset)   ### SEXTRACTOR CALL
        segimg = pyfits.getdata(Dsex['SEGMENTATION'])
        cathdu = pyfits.open(Dsex['CATALOG'])[1]
        
        logging.info("Segmentation done.")
        logging.info("Files %s,%s,%s created",Dsex['SEGMENTATION'],Dsex['OBJECTS'],Dsex['CATALOG'])
        
        objID = segobjs.centroid2id(segimg,(x_i,y_i))

        # The following block searchs for the nearest neighbour and check if the IDs are
        # in accordance. If the given centroid is not part of an object (i.e, ID==0), 
        # the nearest object is taken in place; and if the IDs do not match and non-zero
        # valued just an warning is printed to "debug" level.
        #
        centroids_sex = zip(cathdu.data.field('X_IMAGE'),cathdu.data.field('Y_IMAGE'))
        [ logging.debug("SE Centroids: %s",_c) for _c in centroids_sex ]

        nrst_indx_dist = tbmtch.nearest_neighbour([centroids[i]],centroids_sex)
        logging.debug("Nearest point index|distance: %s (length:%s)",nrst_indx_dist,len(nrst_indx_dist))

        _indx,_dist = nrst_indx_dist[0]
        nrst_obj = cathdu.data[_indx]

        _id = nrst_obj.field('NUMBER')
        if objID == 0:
            objID = _id
            logging.warning("Object ID is 0 for Centroid:%s. Using object (ID) %s, instead.",centroids[i],_id)
        if objID != _id:
            logging.debug("Oops, two different objects were found matching point %s, objIDs %s and %s",centroids[i],objID,_id)
        
        logging.info("ObjectID: %s being readout",objID)
        tbList.append(fts.select_entries(cathdu,'NUMBER',objID))
        
        arc_i = rootname+"_nrstArc_"+str(X[i])+"_"+str(Y[i])+".fits"
        arc = imcp.select(segimg,obj,objID)
        pyfits.writeto(arc_i,arc,hdr,clobber=True)


    # Create output table
    new_tbhdu = tbList[0]
    for i in xrange(1,len(tbList)):
        new_tbhdu = fts.extend_tbHDU(new_tbhdu,tbList[i])

    tb_flList = fts.dict_to_tbHDU({'filename':flList})
    new_tbhdu = fts.merge_tbHDU(new_tbhdu,tb_flList)

    tablefile = rootname+"_arcscat.fits"
    new_tbhdu.writeto(tablefile,clobber=True)

    return

# ---
if __name__ == '__main__':
    if len(sys.argv) < 3:
        print "Usage:\n   %s  <RegionFile.reg> <Image.fits>" % (os.path.basename(sys.argv[0]));
        sys.exit(0);
    run(sys.argv[1],sys.argv[2]);
    print "Done."
    
