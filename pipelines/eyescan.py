#!/usr/bin/env python
import sys;

"""Module to segment and filter a list of positions on given image file."""

##@ eyescan

import pyfits;
from numpy import asarray;
import logging;
import os;
import re;

import sltools;
from sltools.catalog import ascii_data as asc;
from sltools.catalog import fits_data as fts;
from sltools.image import imcp,segobjs;
from sltools.image import sextractor as SE;
from sltools.pipelines import table_matching as tbmtch;

# ---

def run(regionfile, imagefile, args={}, preset='', shape=(100,100), tablefile="objectsTb",stampsfile="pstamp_"):
    """
    Runs the pipeline for postamps creation and segmented
    
    run( regionfile, imagefile [,...] )
    
    For each point listed inside DS9 'regionfile', considered object 
    centroids in 'imagefile', a snapshot and object segmentation are 
    done; shape and rootname (suffixed by "objID", fits extension) 
    are passed through 'shape' and 'stampsfile', resp.
    
    Sextractor config parameters can be passed using 'args'. So far, 
    catalog output parameters are hardcoded.
    Use 'preset' for optimized settings. See sltools.image.sextractor
    for more info.
    
    Input:
     - regionfile : str
        DS9 region file
     - imagefile : str
        FITS filename
     - args : {'option':'value',}
        Sextractor config params
     - preset : str
        See sltools.image.sextractor
     - shape : (int,int)
        Object cutouts shape (pixels,pixels)
     - tablefile : str
        FITS table output filename
     - stampsfile : str
        Output file rootnames
    
    Output:
     - tbhdu : pyfits.BinTableHDU
        Output (fits) table with objects' computed parameters
    
    ---
    """
    
    # Sextractor Params:
    PARAMS=['NUMBER','ELONGATION','ELLIPTICITY','ISOAREA_IMAGE','A_IMAGE','B_IMAGE','THETA_IMAGE','X_IMAGE','Y_IMAGE','X2_IMAGE','Y2_IMAGE','XY_IMAGE']

    # Object's centroid will be placed over postamp's central point
    x_i = shape[0]/2;
    y_i = shape[1]/2;

    _tbList = [];
    _flList = [];
    
    # Read DS9 regionfile and input imagefile..
    #
    D_in = asc.read_ds9cat(regionfile);

    centroids = zip(asarray(D_in['x']),asarray(D_in['y']));
    if len(centroids) == 0:
        logging.warning("No objects in given Centroids list. Finishing.");
        return;

    img,hdr = pyfits.getdata(imagefile,header=True);
    logging.debug("Centroids: %s",centroids);

    # For each centroid...
    #
    for i in xrange(len(centroids)):

        logging.info("Process %s, point: %s",i,centroids[i]);
        
        # Take a snapshot
        #
        _obj,_hdr = imcp.snapshot( img, hdr, centroid=centroids[i], shape=shape );
        file_i = stampsfile+str(i)+".fits";
        try:
            pyfits.writeto(file_i,_obj,_hdr)
        except:
            os.system('rm -f %s' % (file_i));
            pyfits.writeto(file_i,_obj,_hdr)

        logging.info("Snapshot %s created",file_i);
        _flList.append(file_i);
        
        # Run Sextractor over newly created sanpshot..
        #
        _Dsex = SE.run_segobj(file_i, PARAMS,preset);   ### SEXTRACTOR CALL
        segimg = pyfits.getdata(_Dsex['SEGMENTATION']);
        cathdu = pyfits.open(_Dsex['CATALOG'])[1];
        
        logging.info("Segmentation done.");
        logging.info("Files %s,%s,%s created",_Dsex['SEGMENTATION'],_Dsex['OBJECTS'],_Dsex['CATALOG']);
        
        objID = segobjs.centroid2ID(segimg,(x_i,y_i));

        # The following block searchs for the nearest neighbour and check if the IDs are
        # in accordance. If the given centroid is not part of an object (i.e, ID==0), 
        # the nearest object is taken in place; and if the IDs do not match and non-zero
        # valued just an warning is printed to "debug" level.
        #
        centroids_sex = zip(cathdu.data.field('X_IMAGE'),cathdu.data.field('Y_IMAGE'));
        [ logging.debug("SE Centroids: %s",_c) for _c in centroids_sex ];

        nrst_indx_dist = tbmtch.nearest_neighbour([centroids[i]],centroids_sex);
        logging.debug("Nearest point index|distance: %s (length:%s)",nrst_indx_dist,len(nrst_indx_dist));

        _indx,_dist = nrst_indx_dist[0];

        nrst_obj = fts.select_rows(cathdu,_indx);
        _id = nrst_obj.data.field('NUMBER');

        if objID == 0:
            objID = _id;
            logging.warning("Object ID is 0 for Centroid:%s. Using object (ID) %s, instead.",centroids[i],_id);
        if objID != _id:
            logging.debug("Oops, two different objects were found matching point %s, objIDs %s and %s",centroids[i],objID,_id);
            
        logging.info("ObjectID: %s being readout",objID);
        _tbList.append(fts.select_entries(cathdu,'NUMBER',objID));

    
    new_tbhdu = _tbList[0];
    for i in xrange(1,len(_tbList)):
        new_tbhdu = fts.extend_tbHDU(new_tbhdu,_tbList[i]);

    tb_flList = fts.dict_to_tbHDU({'filename':_flList});
    new_tbhdu = fts.merge_tbHDU(new_tbhdu,tb_flList);

    new_tbhdu.writeto(tablefile+'.fit',clobber=True);

    return;

# ---

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print "Usage:\n   %s  <RegionFile.reg> <Image.fits>" % (os.path.basename(sys.argv[0]));
        sys.exit(0);
    run(sys.argv[1],sys.argv[2]);
    print "Done."
    
