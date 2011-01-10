#!/usr/bin/env python
import sys;

"""Module to segment and filter a list of positions on given image file"""

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
from sltools.pipelines import table_matching as tm;

# ---

def run(regionfile, imagefile, shape=(100,100), stampsfile="pstamp_"):
    """
    Input:
     - regionfile : str
        DS9 region file
     - imagefile : str
        FITS filename
     
    Output:
     - tbhdu
    
    ---
    """
    
    # Sextractor Params:
    PARAMS=['NUMBER','FWHM_IMAGE','ELONGATION','ELLIPTICITY','X_IMAGE','Y_IMAGE','ALPHA_SKY','DELTA_SKY','X2_IMAGE','Y2_IMAGE','XY_IMAGE'];

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

    for i in xrange(len(centroids)):

        logging.info("Process %s, point: %s",i,centroids[i]);
        
        # Take a snapshot for each centroid..
        #
        _obj,_hdr = imcp.snapshot( img, hdr, centroid=centroids[i], shape=shape );
        file_i = stampsfile+str(i)+".fits";
        pyfits.writeto(file_i,_obj,_hdr)

        logging.info("Snapshot %s created",file_i);
        _flList.append(file_i);
        
        # Run Sextractor over newly created sanpshot..
        #
        _Dsex = SE.run_segobj(file_i, PARAMS,preset='HST');   ### SEXTRACTOR IS CALLED HERE!
        segimg = pyfits.getdata(_Dsex['SEGMENTATION']);
        cathdu = pyfits.open(_Dsex['CATALOG'])[1];
        logging.info("Segmentation done.");
        logging.debug("Files %s,%s,%s",_Dsex['SEGMENTATION'],_Dsex['OBJECTS'],_Dsex['CATALOG']);
        
        # The following block, search for the nearest-neighbour, could be done only
        # under the exception "objID(=segimg[y_i,x_i])"
        #
        centroids_sex = zip(cathdu.data.field('X_IMAGE'),cathdu.data.field('Y_IMAGE'));
        nrst_indx_dist = tm.nearest_neighbour([centroids[i]],centroids_sex);
        logging.debug("Nearest-point index/distance: %s (length:%s)",nrst_indx_dist,len(nrst_indx_dist));
        _indx,_dist = nrst_indx_dist[0];

        nrst_obj = fts.select_rows(cathdu,_indx);
        objID = segobjs.centroid2ID(segimg,(x_i,y_i));

        _id = nrst_obj.data.field('NUMBER');
        _xy = zip([float(nrst_obj.data.field('X_IMAGE'))],[float(nrst_obj.data.field('Y_IMAGE'))]);
        if objID == 0:
            objID = _id;
            logging.warning("Object ID is 0 for Centroid:%s. Using objID %s, centered:%s, instead.",centroids[i],_id,_xy);
        if objID != _id:
            logging.debug("Oops, two different objects were found matching point %s, objIDs %s and %s",centroids[i],objID,_id);
            
        logging.info("ObjectID: %s, Position: (%s,%s), Value(id|seg): %s",objID,x_i,y_i,segimg[y_i,x_i]);
        _tbList.append(fts.select_entries(cathdu,'NUMBER',objID));

    _dtmp = {'filename':_flList};
    tb_flList = fts.dict_to_tbHDU(_dtmp);
    
    new_tbhdu = _tbList[0];
    for i in xrange(1,len(_tbList)):
        new_tbhdu = fts.extend_tbHDU(new_tbhdu,_tbList[i]);

    new_tbhdu = fts.merge_tbHDU(new_tbhdu,tb_flList);

    rootname = re.sub(".fits","",imagefile);
    try:
        os.mkdir('stamps_'+rootname);
    except:
        os.system('rm -f stamps_%s/*.*'%(rootname));
    os.system('mv %s* stamps_%s/.'%(stampsfile,rootname));
    
    new_tbhdu.writeto(rootname+'_tbarcs.fit');

    return;

# ---

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print "Usage:\n   %s  <RegionFile.reg> <Image.fits>" % (os.path.basename(sys.argv[0]));
        sys.exit(0);
    run(sys.argv[1],sys.argv[2]);
    print "Done."
