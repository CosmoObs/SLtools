#!/usr/bin/env python
import sys;

"""Module to segment and filter a list of positions on given image file"""

##@ eyescan

import pyfits;
import numpy as np;
from numpy import asarray;
import logging;
import os;
import re;

import sltools;
from sltools.catalog import ascii_data as asc;
from sltools.catalog import fits_data as fts;
from sltools.image import imcp,segobjs,combine;
from sltools.image import sextractor;
from sltools.pipelines import table_matching as tbmtch;

# ---

def whole_image(X,Y,imagefile,PARAMS,preset):

    centroids = zip(X,Y);
    if len(centroids) == 0:
        logging.warning("No objects in given (X,Y lists). Finishing.");
        return False;

    # Run Sextractor over the whole image:
    _Dsex = sextractor.run_segobj(imagefile,PARAMS,preset=preset);
    segimg = pyfits.getdata(_Dsex['SEGMENTATION']);
    objimg = pyfits.getdata(_Dsex['OBJECTS']);

    objIDs = [ segimg[o_o[1],o_o[0]]  for o_o in centroids ];
    objIDs = list(set(objIDs)-set([0]));

    # Take only the objects of interest...
    # from image:
    selobjimg = imcp.copy_objects(objimg,segimg,objIDs);
    del segimg,objimg;
    # and catalog:
    cathdu = pyfits.open(_Dsex['CATALOG'])[1];
    selobjhdu = fts.select_entries(cathdu,'NUMBER',*objIDs);
    del cathdu;

    # Re-identify each object:
    for i in range(len(objIDs)):
        selobjhdu.data.field('NUMBER')[i] = i;
    
    return (selobjimg,selobjhdu);


# ---

def run(regionfile, imagefile, args={}, preset='HST_Arcs', shape=(100,100)):# , tablefile="objectsTb"):# ,stampsfile="pstamp_"):
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
    _imList = [];
    
    # Read DS9 regionfile and input imagefile..
    #
    D_in = asc.read_ds9cat(regionfile);
    X = asarray(D_in['x']);
    Y = asarray(D_in['y']);

    centroids = zip(X,Y);

    ####################################################################
    # Segment the whole image at once
    # Needs:
    # X , Y , image(FITS) filename and Sextractor arguments/inputs
    selobjimg_W,selobjhdu_W = whole_image(X,Y,imagefile,PARAMS,preset);
    selobjhdu_W.writeto('Whole_image.fit',clobber=True);
    pyfits.writeto('Whole_image.fits',selobjimg_W,clobber=True);
    ####################################################################
        
    img,hdr = pyfits.getdata(imagefile,header=True);
    logging.debug("Centroids: %s",centroids);

    rootname = re.sub(".fits","",imagefile);

    # For each centroid, do:
    #
    i=-1;
    xy = centroids[:];
    for o_o in centroids:
        i+=1;
        logging.info("Process %s, point: %s",i,o_o);
        
        # Take a snapshot
        #
        _obj,_hdr = imcp.snapshot( img, hdr, centroid=o_o, shape=shape );
        file_i = rootname+"_ps"+str(i)+".fits";
        pyfits.writeto(file_i,_obj,_hdr,clobber=True);
        del _obj,_hdr;
        
        logging.info("Poststamp %s created",file_i);
        
        # Run Sextractor over newly created sanpshot..
        #
        _Dsex = sextractor.run_segobj(file_i, PARAMS,preset=preset);   ### SEXTRACTOR CALL
        segimg = pyfits.getdata(_Dsex['SEGMENTATION']);
        objimg = pyfits.getdata(_Dsex['OBJECTS']);
        cathdu = pyfits.open(_Dsex['CATALOG'])[1];
        
        objID = segimg[y_i,x_i];
        
        if not objID:
            lixo = xy.remove(o_o);
            continue;
            
        logging.info("ObjectID: %s being readout",objID);

        _tbList.append(fts.select_entries(cathdu,'NUMBER',objID));
        _flList.append(file_i);
        _imList.append(imcp.copy_objects(objimg,segimg,[objID]));

    # Initialize output table and image..
    #
    selobjhdu_S = _tbList[0];
    selobjimg_S = np.zeros(img.shape,dtype=img.dtype);    
    selobjimg_S = combine.add_images(selobjimg_S,_imList[0],x=xy[0][0],y=xy[0][1]);
    #
    # and do the same for each object
    #
    for i in xrange(1,len(_tbList)):
        selobjhdu_S = fts.extend_tbHDU(selobjhdu_S,_tbList[i]);
        selobjimg_S = combine.add_images(selobjimg_S,_imList[i],x=xy[i][0],y=xy[i][1]);

    # Write down the FITS catalog..
    #
    tb_flList = fts.dict_to_tbHDU({'filename':_flList});
    selobjhdu_S = fts.merge_tbHDU(selobjhdu_S,tb_flList);

    # And the FITS image, composed by the well segmented objects..
    #
    selobjhdu_S.writeto('Stamps_compose.fit',clobber=True);
    pyfits.writeto('Stamps_compose.fits',selobjimg_S,clobber=True)

    return;

# ---

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print "Usage:\n   %s  <RegionFile.reg> <Image.fits>" % (os.path.basename(sys.argv[0]));
        sys.exit(0);
    run(sys.argv[1],sys.argv[2]);
    print "Done."
    
