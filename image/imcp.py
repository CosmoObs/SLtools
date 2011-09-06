#!/usr/bin/env python
# ==================================
# Authors:
# Carlos Brandt - chbrandt@lncc.br
# ==================================

"""Package to make FITS image copy/cuts and update their header"""

##@package image
##@file imcp
#
#
# This package contains functions to cutout smaller regions
# from a original image, paste the image into a bigger array,
# as well as any other function for image copy manipulation.
#
# Executable package: YES
#
# To see the package help message, type:
#
# > python imcp.py --help
#
# 'imcp' has one mandatory argument: the (input) image FITS file.
# If none of the options are given, imcp outputs a file named
# 'output.fits' with 1/4th of the original image area, with the 
# same central point/pixel.
#
# > python imcp.py input_file.fits
#
# The options available make possible to define the central point
# where the cut will be centered, as well as the side lengths for
# output image. These values can be passed in 'pixel' or 'degrees'
# units.

import sys;
import logging;

import os
import pyfits
import string
import commands
import numpy as np;

import sltools;
from sltools.string import regexp;
from sltools.image import segmentation_ids;
from sltools.image import header_funcs as hf;


# =================================================================================================================
def cutout( img, hdr=None, coord_unit='pixel', xo=0, yo=0, size_unit='pixel', x_size=0, y_size=0, mask=None ):
    """
    Do a snapshot from given fits image.

    cutout( ndarray, ...) -> (ndarray, header)

    The function can deal with 'pixel' and 'degrees' cordinate units. As input, besides the image FITS 
    filename (input_file), it receives a central position (xo,yo) and side lengths (x_size,y_size) for
    output image position and dimensioning. Function returns two values: a numpy.array with resultant 
    image pixel intensities and respective header object.

    In addition, an image mask information can be given to use as interest region selection. 'mask' is
    expected to be a numpy.where like structure (i.e, mask=numpy.where()). If given, returned image will
    have all pixels null except the ones listed in 'mask' parameter.

    If 'xo=0' and 'yo=0', given image (input_file) central pixel will be chosen as (xo,yo).
    If 'x_size=0' and 'y_size=0', half length of each side will be used for output dimensions.


    Input:
     - img        numpy.ndarray : Image array (ndim=2,dtype=float)
     - hdr        pyfits.header : Image (FITS) header
     - coord_unit           str : Position (xo,yo) units ('pixel','degrees')
     - xo                   int : Centroid (horizontal) position for cutting window
     - yo                   int : Centroid (vertical) position for cutting window
     - size_unit            str : Window sizes (x_size,y_size) units ('pixel','degrees')
     - x_size               int : Horizontal cutting window size
     - y_size               int : Vertical cutting window size
     - mask   (ndarray,ndarray) : Tuple with index arrays (Output like numpy.where())

    Output:
     - (ndarray, header) : Resultant image array and (updated) header instance
    
    ---
    """

    image = img;
    
    logging.debug("Input parameters (type(image),header,coord_unit,xo,yo,size_unit,x_size,y_size,mask): %s",(type(image),hdr,coord_unit,xo,yo,size_unit,x_size,y_size,mask));

    xo=float(xo);
    yo=float(yo);
    x_size=float(x_size);
    y_size=float(y_size);

    imagem = image;
#    if hdr:
#        hdr = hdr.copy();
    
    # Initialize some variables..
    #
    x_diff = 0;   x_edge = 0;
    y_diff = 0;   y_edge = 0;
    x_fin = 0;   x_ini = 0;
    y_fin = 0;   y_ini = 0;


    if ( hdr ):
        dimpix = hf.get_pixelscale( hdr );
        logging.info("Pixel_scale: %s",dimpix);

    y_img_size, x_img_size = imagem.shape;
    logging.debug("Input image shape: %s",(x_img_size,y_img_size));

    # Get the side sides (at least, 1!) and transform for the size_unit if necessary..
    #
    x_cut_size = max( 1, int(float(x_size)) );
    y_cut_size = max( 1, int(float(y_size)) );

    if ( size_unit == 'degrees' ):
        x_cut_size = int(float(x_cut_size)/dimpix);
        y_cut_size = int(float(y_cut_size)/dimpix);


    # And if no side size was given, define a default value correspondig to half of original image..
    #
    if not ( x_size ):
        x_cut_size = int(x_img_size/2);
        logging.warning("'x_size' not given. Using half of image x side(%d)", x_cut_size);

    if not ( y_size ):
        y_cut_size = int(y_img_size/2);
        logging.warning("'y_size' not given. Using half of image y side(%d)", y_cut_size);

    logging.debug("Output image shape: %s",(x_cut_size,y_cut_size));

    # Check requested output.vs.input image sizes..
    #
    if ( x_cut_size == x_img_size and y_cut_size == y_img_size ):
        logging.warning("Requested output sizes are the same as input image. Returning image and header unchanged.");
        return (imagem,hdr);

    # Verify central coordinates values..
    #
    if ( coord_unit == 'pixel' and (xo != 0 and yo != 0) ):
        x_halo = int(float(xo));
        y_halo = int(float(yo));

    elif ( coord_unit == 'degrees' and (xo != 0 and yo != 0) ):
        _out = commands.getoutput( 'sky2xy %s %s %s' % (input_file, xo, yo) );
        _out = string.split(_out);
        x_halo = int(float(_out[-2]))-1;
        y_halo = int(float(_out[-1]))-1;

    elif ( xo == 0 and yo == 0 ):
        x_halo = int(x_img_size/2);
        y_halo = int(y_img_size/2);
        logging.warning("No central coordinates were given for snapshot. Using image central pixel as the cut center.");
    else:
        logging.error("Central positioning is out of valid values.");
        return (False,False);

    logging.debug("Central point for output image: %s",(x_halo,y_halo));

    # Define the images (in/out) slices to be copied..
    #
    x_ini = x_halo - int(x_cut_size/2) #-1;
    x_fin = x_ini + x_cut_size;
    y_ini = y_halo - int(y_cut_size/2) #-1;
    y_fin = y_ini + y_cut_size;

    x_ini_old = max( 0, x_ini );   x_fin_old = min( x_img_size, x_fin );
    y_ini_old = max( 0, y_ini );   y_fin_old = min( y_img_size, y_fin );

    x_ini_new = abs( min( 0, x_ini ));   x_fin_new = x_cut_size - (x_fin - x_fin_old);
    y_ini_new = abs( min( 0, y_ini ));   y_fin_new = y_cut_size - (y_fin - y_fin_old);

    # If header, update pixel<->sky information..
    #
    if ( hdr ):
        hdr = hf.update_coordinates(hdr.copy(), x_ini, y_ini);
        hdr.update('NAXIS1',x_cut_size);
        hdr.update('NAXIS2',y_cut_size);

    # Initialize new image, and take all index list..
    #
    imagemnova = np.zeros( (y_cut_size,x_cut_size), dtype=imagem.dtype );
    ind_z = np.where(imagemnova == 0);

    # Copy requested image slice..
    #
    imagemnova[ y_ini_new:y_fin_new, x_ini_new:x_fin_new ] = imagem[ y_ini_old:y_fin_old, x_ini_old:x_fin_old ];

    # If 'mask', maintain just "central" object on it..
    #
    if ( mask ):
        msk = ( mask[0]-y_ini, mask[1]-x_ini )

        zip_m = zip( msk[0], msk[1] );
        zip_z = zip( ind_z[0], ind_z[1] );

        L = list(set(zip_z) - set(zip_m));

        try:
            ind_0, ind_1 = zip(*L);
            indx = ( np.array(ind_0), np.array(ind_1) );
            imagemnova[ indx ] = 0;
        except:
            pass;

    return (imagemnova, hdr);
    
# ---

# ==========================================================================================
def segstamp(segimg, objID, objimg=None, hdr=None, increase=0, relative_increase=False, connected=False):
    """
    Identify objects on given images by their IDs and return object images

    segstamp( segimg, objID, ...) -> (ndarray,header)

    By default, if 'objIDs' is not given, postamp will scan segmentation image 
    'seg_img' for the list of object ID numbers. If 'objIDs' is given, those IDs 
    will be used for object poststamps creation.

    'increase' and 'relative_increase' define whether the poststamps will have a 
    size different from object-only dimensions or just the object (pixels) itself, 
    and, also, if this increasement value is a multiplicative factor (True) or an 
    additive one (False).

    Since a list with object IDs is given, a list with arrays, with each IDentified
    object, is returned.

    Input:
     - segimg          ndarray : Segmented image (ndim=2,dtype=int)
     - objIDs            [int] : List with object IDs of interest in 'segimg'.
     - objimg          ndarray : Objects image (ndim=2,dtype=float)
     - hdr       pyfits.header : FITS header to be updated for each poststamp
     - increase          float : Value for poststamp resizing (> 0)
     - relative_increase  bool : Is 'increase' a additive value (default,False) or multiplicative one(True)?
        

    Output:
     - (ndarray,header)  : Image (poststamp) array and corresponding header

    ---
    """


    _id = objID;
    ind = segobjs.create_IDmask(segimg, objID);

    y_min = min( ind[0] );
    x_min = min( ind[1] );

    y_idx = ind[0] - y_min;
    x_idx = ind[1] - x_min;

    y_size = max( y_idx ) + 1;
    x_size = max( x_idx ) + 1;
    
   
    if (connected == True ):		
      ind=elliminate_disconected(ind)
 
    # Central pixel on original image:
    yo = y_size/2 + y_min;
    xo = x_size/2 + x_min;

    if ( increase != 0 ):		
        if (relative_increase == True):
            x_size = x_size*increase;
            y_size = y_size*increase;
        else:
            x_size = x_size + 2*increase;
            y_size = y_size + 2*increase;

    if objimg!=None:
        image_out, hdr = cutout( objimg, hdr, xo=int(xo), yo=int(yo), x_size=int(x_size), y_size=int(y_size), mask=ind );
    else:
        image_out, hdr = cutout( segimg, hdr, xo=int(xo), yo=int(yo), x_size=int(x_size), y_size=int(y_size), mask=ind );
    
    return ( image_out, hdr );

# ---
def select(segimg,objimg,objID):
    """
    Select pixels from a specific ID and outputs an array with the same input size
    
    Input:
     - segimg  ndarray : Segmented image array (ndim=2,dtype=int)
     - objimg  ndarray : Objects (observed) image array (ndim=2,dtype=float)
     - objID       int : Object ID to request output for
    
    Output:
     - outimg  ndarray : Same size (segimg) array with non-null pixels only for objID
    
    ---
    """

    outimg = np.zeros(segimg.shape)
    
    indxs = np.where(segimg==objID)
    if not len(indxs):
        return None
    
    outimg[indxs] = objimg[indxs]
    
    return outimg

# ---

def separate_disconected(ind, high_area=False):

    """
    From a list of points that belongs to objects it separates in groups of connected objects. If high_area=False return a list which first index represents each object , the second and third its coordinates. if high_area=True  returns the only a list with the coordinates of the object that has greater area.

   
    Input:
     - ind : array like
        the list of points where ind[0] is a list with the first coordinates (e.g. x) and ind[1] the second coordinates (e.g. y)
     - high_area : bool
        enables a criteria to return only the coordinates of the object that has greater area
   
    Output:
     - (nparray)  : if high_area=False  is list which first index represent the  object.  The second  and third index represents  the first and second coordinates. if high_area=True  is list which first and second index represents lists the first and second coordinates

    ---
    """

    p=[]
    for i in range(0,len(ind[0])):
        p_aux=[ind[0][i],ind[1][i]]
        p.append(p_aux)
   
  
    Objects=[]
    objects=[]
    while len(p)>0:
        p_test=[p[0]]
        while len(p_test)>0:
            p_center=p_test[0]
            p_neighbors=[[p_center[0],p_center[1]],[p_center[0]+1,p_center[1]],[p_center[0]+1,p_center[1]+1],[p_center[0]+1,p_center[1]-1],[p_center[0]-1,p_center[1]],[p_center[0]-1,p_center[1]+1],[p_center[0]-1,p_center[1]-1],[p_center[0],p_center[1]-1],[p_center[0],p_center[1]+1]]
            for i in range(0,len(p_neighbors)):
                if (p_neighbors[i] in p) and not(p_neighbors[i] in objects):
                    objects.append(p_neighbors[i])
                    p_test.append(p_neighbors[i])
                    p.remove(p_neighbors[i])  
            p_test.remove(p_test[0])
        Objects.append(objects)
        objects=[]

    if high_area==False:
        return Objects
    else:
        # criteria to select an interest object
   
        Area_max=len(Objects[0])
        id_max=0
        for i in range(0,len(Objects)):
            if len(Objects[i])>Area_max:
                Area_max=len(Objects[i])
                id_max=i
        ind_new=[[-100],[-1000]]
        for i in range(0,len(Objects[id_max])):
            ind_new[0].append(Objects[id_max][i][0])
            ind_new[1].append(Objects[id_max][i][1])
        ind_new[0].remove(-100)
        ind_new[1].remove(-1000)
        ind_n=[np.array(ind_new[0]),np.array(ind_new[1])]
        return ind_n          

# ---

# \cond
# ============================================================================
if __name__ == "__main__" :

    from optparse import OptionParser;
    from sltools.io import log;
    
    logging = log.init(debug=True,verbose=False);
    
    usage="\n  %prog image.fits x y  [options]"
    parser = OptionParser(usage=usage);

    parser.add_option('-j', action='store_true',
                dest='is_coord_deg', default=False,
                help="If given (x,y) coordinates are in degrees, default is False");
    parser.add_option('-J', action='store_true',
                dest='is_shape_deg', default=False,
                help="If given 'shape' (x,y sizes) are in degrees, default is False");
    parser.add_option('-s',
                dest='dim', default='250,250',
                help="Output stamps shape default is 250,250");
    parser.add_option('--increase',
                dest='incr', default=0,
                help="Amount of border/enlarge of object poststamp. It is an additive or multiplicative factor. See 'relative_increase' flag");
    parser.add_option('--relative_increase', action='store_true',
                dest='rel_incr', default=False,
                help="Turn on the multiplicative factor for increase image, otherwise, 'increase' will be treated as an absolute value.");
    parser.add_option('-o',
                dest='sufx', default='pstamp',
                help="Output images (stamps) filename suffix [pstamp_]");
    parser.add_option('--no-header', action='store_false',
                dest='use_header', default=True,
                help="Do not use input image header.");
    parser.add_option('--outdir',
                dest='dir_name', default='out_imcp',
                help="Directory name to store output images(poststamps)");

    parser.add_option('--segimg',
                dest='segimg', default=None,
                help="Segmented image (e.g, SE's SEGMENTATION output)");
    parser.add_option('--objIDs',
                dest='IDs',default='',
                help="Comma-separated list of object IDs found on 'segimg' to extract");

    
    (opts,args) = parser.parse_args();

    segimg = opts.segimg;
    IDs = opts.IDs;
    incr = opts.incr;
    rel_incr = opts.rel_incr;
    use_header = opts.use_header;
    dir_name = opts.dir_name;

    if ( len(args) < 3 ):
        parser.print_help();
        sys.exit(0);

    infits = args[0];
    x = args[1];
    y = args[2];
    coord_deg = opts.is_coord_deg;
    shape_deg = opts.is_shape_deg;
    dx,dy = string.split(opts.dim,sep=',');
    outfits = opts.sufx;

    coord = 'degrees' if coord_deg else 'pixel';
    size = 'degrees' if shape_deg else 'pixel';        

    hdulist = pyfits.open(infits,memmap=True);
    if (use_header):
#        image, hdr = pyfits.getdata(infits,header=True);
        image = hdulist[0].data;
        hdr = hdulist[0].header;
    else:
#        image = pyfits.getdata(infits, header=False);
        image = hdulist[0].data;
        hdr = None;

    try:
        os.mkdir( dir_name );
    except OSError:
        pass;

    owd = os.getcwd()+"/";
    os.chdir( dir_name );

    if (segimg):
        IDs = regexp.str2lst(IDs);
        imglst, hdrlst = segstamp(image, segimg, hdr=hdr, increase=incr, relative_increase=rel_incr, objIDs=objIDs)
    else:
        imgcut, hdr = cutout( image, hdr=hdr, coord_unit=coord, xo=x, yo=y, size_unit=size, x_size=dx, y_size=dy );
        IDs = ['0'];
        imglst = [imgcut];
        hdrlst = [hdr];

    for i in range(len(IDs)):
        outfits = outfits+"%s.fits" % (IDs[i]);
        os.system('[ -f %s ] && rm -f %s' % (outfits,outfits));
        pyfits.writeto( outfits, imglst[i], hdrlst[i] );

    print >> sys.stdout, "Done.";
    os.chdir( owd );
    sys.exit(0);





# ------------------
# \endcond
