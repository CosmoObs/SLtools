#!/usr/bin/env python

import sys;
import os;
import pyfits;
from pylab import *;
from pywcs import WCS;

# SLTOOLS:
from sltools.image import imcp;
from sltools.image import sextractor as SE;


# -----

class CoordPoints:
    """Data structure to store and transform px/sky coordinates"""
    def __init__(self,x=[],y=[]):
        self.xlist = x[:]
        self.ylist = y[:]
    def __call__(self,event):
        self.xlist.append(event.xdata)
        self.ylist.append(event.ydata)
    def xy_coords(self):
        return zip(self.xlist, self.ylist)
    def dg_coords(self,hdr):
        pw = WCS(header=hdr)
        _x, _y = pw.wcs_pix2sky(self.xlist, self.ylist, 1)
        return zip(list(_x), list(_y))

# -----

def click(image,hdr):
    """Plot given image array and collect clicked coords

    Input:
     - image  <ndarray> : Image data
     - hdr    <header>  : Header object

    Output:
    - {'x','y','ra','dec'}  <dict> : selected points
    
    """

    imshow(image);
    clk = CoordPoints();
    cid = connect('button_press_event',clk);
    show();
    disconnect(cid);

    x, y = zip(*clk.xy_coords());
    if (hdr):
        ra, dec = zip(*clk.dg_coords(hdr));
    else:
        ra = dec = None;


    return({'x':x, 'y':y, 'ra':ra, 'dec':dec});

# ---

def select_objects(segimg, objimg, centroids, header=None):
    """Select objects based on given positions

    Input:
     - segimg  <ndarray> :
       Numpy array with SEGMENTATION image information
     - objimg  <ndarray> :
       Numpy array with OBJECTS/original image information
     - centroids  <list:(int,int)> :
       List of tuples with positions to search for objects [(x0,y0),(x1,y1),...]

    """

    # Initialize object IDs list..
    #
    objIDs = [];
    
    # Is there an identified object on given (xo,yo) point? If yes, store its ID..
    #
    for o_o in centroids:
        xo, yo = o_o;
        objid = segimg[yo,xo];
        if ( objid != 0 ):
            objIDs.append(objid);
        else:
            print >> sys.stdout, "No objects were identified for (x=%s,y=%s) position." % (xo,yo);

    # Read out the objects recognized from images..
    #
    objs, hdrs = imcp.sextamp( segimg, objimg, header, increase=2, relative_increase=True, objIDs=objIDs );

#    image_out, hdr = cutout( obj_img, header=hdr, xo=int(xo), yo=int(yo), x_size=int(x_size), y_size=int(y_size), mask=ind );


    return ({'IDs' : objIDs, 'images' : objs, 'headers' : hdrs});

# ---

def write_results(images=[], headers=[], IDs=[], out_rootname='out_'):

    for i in range( len(IDs) ):

        outname = out_rootname+"%03d.fits" % (IDs[i]);
        print >> sys.stdout, "Writing object (id: %s) image:  %s ..." % (IDs[i],outname);
        os.system( 'rm %s &> /dev/null' % (outname) );
        pyfits.writeto( outname, images[i], headers[i] );

# ---

def generate_catalogues(fits_image, objects=[], headers=[], objIDs=[], cat_dir=''):
    """Generates catalogues for poststamps from given image"""
    import string;

    def csv(cat_name, *keys, **key_lists):
        """Write a CSV catalog"""
        catFile = open(cat_name,'w');
        catObj = csv.writer(catFile, delimiter=',', quotechar='\"');


    # Check basic stuff for catalogues creation..
    #
    if (objects==[]):
        return None;
    if (len(objects) != len(headers)  or  len(objects) != len(objIDs)):
        return False;

    img_rootname = string.join( string.split(fits_image, sep=".")[:-1], sep="." );
    
    # Set the output directory based on given image name if not('') was given..
    #
    if (cat_dir==''):
        cat_dir = img_rootname;
    try:
        os.system('mkdir %s' % (cat_dir))
    except:
        pass;

    cat_name = "Catobjs_"+img_rootname+".dat";
#    fp = open(cat_name,'w');
#    fp.write("ID,filename,source_image,x,y,ra,dec")
#    for _i in range(len(objects)):
#        _name = "obj_%03d" % ();
#        fp.write("%s," % (objIDs[_i],))
    
# ---

def run(fits_image, preset='', use_header=True):
    """Present image for interactive choice of the object to segment

    dic = run( fits_image )

    Given FITS image is presented on an interacive window where the
    user can click on objects for selection. A segmentation algorithm
    will run on the image and selected points will be cutout if they
    correspond to a segmented object. 'preset' argument can be used
    to use pre-set Sextractor(SE) configuration (see sltools.Package.sextractor
    for more information). If we don't want to use image's header info
    'use_header' argument can be set to "False"; default is "True"

    Input:
     - fits_image   <str> : FITS file containg the image to use
     - preset       <str> : SE's pre-set configuration
     - use_header  <bool> : whether to use (or not) image's header

    Output:
     - {'IDs', 'images', 'headers'}
     
    """

    # Open given image..
    #
    if (use_header):
        img, header = pyfits.getdata(fits_image, header=use_header);
    else:
        img = pyfits.getdata(fits_image, header=use_header);
        header = None;

    # Show the image and make "clickable" for objects selection..
    #
    _dic = click(img, header);
    x = _dic['x'];
    y = _dic['y'];
    ra = _dic['ra'];
    dec = _dic['dec'];
    del _dic;
    
    centroids = zip(x,y);
    
    # Segment image with Sextractor..
    #
    if preset=='none':
        preset='';
    _dic = SE.run_segobj(fits_image, preset=preset);
    if not _dic:
        return False;

    objimg = pyfits.getdata( _dic['OBJECTS'] );
    segimg = pyfits.getdata( _dic['SEGMENTATION'] );
    cat = pyfits.open( _dic['CATALOG'] )[1].data;
    del _dic;


    _dic = select_objects(segimg, objimg, centroids, header);
    if not _dic:
        return False;
    
    return (_dic);

# ---

if __name__ == "__main__" :

    if (len(sys.argv) < 3):
        print "Usage: \n   %s  <fits_image_file>  {HST,DC4,DC5,CFHT,none}" % (sys.argv[0]);
        sys.exit(1);

    out = run(sys.argv[1],sys.argv[2]);

    write_results(images=out['images'], headers=out['headers'], IDs=out['IDs']);
    
    print >> sys.stdout, "Done.";
