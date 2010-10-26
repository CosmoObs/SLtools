#!/usr/bin/env python

import sys;
import os;
import pyfits;
from pylab import *;
from pywcs import WCS;

# SLTOOLS:
from sltools.image import imcp;
from sltools.Packages import sextractor as SE;


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

def fits_2_img(fits_image, preset=''):
    """Run Sextractor and read images to arrays"""

    if preset=='none':
        preset='';
    _dic = SE.run_segobj(fits_image, preset=preset);

    if not _dic:
        return False;
    
    objname = _dic['OBJECTS_file'];
    segname = _dic['SEGMENTATION_file'];
    catalog = _dic['CATALOG_file'];
    
    # Read FITS images and catalog..
    #
    objimg = pyfits.getdata( objname );
    segimg = pyfits.getdata( segname );
    cat = pyfits.open(catalog)[1].data;


    return (segimg,objimg,cat);

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
            print >> sys.stdout, "No objects were identified on given (x=%s,y=%s) position." % (xo,yo);
            return None;

    # Read out the objects recognized from images..
    #
    objs, hdrs = imcp.sextamp( segimg, objimg, header, increase=2, relative_increase=True, objIDs=objIDs );

#    image_out, hdr = cutout( obj_img, header=hdr, xo=int(xo), yo=int(yo), x_size=int(x_size), y_size=int(y_size), mask=ind );


    return ({'IDs' : objIDs, 'images' : objs, 'headers' : hdrs});

# ---

def run(fits_image, preset='', use_header=True):

    # Open given image..
    #
    if (use_header):
        img, header = pyfits.getdata(fits_image, header=use_header);
    else:
        img = pyfits.getdata(fits_image, header=use_header);
        header = None;

    # Show the image and make "clickable" for objects selection..
    #
    dic = click(img, header);
    x = dic['x'];
    y = dic['y'];
    ra = dic['ra'];
    dec = dic['dec'];
    
    print "%s:" % (fits_image);
    print "%s  X     :  Y     :  RA       :  DEC" % (len(fits_image)*' ');
    if (header):
        for i in range(len(x)):
            print "%s %4.2f : %4.2f : %3.5f : %2.5f" % (len(fits_image)*' ',x[i],y[i],ra[i],dec[i]);
    else:
        for i in range(len(x)):
            print "%s %.2f : %.2f : --- : ---" % (len(fits_image)*' ',x[i],y[i]);

    # Segment image? : Yes
    #
    _out = fits_2_img(fits_image, preset);
    if not _out:
        return False;
    segimg, objimg, cat = _out[0], _out[1], _out[2];

    centroids = zip(x,y);
    
    dic = select_objects(segimg, objimg, centroids, header);
    if not dic:
        return False;
    objIDs = dic['IDs'];
    objs = dic['images'];
    hdrs = dic['headers'];
    
    imshow(objs[0]);
    show()
    
    for i in range( len(objIDs) ):

        outname = "out_segmented_%02d.fits" % (objIDs[i]);
        print >> sys.stdout, "Writing object (id: %s) image:  %s ..." % (objIDs[i],outname);
        os.system( 'rm %s &> /dev/null' % (outname) );
        pyfits.writeto( outname, objs[0], hdrs[0] );
#        pyfits.writeto( outname, objs[i] );

    print >> sys.stdout, "Done.";

    return (dic);

# ---

if __name__ == "__main__" :

    if (len(sys.argv) < 3):
        print "Usage: \n   %s  <fits_image_file>  {HST,DC4,DC5,CFHT,none}" % (sys.argv[0]);
        sys.exit(1);

    run(sys.argv[1],sys.argv[2]);
