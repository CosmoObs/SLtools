#!/usr/bin/env python

import sys;
import os;

import re;
import string;
import numpy as np;
import pyfits;

import sltools;
from sltools.image import sextractor,imcp;
from sltools.catalog. import fits_data as fts;

# ---

def correlate_SEtables(filenames,centroids):
    """
    Look for SE features for each centroid on seg. images
    
    A list of segmented image filenames is read from 'filenames'
    together with the list of 'centroids' (x,y). The (x,y) information
    is used to read the corresponding (segmented) object on seg. image
    and then - with the ID just read -, the particular object features
    is read from SE' catalog.
    
    So far, the catalog filename corresponding to the image filename is
    recovered substituting the ".fits" image filename extension with 
    "_cat.fit" extension.
    
    Input:
     - filenames : [str,]
        List of original image filenames
     - centroids : [(int,int),]
        List of objects centroids
    
    Output:
     - tbhdu : pyfits.BinTableHDU
        Table (HDU) with object features
    
    ---
    """
    imgname = '';
    _tbList = [];
    # Detect objects corresponding to given centroids..
    #
    for i in xrange(len(centroids)):
        if filenames[i] != imgname:
            imgname = filenames[i];
            segname = re.sub(".fits","_seg.fits",imgname);
            catname = re.sub(".fits","_cat.fit",imgname);
            segimg = pyfits.getdata(segname);
            tbhdu_sex = pyfits.open(catname)[1];
            
        objID = segobjs.centroid2ID(segimg,centroids[i]);
        _tbList.append(fts.select_entries(tbhdu_sex,'NUMBER',objID));

    new_tbhdu = _tbList[0];
    for i in xrange(1,len(_tbList)):
        new_tbhdu = extend_tbHDU(new_tbhdu,_tbList[i]);

    return new_tbhdu;
    
# ---

def tabled_images(tbhdu_tbm, x='X', y='Y', filename='filename'):
    """ Read objects from inside the given table HDU
    
    Already segmented, and object images filenames given by column
    'filename', 'x' and 'y' values will be used to search for the
    corresponding object pixels and ID in "segmentation" image.
    """

    filenames = tbhdu_tbm.data.field(filename);
    centroids = zip(tbhdu_tbm.data.field(x),tbhdu_tbm.data.field(y));

        # Access the objID (only) image,
#        _otmp, _htmp = imcp.segstamp( segimg, objID=objID, hdr=hdr, objimg=None, increase=4, relative_increase=True );

    return
    
#    return {'IDs' : objIDs, 'images' : objs, 'headers' : hdrs};


# \cond
# =======================================================================================================

if __name__ == "__main__" :

    import optparse;
    from optparse import OptionParser;

    # Read input options..
    usage="Usage:  %prog [options] <objects_image.fits>";
    parser = OptionParser(usage=usage);

    parser.add_option('-p',
                      dest='params', default='default.param',
                      help="File with parameters to output (SE's default.param)");
    parser.add_option('-c',
                      dest='config', default=None,
                      help="Config file for (S)extraction");
    parser.add_option('-o',
                      dest='file_name', default='pstamp_',
                      help="Name component for output images, before '#id.fits'");
    parser.add_option('--increase',
                      dest='incr', default=0,
                      help="Multiplicative factor (relative to object size) for image border");
    parser.add_option('--no-header', action='store_false',
                      dest='use_header', default=True,
                      help="Do not use existing image header.");
    parser.add_option('--xo',
                      dest='xo', default=0,
                      help="Postion (X) to search for an identified object");
    parser.add_option('--yo',
                      dest='yo', default=0,
                      help="Postion (Y) to search for an identified object");
    parser.add_option('--preset',
                      dest='preset', default='',
                      help="DC4/DC5/HST/CFHT");
    parser.add_option('--outdir',
                      dest='dir_name', default='out_objsshot',
                      help="Directory name to put output images");

    (opts,args) = parser.parse_args();

    if ( len(args) != 1 ):
        parser.error("Wrong number of arguments. Try option '--help' for usage and help messages.");

    fits_image = args[0];
    configfile = opts.config;
    paramsfile = opts.params;
    use_header = opts.use_header;
    increase   = opts.incr;
    xc = opts.xo;
    yc = opts.yo;
    dir_name   = opts.dir_name;
    out_name   = opts.file_name;
    instrument = opts.preset;


    # Check if config file and then ewad it..
    #
    if ( configfile != None ):

        # Check if config file is valid:
        #
        stat = _check_config(configfile)
        if not (stat):
            print >> sys.stderr, "\nError: Given config file, %s, seems not to be a valid one." % (configfile);
            print >> sys.stderr, "        Check if it has a [default] section and valid <key : value> lines.\n";
            sys.exit(2);

        # Try to read it:
        #
        try:
            config = cp.read_config(configfile, 'default');
            if ( config == {} ):
                raise ValueError;
        except:
            print >> sys.stderr, "Error: Check your config file. Remember to place SE's options below [default] section.";
            sys.exit(2);

    else:
        config = {'default' : {}};


    # Read parameters from given file..
    #
    try:
        params = open( paramsfile ).readlines();
        params = regexp.line_filter( params, comments="#" );
    except:
        print >> sys.stderr, "Error: No SE's parameters file (.param). Use -p option to pass the filename.";
        sys.exit(1);

    # ===================================================================================


    # OK.. Now lets run the objects identification procedure..
    # I'll first create a directory to put the output images in there
    #
    try:
        os.mkdir( dir_name );
    except OSError:
        pass;

    owd = os.getcwd()+"/";
    os.chdir( dir_name );


    # Run ObjSShot main function. It passes fits filename and read params and config settings.
    # The output is a list of image arrays (numpy) with identified objects on them.
    #
    if ( xc != 0 and yc != 0 ):
        outdic = pickout_obj( owd+fits_image, xo=xc, yo=yc, params=params, args=config['default'], header=use_header, increase=float(increase), preset=instrument );
    else:
        outdic = readout_objs( owd+fits_image, params, config['default'], header=use_header, increase=float(increase), preset=instrument );

    if ( outdic == False ):
        print >> sys.stderr,"Finished.";
        sys.exit(1);
    if ( outdic == None ):
        print >> sys.stdout,"Finished.";
        sys.exit(0);

    objIDs = outdic['IDs'];
    objs_list = outdic['images'];
    hdrs_list = outdic['headers'];


    # And for each object array in list, write down their pixels ;-)
    #
    for i in xrange( len(objIDs) ):

        print >> sys.stdout, "Writing object (id:) %s image..." % (objIDs[i]);
        outname = out_name+"%04d.fits" % (objIDs[i]);
        os.system( 'rm %s &> /dev/null' % (outname) );
        pyfits.writeto( outname, objs_list[i], hdrs_list[i] );


    print >> sys.stdout, "Done.";
    os.chdir( owd );
    sys.exit(0);

# \endcond
