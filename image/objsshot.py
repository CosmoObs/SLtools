#!/usr/bin/env python

"""Program for automate SExtractor run and produce objects postamps"""

##@package objsshot
#
# Objects in given image are identified using SExtractor and poststamps
#are created for each of the (segmented) objects.
#
# Executable: YES
#
# To run the code, just type:
#
# > python objsshot.py 
#
# Use "--help" for a list of program options


import sys;
import os;

import gc;
gc.enable();

import optparse;
from optparse import OptionParser;

import re;
import string;
import numpy as np;
import pyfits;

from sltools.io import *;


#---------------------
def _valid_line(wort):

    wort.rstrip("\n");
    wort = re.sub("^\s*","",wort);
    wort = re.sub("\s+"," ",wort);
    wort = re.sub("#.*","",wort);

    return ( not (bool(re.search("^$",wort))  or  bool(re.search("^#",wort))) );

# ---

#---------------------------
def _check_config(cfg_file):
    """Check if given config file has the necessary structure.

    In particular, if checking fails, a fixing procedure based on
    Sextractor config file (default.sex) design.
    """

    fp = open(cfg_file, 'r');

    file = filter( _valid_line, fp.readlines() );

    for line in file:

        if ( re.search("\[.*\]",line) ) :
            continue;

        if not ( bool(re.search("=",line))  or  bool(re.search(":",line)) ) :
            return (False);

        if not ( len(string.split(line,sep="=")) >= 2  or  len(string.split(line,sep=":")) >= 2 ) :
            return (False);

    fp.close();


    return (True);

# ---

# ===========================================================================
def run(fits_image, params=[], args={}, header=False, increase=0, custom=''):
    """
    Function to run SExtractor for objects segmentation and poststamp creation

    objsshot.run( image.fits [,...] )

    Function receives a FITS image, and possibly a list of SEx output params,
    a dictionary with SEx command-line(/config) arguments. It can be chosen to
    update header (sky coordinates info).

    'increase' is used for dimensioning the output poststamp of identified
    objects. It is a multiplicative factor for resizing poststamp relative to
    object dimensions. (Notice that 0 < increase < 1 is possible and will
    return a poststamp with part of the object)

    Custom SExtractor configuration parameters can be used, chosen to be *good*
    extraction parameters for different telescopes.
    custom : 'HST'
             'DC4'
             'DC5'
    Notice that this custom choices do not disable other parameters.

    Input:
     - fits_image : FITs image file containing objects to be segmented
     - params     : list of parameters to output on sextractor catalogues
     - config     : dictionary with [default] section contents.
     - header     : Whether or not to use (and update) header info.
     - increase   : multiplicative factor for poststamp sizing
     - custom     : HST, DC4, DC5

    Output:
     - (dic) : dictionary structure with the keys 'objIDs', 'images' and 'headers'.
               Each of them is a list with corresponding information.

    """


    # Now we set some Sextractor arguments for this script purposes..
    #
    rootname = string.split( string.replace( fits_image,".fits","" ), sep="/" )[-1];

    objname = rootname+'_obj.fits';
    segname = rootname+'_seg.fits';
    catalog = rootname+'_cat.fit';

    # And update given arguments (via 'args.dic') with necessary parameters..
    #
    args.update( {'checkimage_type' : 'OBJECTS,SEGMENTATION', 'checkimage_name' : objname+','+segname, 'catalog_type' : 'FITS_1.0', 'catalog_name' : catalog} );


    # Run sextractor module; output catalogues to "catalog" files..
    #
    out = segment.run_sex(fits_image, params, args, custom);

    # Check if Sextractor ran OK. It is suppose to output an dictionary with fits filenames if all OK,
    # otherwise, if not, return False:
    #
    if (out == False):
        print >> sys.stderr, "Error: Sextractor raised and error coded during segmentation. Finishing run."
        return (False);


    # Read fresh images just created by sex run..
    #
    obj = pyfits.getdata( objname );
    seg = pyfits.getdata( segname );


    # If 'header' argument not False/None, read it..
    #
    if (header):
        header = pyfits.getheader( fits_image );
    else:
        header = None;


    # Open "segmentation" output catalog from SE. Notice that object IDs are supposed to lie on 1st column..
    #
# ASCII version: CATALOG_TYPE = ASCII_HEAD
#    cat = np.loadtxt( catalog, comments="#" );
#    try:
#        objIDs = list( cat[:,0].astype(np.int32) );
#    except:
#        objIDs = list( cat.astype(np.int32) );
#
# FITS version: CATALOG_TYPE = FITS_1.0
    cat = pyfits.open(catalog)[1].data;
    objIDs = cat.field('NUMBER');


    ## objects_IDfy return a list of arrays; result from the objIDs objects located in seg image merged with obj image.
    # So, each object in objIDs is returned (also inside a/the list) as array with their corresponding pixels..
    #
    objs, hdrs = imcp.sextamp( seg, obj, header, increase=increase, relative_increase=True, objIDs=objIDs );


    return ({'IDs' : objIDs, 'images' : objs, 'headers' : hdrs});

# ---

# \cond
# =======================================================================================================

if __name__ == "__main__" :


    # Checking dependencies..
    stat = check_dependencies( { 'binaries':['sex'] } );
    if not stat:
        print >> sys.stderr, "Error: Dependencies check failed. Fix it!";
        sys.exit(2);

    # Read input options..
    usage="Usage:  %prog [options] <objects_image.fits>";
    parser = OptionParser(usage=usage);

    parser.add_option('-p', '--params',
                      dest='params', default='default.param',
                      help="File with parameters to output [default.param]");
    parser.add_option('--config',
                      dest='config', default=None,
                      help="Config file for (S)extraction");
    parser.add_option('--header', action='store_true',
                      dest='use_header', default=False,
                      help="Use image header for WCS adjustment(?)");
    parser.add_option('--increase',
                      dest='incr', default=0,
                      help="Multiplicative factor (relative to image size) for image border");
    parser.add_option('--instrument',
                      dest='telesc', default='',
                      help="DC$/DC5/HST");
    parser.add_option('--outdir',
                      dest='dir_name', default='out_stamps',
                      help="Directory name to put output images");
    parser.add_option('--outname',
                      dest='file_name', default='pstamp_',
                      help="Name component for output images, before '#id.fits'");

    (opts,args) = parser.parse_args();

    if ( len(args) != 1 ):
        parser.error("Wrong number of arguments. Try option '--help' for usage and help messages.");

    fits_image = args[0];
    configfile = opts.config;
    paramsfile = opts.params;
    use_header = opts.use_header;
    increase   = opts.incr;
    dir_name   = opts.dir_name;
    out_name   = opts.file_name;
    instrument = opts.telesc;


    # Check if params file was given..
    #
    if ( paramsfile == None ):
        parser.error("No params file was given. Try option '--help' for usage and help messages.\n");

    # Read parameters from given file..
    #
    params = open( paramsfile ).readlines();
    params = regexp.line_filter( params, comments="#" );


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


    # ===================================================================================


    # OK.. Now lets run the objects identification procedure..
    # I'll first create a directory to put the output images in there
    #
    try:
        os.mkdir( dir_name );
    except OSError:
        pass;

    os.chdir( dir_name );


    # Run ObjSShot main function. It passes fits filename and read params and config settings.
    # The output is a list of image arrays (numpy) with identified objects on them.
    #
    outdic = run( "../"+fits_image, params, config['default'], header=use_header, increase=float(increase), custom=instrument );

    objIDs = outdic['IDs'];
    objs_list = outdic['images'];
    hdrs_list = outdic['headers'];


    # And for each object array in list, write down their pixels ;-)
    #
    for i in xrange( len(objIDs) ):

        outname = out_name+"%03d.fits" % (i);
        os.system( 'rm %s &> /dev/null' % (outname) );
        pyfits.writeto( outname, objs_list[i], hdrs_list[i] );


    sys.exit(0);

# \endcond
