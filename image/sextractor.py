#!/usr/bin/env python

"""Module to deal with Sextractor"""

##@ sextractor
#
#
# This package runs Sextractor.
#
# Executable package : No


import sys;

import os;
import re;
import string;
import pyfits;
import numpy as np;
import sltools;


# PRESET CONFIG FOR SEXTRACTOR/INSTRUMENTS #
# ------------------------------------------
def _cstm_args(instrument):

    if instrument == 'DC4':
        _dic = {
            'FILTER_NAME' : 'default.conv',
            'STARNNW_NAME' : 'default.nnw',
            'DETECT_TYPE' : 'CCD',
            'DETECT_IMAGE' : 'SAME',
            'FLAG_IMAGE' : 'NONE',
            'DETECT_MINAREA' : '8',
            'THRESH_TYPE' : 'RELATIVE',
            'DETECT_THRESH' : '2.5',
            'ANALYSIS_THRESH' : '2.5',
            'FILTER' : 'Y',
            'DEBLEND_NTHRESH' : '20',
            'DEBLEND_MINCONT' : '0.0005',
            'CLEAN' : 'Y',
            'CLEAN_PARAM' : '1.0',
            'BLANK' : 'Y',
            'PHOT_APERTURES' : '5',
            'PHOT_AUTOPARAMS' : '2.5,',
            'SATUR_LEVEL' : '50000.0',
            'MAG_ZEROPOINT' : '31.0414',
            'MAG_GAMMA' : '4.0',
            'GAIN' : '0.0',
            'PIXEL_SCALE' : '0.27',
            'SEEING_FWHM' : '1.2',
            'BACK_SIZE' : '64',
            'BACK_FILTERSIZE' : '3',
            'BACKPHOTO_TYPE' : 'GLOBAL',
            'BACKPHOTO_THICK' : '24',
            'MEMORY_OBJSTACK' : '2000',
            'MEMORY_PIXSTACK' : '100000',
            'MEMORY_BUFSIZE' : '512',
            'SCAN_ISOAPRATIO' : '0.6',
            'VERBOSE_TYPE' : 'NORMAL'
            };

    elif instrument == 'DC5':
        _dic = {
            'FILTER_NAME' : 'default.conv',
            'STARNNW_NAME' : 'default.nnw',
            'DETECT_TYPE' : 'CCD',
            'DETECT_IMAGE' : 'SAME',
            'FLAG_IMAGE' : 'NONE',
            'DETECT_MINAREA' : '10',
            'THRESH_TYPE' : 'RELATIVE',
            'DETECT_THRESH' : '2.5',
            'ANALYSIS_THRESH' : '2.5',
            'FILTER' : 'Y',
            'DEBLEND_NTHRESH' : '20',
            'DEBLEND_MINCONT' : '0.0005',
            'CLEAN' : 'Y',
            'CLEAN_PARAM' : '1.0',
            'BLANK' : 'Y',
            'PHOT_APERTURES' : '5',
            'PHOT_AUTOPARAMS' : '2.5,3.5',
            'SATUR_LEVEL' : '50000.0',
            'MAG_ZEROPOINT' : '31.0414',
            'MAG_GAMMA' : '4.0',
            'GAIN' : '0.0',
            'PIXEL_SCALE' : '0.27',
            'SEEING_FWHM' : '0.8',
            'BACK_SIZE' : '64',
            'BACK_FILTERSIZE' : '3',
            'BACKPHOTO_TYPE' : 'GLOBAL',
            'BACKPHOTO_THICK' : '24',
            'MEMORY_OBJSTACK' : '2000',
            'MEMORY_PIXSTACK' : '100000',
            'MEMORY_BUFSIZE' : '512',
            'SCAN_ISOAPRATIO' : '0.6',
            'VERBOSE_TYPE' : 'NORMAL'
            };

    elif instrument == 'HST':
        _dic = {
            'FILTER_NAME' : 'default.conv',
            'STARNNW_NAME' : 'default.nnw',
            'DETECT_TYPE' : 'CCD',
            'DETECT_IMAGE' : 'SAME',
            'FLAG_IMAGE' : 'NONE',
            'DETECT_MINAREA' : '10',
            'THRESH_TYPE' : 'RELATIVE',
            'DETECT_THRESH' : '2.7',
            'ANALYSIS_THRESH' : '2.7',
            'FILTER' : 'Y',
            'DEBLEND_NTHRESH' : '20',
            'DEBLEND_MINCONT' : '0.0005',
            'CLEAN' : 'Y',
            'CLEAN_PARAM' : '1.0',
            'BLANK' : 'Y',
            'PHOT_APERTURES' : '5',
            'PHOT_AUTOPARAMS' : '2.5,3.5',
            'SATUR_LEVEL' : '50000.0',
            'MAG_ZEROPOINT' : '21.59',
            'MAG_GAMMA' : '4.0',
            'GAIN' : '7.0',
            'PIXEL_SCALE' : '0.1',
            'SEEING_FWHM' : '0.22',
            'BACK_SIZE' : '64',
            'BACK_FILTERSIZE' : '3',
            'BACKPHOTO_TYPE' : 'GLOBAL',
            'BACKPHOTO_THICK' : '24',
            'MEMORY_OBJSTACK' : '3000',
            'MEMORY_PIXSTACK' : '300000',
            'MEMORY_BUFSIZE' : '1024',
            'SCAN_ISOAPRATIO' : '0.6',
            'VERBOSE_TYPE' : 'NORMAL'
            };

    elif instrument == 'CFHT':
        _dic = {
            'DETECT_TYPE' : 'CCD',
            'DETECT_MINAREA' : '3',
            'DETECT_THRESH' :  '3',
            'ANALYSIS_THRESH' : '3',
            'FILTER' : 'Y',
            'FILTER_NAME' : 'default.conv',
            'DEBLEND_NTHRESH' : '32',
            'DEBLEND_MINCONT' : '0.005',
            'CLEAN' : 'Y',
            'CLEAN_PARAM' : '1.0',
            'MASK_TYPE' : 'CORRECT',
            'PHOT_APERTURES' : '10',
            'PHOT_AUTOPARAMS' : '2.5,3.5',
            'SATUR_LEVEL' : '50000.0',
            'MAG_GAMMA' : '4.0',
            'GAIN' : '0.0',
            'PIXEL_SCALE' : '1.0',
            'SEEING_FWHM' : '1.2',
            'STARNNW_NAME' : 'default.nnw',
            'BACK_SIZE' : '64',
            'BACK_FILTERSIZE' : '3',
            'BACKPHOTO_TYPE' : 'GLOBAL',
            'CHECKIMAGE_TYPE' : 'NONE',
            'MEMORY_OBJSTACK' : '3000',
            'MEMORY_PIXSTACK' : '300000',
            'MEMORY_BUFSIZE' : '1024',
            'VERBOSE_TYPE' : 'NORMAL'
            };

    else:
        _dic = {};


    return (_dic);

# ---

#=======================================================================
def run(fits_image, params=[], args={}, preset=''):
    """SExtract the Image and read the outputs to array

    Input:
    - fits_image <str> : filename
    - use_header <bool> : True/False
    - params <list str> : See SE's default.param file. Those params
    - args <dic> : SE's command-line arguments
    - preset <str> : preset configs (DC4,DC5,...)

    """


    # Object features to output..
    #
    if ( params != [] ):

        param_file = "run_sex.param";
        args.update( {'parameters_name' : param_file} );

        fparam = open(param_file,'w');
        for _param in params:
            fparam.write( "%s\n" % (_param) );
        fparam.close();


    # ==========================================================================
    # WRITE DOWN "default.conv" for sextractor if necessary..
    #
    if (not args.has_key('filter_name')  or  args['filter_name']=='default.conv'):
        fp = open('default.conv','w');
        fp.write("CONV NORM\n");
        fp.write("# 3x3 ``all-ground'' convolution mask with FWHM = 2 pixels.\n");
        fp.write("1 2 1\n");
        fp.write("2 4 2\n");
        fp.write("1 2 1\n");
        fp.close();
    # ==========================================================================


    cargs = _cstm_args( preset );
    cargs.update( args );

    # Build up sextractor command line..
    #
    cmd_line = '';
    for key in cargs:
        cmd_line = cmd_line + ' -' + string.upper(key) + ' ' + re.sub( "\s+","",str(cargs[key]));

    # Run sex..
    #
    status = os.system( 'sex %s %s' % (fits_image,cmd_line));
    if ( status != 0 ):
        print >> sys.stderr, "Error: Sextractor raised and error code '%s' during the run." % (status);


    if (status):
        return (False);

    return (True);

# -

#---------------------------------------------------------
def run_segobj(fits_image, params=[], args={}, preset=''):
    """
    Run Sextractor on given (fits) image, with optional parameters

    run( 'image.fits' [ [...],{...} ] )

    Function runs Sextractor on a given FITS image (fits_image) to generate 
    sex's "OBJECTS" and "SEGMENTATION" outputs. Output parameters (catalogs) 
    can be passed in a list of strings. SE command-line arguments can
    be passed as a dictionary {'key' = 'value'}.

    Two catalogue files are created during sex's runs: cats_name+'_obj.cat' and 
    cats_name+'_seg.cat', with "objects" and "segmentation" features output 
    values, respectively.
    Two FITS files - images - are also created during the runs, one for each
    run type: "objects" and "segmentation". File names are just the given 
    filename with extension ".fits" changed by "_obj.fits" and "_seg.fits", 
    respectively.

    Custom SExtractor configuration parameters can be used, chosen to be *good*
    extraction parameters for different telescopes.
    preset : 'HST'
             'DC4'
             'DC5'
             'CFHT'
    Notice that this custom choices do not disable other parameters.

    Input:
     - fits_image : FITS image filename
     - params     : List of parameters to output on catalogues
     - args       : Dictionary with sextractor line arguments. If not given,
                    SE defaults will be used.
     - preset     : 'HST', 'DC4', 'DC5' or 'CFHT'

    Output:
    - {'OBJECTS':str, 'SEGMENTATION':str, CATATLOG:str}
        SE's OBJECTS, SEGMENTATION, CATALOG output filenames
      
    """


    # Now we set some Sextractor arguments for this script purposes..
    #
    _rootname = string.split( string.replace( fits_image,".fits","" ), sep="/" )[-1];
    objimgname = _rootname+'_obj.fits';
    segimgname = _rootname+'_seg.fits';
    catfilename = _rootname+'_cat.fit';

    # Update given parameters, for output; we need at least 'NUMBER'..
    #
    if not ('NUMBER' in params):
        params.append('NUMBER');
        
    # And update given arguments (via 'args.dic') with necessary parameters..
    #
    args.update( {'checkimage_type' : 'OBJECTS,SEGMENTATION', 'checkimage_name' : objimgname+','+segimgname, 'catalog_type' : 'FITS_1.0', 'catalog_name' : catfilename} );

    # Run sextractor module; output catalogues to "catalog" files..
    #
    out = run(fits_image, params, args, preset);

    if (out == False):
        return (False);

    return ({'OBJECTS':objimgname, 'SEGMENTATION':segimgname, 'CATALOG':catfilename});

# ---

# ================================================================================
def _valid_line(wort):
    wort.rstrip("\n");
    wort = re.sub("^\s*","",wort);
    wort = re.sub("\s+"," ",wort);
    wort = re.sub("#.*","",wort);
    return ( not (bool(re.search("^$",wort))  or  bool(re.search("^#",wort))) );

def _check_config(cfg_file):
    """Check if given config file has the necessary structure.

    In particular, if checking fails, a fixing procedure based on
    Sextractor config file (default.sex) design.
    """

    fp = open(cfg_file, 'r');
    file = filter( _valid_line, fp.readlines() );

    cnt_fail = 0;
    cnt_succ = 0;
    for line in file:
        if ( bool(re.search("\[.*\]",line)) ) :
            continue;
        else:
            cnt_fail += 1;
        if not ( bool(re.search("=",line))  or  bool(re.search(":",line)) ) :
            cnt_fail += 1;
        if not ( len(string.split(line,sep="=")) >= 2  or  len(string.split(line,sep=":")) >= 2 ) :
            cnt_fail += 1;
        cnt_succ += 1;

    fp.close();

    if (bool(cnt_fail)  and  not bool(cnt_succ)):
        return (False);
    elif (bool(cnt_succ)  and  not bool(cnt_fail)):
        return (True);
    else:
        return (None);

def _fix_config(config_file):

    fp = open(config_file,'r');
    fout = config_file+'.mod';
    fpfix = open(fout,'w');

    print >> fpfix,"[default]";

    file = filter( _valid_line,fp.readlines() );
    for line in file:
        laine = re.sub("\s","=");
        laine = re.sub("\s",":");
        print >> fpfix,"%s" % (laine);

    return (fout);

# ---

# \cond
#==========================
if __name__ == "__main__" :

    from sltools.io import config_parser as CP;
    
    from optparse import OptionParser;

    usage="Usage: %prog [options] <image_fits_file>"
    parser = OptionParser(usage=usage);

    parser.add_option('-p', '--params_file',
                        dest='params', default=None,
                        help="Sextractor params file");
    parser.add_option('-c', '--config_file',
                        dest='config', default=None,
                        help="Sextractor config file");
    parser.add_option('--segment', action='store_true',
                        dest='run_seg', default=False,
			  help="Run sextractor in SEGMENTATION mode");
    parser.add_option('--preset',
                        dest='preset', default='',
                        help="Use preset SE's configuration for different instruments: DC4, DC5, HST, CFHT");
        
    (opts,args) = parser.parse_args();

    run_seg = opts.run_seg;
    params_file = opts.params;
    config_file = opts.config;
    preset = opts.preset;

    if ( len(args) != 1 ):
        parser.error("Wrong number of arguments. Try '--help' for usage and help messages.");

    # FITS image file:
    #
	infits = args[0];

    # SE's config file:
    #
    if (config_file == None):
        print >> sys.stdout, "Warning: No SE's config file given. Using program defaults.";
        config = {};
    else:
        _status = _check_config(config_file);
        if ( _status == False ):
            config_file = _fix_config(config_file);
        elif ( _status == None ):
            print >> sys.stderr, "Error: Given config file does not looks like a INI file neither a plain key-value one.";
            sys.exit(1);
        else:
            config = CP.read_ini(config_file);

    # SE's param file:
    #
    if (params_file):
        fp_params = open(params_file);
        params = [ line.rstrip("\n")   for line in fp_params.readlines() ]
    else:
        params = [];

    # Run SE:
    #
    if ( run_seg ):
        _out = run_segobj(infits, params=params, args=config, preset=preset);
        if (_out):
            print "Segmentation output files:";
            print "CATALOG: %s" % (_out['CATALOG']);
            print "OBJECTS: %s" % (_out['OBJECTS']);
            print "SEGMENT: %s" % (_out['SEGMENTATION']);
        else:
            sys.exit(1);
    else:
        _out = run(infits, params=params, args=config, preset=preset);
        if not (_out):
            sys.exit(1);

    sys.exit(0);

# ------------------
# \endcond

"""
    # Read fresh images just created by sex run..
    #
    objimg = pyfits.getdata( objimgname );
    segimg = pyfits.getdata( segimgname );

    # FITS version: CATALOG_TYPE = FITS_1.0
    cat = pyfits.open(catalog)[1].data;

    # If 'header' argument not False/None, read it..
    #
    if (use_header):
        header = pyfits.getheader( fits_image );
    else:
        header = None;


    return (objimg, segimg, header, cat);
"""
