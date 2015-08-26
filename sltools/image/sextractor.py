#!/usr/bin/env python

"""Module to deal with Sextractor"""

##@ sextractor
#
#
# This package runs Sextractor.
#
# Executable package : Yes


import sys;

import logging;
import os;
import re;
import string;
import pyfits;
import numpy as np;
import sltools;


# PRESET CONFIG FOR SEXTRACTOR/INSTRUMENTS #
# ------------------------------------------
def se_presets(instrument):

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
            'PHOT_AUTOPARAMS' : '2.5,3.5',
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
    
    elif instrument == 'HST_Arcs':
       _dic = {
            'FILTER_NAME' : 'default.conv',
            'STARNNW_NAME' : 'default.nnw',
            'DETECT_TYPE' : 'CCD',
            'DETECT_IMAGE' : 'SAME',
            'FLAG_IMAGE' : 'NONE',
            'DETECT_MINAREA' : '20',
            'THRESH_TYPE' : 'RELATIVE',
            'DETECT_THRESH' : '2.7',
            'ANALYSIS_THRESH' : '2.7',
            'FILTER' : 'Y',
            'DEBLEND_NTHRESH' : '10',
            'DEBLEND_MINCONT' : '0.05',
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

    elif instrument == 'CS82':
       _dic = {

            'SEEING_FWHM' : '0.22',
            'THRESH_TYPE' : 'RELATIVE',
            'DETECT_THRESH' : '2.5',
            'FILTER' : 'N',
            'DEBLEND_MINCONT' : '0.5',
            'DEBLEND_NTHRESH' : '32',
            'DETECT_MINAREA'  : '10',
            'ANALYSIS_THRESH' : '2.5',
            'PIXEL_SCALE': '0.187',
            'MAG_ZEROPOINT' : '21.59'
            };

    elif instrument == 'sims':
       _dic = {
            'DETECT_TYPE' : 'CCD',
            'DETECT_MINAREA' : '4',
            'THRESH_TYPE' : 'ABSOLUTE',
            'DETECT_THRESH' :  '0.000004',
            'ANALYSIS_THRESH' : '1.0',
            'FILTER' : 'Y',
            'FILTER_NAME' : 'default.conv',
            'DEBLEND_NTHRESH' : '32',
            'DEBLEND_MINCONT' : '0.15',
            'CLEAN' : 'Y',
            'CLEAN_PARAM' : '1.0',
            'MASK_TYPE' : 'CORRECT',
            'PHOT_APERTURES' : '5',
            'PHOT_AUTOPARAMS' : '2.5,3.5',
            'SATUR_LEVEL' : '50000.0',
            'MAG_GAMMA' : '4.0',
            'GAIN' : '0.0',
            'PIXEL_SCALE' : '1.0',
            'SEEING_FWHM' : '0.3',
            'BACK_TYPE' : 'MANUAL',
            'BACK_VALUE' : '0.0,0.0',
            'STARNNW_NAME' : 'default.nnw',
            'BACK_SIZE' : '64',
            'BACK_FILTERSIZE' : '3',
            'BACKPHOTO_TYPE' : 'GLOBAL',
            'CHECKIMAGE_TYPE' : 'NONE',
            'MEMORY_OBJSTACK' : '3000',
            'MEMORY_PIXSTACK' : '1000000',
            'MEMORY_BUFSIZE' : '1024',
            'VERBOSE_TYPE' : 'NORMAL'
            };

    else:
       _dic = {};


    return (_dic);

# ---

#=======================================================================
def run(filename, params=[], args={}, preset='', temp_dir='', quiet=False):
    """
    Run Sextractor over given image with optional parameters

    run( 'image.fits' )  ->  <bool>

    Default Sextractor(SE) run is done on given FITS image, 'filename',

    $ sex -d > default.sex  &&  sex 'filename' -c default.sex

    This function is just an interface to run SE (v2.8.6) with default
    convolution matrix (default.conv) if no extra arguments are given,

    $ cat default.conv
    CONV NORM
    # 3x3 ``all-ground'' convolution mask with FWHM = 2 pixels.
    1 2 1
    2 4 2
    1 2 1

    (That is the only different when running "sex" from here instead of
    the system's shell: the creation of 'default.conv' file during run)

    It is optional to give SE's command-line (config) arguments, as well 
    as the parameters to compute (and output) from image. See 'args' and 
    'params' (respectively) below. 
    'preset' are some pre-set SE configurations to particular images 
    from instrumentation commonly used in our group.

    'args' is a dictionary datatype interface to SE's command-line arguments.
    The dictionary keys should be the optional argument labels, with the 
    corresponding values.
    E.g.,
    args = {"CATALOG_TYPE" : "FITS_1.0", "PHOT_AUTOPARAMS" : "2.4,3.6"}

    'params' is a list of strings with the parameters (computed from image)
    to output (SE's default.param).
    E.g.,
    params = ['NUMBER', 'X_IMAGE', 'Y_IMAGE']

    If 'params' is not empty, a file called "run_sex.param" will be 
    created in 'temp_dir' directory with params contents and used by 
    "sex" during run. If 'args' *and* 'preset' are used, args values 
    overwrite preset ones in common (and both, SE defaults).
    File (default) with SE convolution matrix (default.conv) will also
    be created inside 'temp_dir'.

    Input:
     - filename : str
        FITS image filename
     - params : [str,]
        Sextractor parameters to output (see SE's default.param)
     - args : {str:str,}
        Sextractor command-line arguments
     - preset : str
        See sextractor.se_presets() documentation for options

    Output:
     - True|False

    """
    
    fits_image = filename;
    if temp_dir != '':
        temp_dir = str(temp_dir)+'/';
    
    # Object features to output..
    #
    if ( params != [] ):
        param_file = temp_dir+'run_sex.param';
        args.update( {'parameters_name' : param_file} );
        fparam = open(param_file,'w');
        for _param in params:
            fparam.write( "%s\n" % (_param) );
        fparam.close();


    # ==========================================================================
    # WRITE DOWN "default.conv" for sextractor if necessary..
    #
    if (not args.has_key('filter_name')  or  args['filter_name']=='default.conv'):
        conv_file = temp_dir+'default.conv';
        args.update( {'filter_name' : conv_file} );
        fp = open(conv_file,'w');
        fp.write("CONV NORM\n");
        fp.write("# 3x3 ``all-ground'' convolution mask with FWHM = 2 pixels.\n");
        fp.write("1 2 1\n");
        fp.write("2 4 2\n");
        fp.write("1 2 1\n");
        fp.close();
    # ==========================================================================


    cargs = se_presets(preset);
    cargs.update( args );

    # Build up sextractor command line..
    #
    cmd_line = '';
    for key in cargs:
        cmd_line = cmd_line + ' -' + string.upper(key) + ' ' + re.sub( "\s+","",str(cargs[key]));

    # Run sex..
    #
    _dev = '';
    if quiet:
        _dev = '&>/dev/null';

    status = os.system( 'sex %s %s %s' % (fits_image,cmd_line,_dev));
    if ( status != 0 ):
        logging.error("Error: Sextractor raised and error code '%s' during the run.",status);
        return (False);

    return (True);

# -

#---------------------------------------------------------
def run_segobj(filename, params=[], args={}, preset='', temp_dir='', quiet=False):
    """ Run Sextractor over given image with specific outputs

    run_segobj( 'image.fits' )  ->  <dict>

    See sextractor.run() help page for description about
    'filename', 'params', 'args' and 'preset' arguments.

    This functions adds to 'params' list the value 'NUMBER' if not yet 
    given. And set some 'args' options to specific values with the goal 
    to have The SEGMENTATION and OBJECTS output versions of the input 
    image, as well as a FITS catalog with 'params' computed features.

    SE's hardcoded argument values:
    CHECKIMAGE_TYPE : "OBJECTS,SEGMENTATION"
    CHECKIMAGE_NAME : objects_name , segmentation_name
    CATALOG_TYPE : "FITS_1.0"
    CATALOG_NAME : catalog_name

    "objects_name", "segmentation_name" and "catalog_name" 
    have values composed by the image 'filename' (without the ".fits" 
    extension) with the suffixes "_obj.fits", "_seg.fits" and "_cat.fit", 
    respectively. These file names are given as output in a dictionary
    with keys ""OBJECTS", "SEGMENTATION" and "CATALOG".


    Input:
     - filename : str
        FITS image filename
     - params : [str,]
        Sextractor parameters to output (see SE's default.param)
     - args : {str:str,}
        Sextractor command-line arguments
     - preset : str
        Options are: 'HST', 'CFHT', 'DC4', 'DC5'

    Output:
     - {'OBJECTS':str, 'SEGMENTATION':str, CATALOG:str}
        SE's OBJECTS, SEGMENTATION, CATALOG output filenames
      
    """

    fits_image = filename;
    
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
    cargs={'checkimage_type' : 'OBJECTS,SEGMENTATION', 'checkimage_name' : objimgname+','+segimgname, 'catalog_type' : 'FITS_1.0', 'catalog_name' : catfilename}
    cargs.update(args);

    # Run sextractor module; output catalogues to "catalog" files..
    #
    out = run(fits_image, params, cargs, preset, temp_dir, quiet);

    if (out == False):
        return (False);

    return ({'OBJECTS':objimgname, 'SEGMENTATION':segimgname, 'CATALOG':catfilename});





def flags_distribuition(hdu_cat):
   
    """ 
    Give number of objects per flag in a fits catalog

    
    input:
    - hdu : <hdulist> it is the hdu fits sextractor catalog, eg. hdu=pyfits.open("catalog_path")

    output:
    - <list> It is a list with the number of object per flag. The '0' index is the flag '0', the index '1' is the flag '1' , the index '2' is the flag '2' , the index '3' is the flag '4', the index '4' the flag '8', the index '5' the flag '16', the index '6' is the flag '32', the index '7' the flag '64', the index '8' the flag '128'.   

    """
    flags_dist=[0,0,0,0,0,0,0,0,0]
    catalogs_data =hdu_cat[1].data
    for i in range(0,len(catalogs_data)):
        flag=flags_counter(hdu_cat,i+1)
        for j in range(0,len(flags_dist)):
            if int(flag[j])==1:
                    flags_dist[j]=flags_dist[j]+1
    return flags_dist     

def flags_counter(hdu_cat,_id):
    """ 
    Returns flags  on a given object

    
    input:
    - hdu : <hdulist> it is the hdu fits sextractor catalog, eg. hdu=pyfits.open("catalog_path")
    - _id : <int> it the object id in the sextractor fits cataloh  

    output:
    - <list> It is a list with information of which flags are 'on' in a given object. If the flag is 'on' the list item has the value '1' and if the flag is 'off' list item value is '0'. The '0' index is the flag '0', the index '1' is the flag '1' , the index '2' is the flag '2' , the index '3' is the flag '4', the index '4' the flag '8', the index '5' the flag '16', the index '6' is the flag '32', the index '7' the flag '64', the index '8' the flag '128'.   

    """

    catalogs_data =hdu_cat[1].data
    i=_id-1
    flag=catalogs_data[i]['FLAGS']
    flag=bin(flag)
    flag=str(flag)
    flag=string.replace(flag,"0b","")
    flags_list=[0,0,0,0,0,0,0,0,0]
    for j in range(0,len(flag)):
            if int(flag)==0:
                flags_list[0]= 1
            else:
                if int(flag[j])==1:
                    flags_list[j+1]=+1
    return flags_list


def read_config(SE_configfile):
    """ Read SE config file (plain: 'key value' ascii) to a dictionary."""
    rs = re.sub;
    
    def valid_line(w):
        return ( not (bool(re.search("^$",w))  or  bool(re.search("^#",w))) );
    lines = filter( valid_line, [ rs("#.*","",rs("\s+"," ",rs("^\s*","",w.rstrip("\n")))) for w in open(SE_configfile).readlines() ]);
    
    config={};
    for each_line in lines:
        _line = string.split(each_line,sep=" ");
        config[_line[0]] = string.join(_line[1:],sep="");
    
    return config;

# \cond
#==========================
if __name__ == "__main__" :

    from sltools.io import config_parser as CP;

    from optparse import OptionParser;

    usage="\n %prog [options] <image.fits>"
    parser = OptionParser(usage=usage);

    parser.add_option('-p',
                    dest='params', default=None,
                    help="Sextractor params file (e,g. default.param)");
    parser.add_option('-c',
                    dest='config', default=None,
                    help="Sextractor config file (e.g, default.sex)");
    parser.add_option('--segment', action='store_true',
                    dest='run_seg', default=False,
                    help="Run sextractor in SEGMENTATION mode");
    parser.add_option('--preset',
                    dest='instrument', default='',
                    help="Use preset SE's configuration for different instruments: DC4, DC5, HST, CFHT");
    parser.add_option('--temp_dir',
                    dest='temp_dir', default='',
                    help="Temporary directory for run-time files.");
    parser.add_option('-q','--quiet', action='store_true',
                    dest='quiet', default=False,
                    help="Make Sextractor run quiet.");
        
    (opts,args) = parser.parse_args();

    run_seg = opts.run_seg;
    params_file = opts.params;
    config_file = opts.config;
    preset = opts.instrument;
    quiet = opts.quiet;
    temp_dir = opts.temp_dir;

    if ( len(args) != 1 ):
        parser.print_usage();
        print "Try '--help' for more info.";
        sys.exit(0)
        
    infits = args[0];

    if (config_file == None):
        config = {};
    else:
        config = read_config(config_file);

    # SE's param file:
    #
    if (params_file):
        params = [ line.rstrip("\n")   for line in open(params_file).readlines() ];
    else:
        params = [];

    # Run SE:
    #
    if ( run_seg ):
        _out = run_segobj(infits, params=params, args=config, preset=preset, temp_dir=temp_dir, quiet=quiet);
        if not (_out):
            sys.exit(1);
        else:
            print "SExtraction output files:";
            print "CATALOG: %s" % (_out['CATALOG']);
            print "OBJECTS: %s" % (_out['OBJECTS']);
            print "SEGMENT: %s" % (_out['SEGMENTATION']);
    else:
        _out = run(infits, params=params, args=config, preset=preset, temp_dir=temp_dir, quiet=quiet);
        if not (_out):
            sys.exit(1);

    sys.exit(0);

# ------------------
# \endcond
