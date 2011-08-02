import sltools.image.sextractor as se
import numpy as np


def run(file_name, params=[], args={}, preset='', temp_dir='', quiet=False, arg_min={}, arg_max={}, arg_step={}, Filter_Name='', plots=False):

    """ Run Sextractor over given image and a range of inputs

    See sextractor.run() help page for description about
    'filename', 'params', 'args' and 'preset' arguments.

    See sextractor.run_segobj() to hardcoded arguments 
    

    "objects_name", "segmentation_name" and "catalog_name" 
    have values composed by the image 'filename' (without the ".fits" 
    extension) with the suffixes "_obj", "_seg" and "_cat", respectively 
    plus a parameter combination id. The table with parameter combination id 
    and the value of those parameters are in the  ascii file with suffix "_variables.txt".


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

    arg_min_t={'DETECT_THRESH' : '2.5', 'DEBLEND_MINCONT' : '0.5',  'DEBLEND_NTHRESH' : '32', 'DETECT_MINAREA'  : '10', 'ANALYSIS_THRESH' : '2.5'}
    arg_max_t={'DETECT_THRESH' : '2.5', 'DEBLEND_MINCONT' : '0.5',  'DEBLEND_NTHRESH' : '32', 'DETECT_MINAREA'  : '10', 'ANALYSIS_THRESH' : '2.5'}
    arg_step_t={'DETECT_THRESH' : '2.5', 'DEBLEND_MINCONT' : '0.5',  'DEBLEND_NTHRESH' : '32', 'DETECT_MINAREA'  : '10', 'ANALYSIS_THRESH' : '2.5'}
    arg_min_t.update(arg_min)
    arg_max_t.update(arg_max)    
    arg_step_t.update(arg_step)
    fits_image = filename;
    _rootname = string.split( string.replace( fits_image,".fits","" ), sep="/" )[-1]
    f_out=open(_rootname+"_variables.txt","w")
    detect_thresh=np.arange(arg_min_t['DETECT_THRESH'],arg_max_t['DETECT_THRESH'], arg_step_t['DETECT_THRESH'])
    deblend_min=np.arange(arg_min_t['DEBLEND_MINCONT'],arg_max_t['DEBLEND_MINCONT'], arg_step_t['DEBLEND_MINCONT'])
    deblend_nthresh=np.arange(arg_min_t['DEBLEND_NTHRESH'],arg_max_t['DEBLEND_NTHRESH'], arg_step_t['DEBLEND_NTHRESH'])
    detec_min=np.arange(arg_min_t['DETECT_MINAREA'],arg_max_t['DETECT_MINAREA'], arg_step_t['DETECT_MINAREA'])
    analysis_thresh=np.arange(arg_min_t['ANALYSIS_THRESH'],arg_max_t['ANALYSIS_THRESH'], arg_step_t['ANALYSIS_THRESH'])
    Filter_=['Y','N']
    
    if preset!='':
        cargs = se_presets(preset);
        args.uptdate(cargs)    
    if Filter_Name=='':
        args.update({'FILTER_NAME' : temp_dir+'default.conv'})
    for x in range(0,len(Filter_)):
        args.update({'FILTER' : Filter[x]})
        for i in range(0,len(detect_thresh)):
            args.update({'DETECT_THRESH' : detect_thresh[i]})
            for j in range(0,len(deblend_min)):
                args.update({'DEBLEND_MINCONT' : deblend_min[j]})
                for k in range(0,len(deblend_nthresh)):
                    args.update({'DEBLEND_NTHRESH' : deblend_nthresh[k]})
                    for w in range(0,len(detec_min)):
                        args.update({'DETECT_MINAREA' : detec_min[w]})
                        for y in range(0,len(analysis_thresh)):
                            args.update({'ANALYSIS_THRESH' : analysis_thresh[y]})
                            objimgname = _rootname+'var'+str(i)+str(j)+str(k)+str(w)+str(y)+'_obj.fits';
                            segimgname = _rootname+'var'+str(i)+str(j)+str(k)+str(w)+str(y)+'_seg.fits';
                            f_out.write(str(i)+str(j)+str(k)+str(w)+str(y)+'        '+'DETECT_THRESH='+str(args['DETECT_THRESH'])+'    '+'DEBLEND_MINCONT='+str(args['DEBLEND_MINCONT'])+'    '+'DEBLEND_NTHRESH='+str(args['DEBLEND_NTHRESH'])+'    '+'DETECT_MINAREA='+str(args['DETECT_MINAREA'])+'ANALYSIS_THRESH='+str(args['ANALYSIS_THRESH'])+'    '+"FILTER="+Filer[x]+'    \n \n');    
                            args.update({'checkimage_name' : objimgname+','+segimgname})
                            run_segobj(filename, params, args, preset, temp_dir, quiet)            
    f_out.close()       
           
            

    return 0

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
                    help="Use preset SE's configuration for different instruments: e.g. DC4, DC5, HST, CFHT");
    parser.add_option('--temp_dir',
                    dest='temp_dir', default='',
                    help="Temporary directory for run-time files.");
    parser.add_option('-q','--quiet', action='store_true',
                    dest='quiet', default=False,
                    help="Make Sextractor run quiet.");
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
        config = se.read_config(config_file);

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



  
