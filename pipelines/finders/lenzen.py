#!/usr/bin/env python
import sys;
import os;

"""Python interface fo Lenzen's arc-finder"""

import pyfits;
from sltools.catalog import ascii_data;
from sltools.catalog import fits_data;

# ---
def init(imagefile,eccentricity=0.7,outputimage='lenzen_out.fits',outputcat='lenzen_out.cat'):

    # Atualizar variaveis:
    # [INPUT] : 'Infile'
    # [SCALE_INPUT] : 'Min' , 'Max'
    # [OUTPUT] : 'Outputfile' , 'Catalogfile'
    #
    _config_large="""
    Temp_dir:      temp
    Data_dir:      data

    [INPUT]
    Infile: %s                 # FITS or PGM file name

    Origin_x: 0
    Origin_y: 0


    [SCALE_INPUT]
    Mode:  4                    # Scaling mode,recommended: 4 or 5   
    Min:   %s                  # 'auto' or value
    Max:   %s

    [SMOOTHING]
    Steps:         1            # Numsteps*Tau define duration of diffusion process
                                # After each step the result can be stored in an outputfile (see below)
    Tau:           13           # typical range[1-15]
    K:             0.0001       # typical range[0.001-.00001]
    Sigma:         2            # typical range [1-5]
    Rho:           15           # typical value for arclets ~7-10, arcs ~10-20

    Outputfile:    smoothing    # 'none' or filename without ending
    Outputformat:  pgm          # 'pgm' or 'fits'

    [DETECTION]
    Minlevel:      0.3          # [0,1] local maxima below are ignored
    Mindiff:       0.3          # [0,1] minimal aquired difference between object and background
    Minecc:        %s           # [0,1] minimal aquired eccentricity
    Maxthick:      10           # maximal thickness 

    [OUTPUT]
    Extra_output:   none        # 'none' or file name
                                # ending is automatically attachend (PGM/PPM)
    Outputfile:     %s          # file name with fits- or PPM-ending
    Data_directory: DATA/       # where to put the tempfiles, e.g. ./

    Catalogfile:    %s        # 'none' or file name
    Catalogtype:    small       # 'small' or 'large'

    """
    
    os.system('mkdir temp data DATA');
    os.system('ln -s ../%s data/.' % (imagefile));
    param_file = 'input.par';
    ecc = eccentricity;
    outimg = outputimage;
    outcat = outputcat;
    fin = open(param_file,'w');
    print >> fin, _config_large % (imagefile,'auto','99%',ecc,outimg,outcat)

    return param_file;

# --

def process_output(imgfile,catfile,clean=False):

    from scipy import ndimage;
    import numpy as np;
    
    # ----------
    # Catalog
    os.system("sed -e 's/^[ ]*//' %s | tr -s ' ' > %s" % ('data/'+catfile,catfile[:-4]+'.txt'));

    columns = ['X','Y','ecc','lambda1','lambda2','A','B','theta','size'];
    Catin = fits_data.dict_to_tbHDU(ascii_data.dict_from_csv(catfile[:-4]+'.txt',columns,header_lines=10,delimiter=' '));
    
    # Lenzen outputs the 'Y' positions inverted. Lets fix that:
    y_size = pyfits.getheader('data/'+imgfile)['naxis1'];
    
#    Catin.data.field('Y') = y_size - Catin.data.field('Y');
    Catout = fits_data.sample_entries(Catin,ecc=0.7);
    Catout.name = "Lenzen_arcsfound";
    finalcat = 'Catalog_arcs_lenzen.out.fit';
    Catout.writeto(finalcat,clobber=True);
    # ----------
    
    # ----------
    # Image
    img = ndimage.imread('data/result.ppm');
    comb = np.maximum(img[...,0],img[...,1]);
    diff = img[...,2] - comb;
    finalimg = 'Image_arcs_lenzen.out.fits';
    pyfits.writeto(finalimg,diff,clobber=True)
    # ----------
    
    if clean:
        os.system('rm -rf temp/ %s %s' % (imgfile,catfile));

    return finalcat;

# --

def run(imgfile='img.fits', inputfile=None, clean=False):

    outputimage = 'lenzen_out.fits';
    outputcat = 'lenzen_out.cat';
    ecc = 0.7;
    
    # Initialize Horesh: read out necessary values from header
    # and write it down to the input file.
    if not inputfile:
        inputfile = init(imgfile,ecc,outputimage,outputcat);

    os.system('lenzen %s' % (inputfile));
    
    finalcat = process_output(outputimage,outputcat);

    return finalcat;
    
# ---
if __name__ == '__main__':

    import optparse;

    usage="\n  %prog <image.fits>"
    parser = optparse.OptionParser(usage=usage);

    parser.add_option('-p',
                        dest='input', default=None,
                        help="Lenzen' arcfinder input parameters file (-create)");
                        
    (opts,args) = parser.parse_args();

    if len(args) == 0:
        parser.error("No input images were given.");

    inputfile = opts.input;

    imgfile = args[0];
    
    # Run Horesh's arcfinder
    finalcat = run(imgfile,inputfile)

    sys.exit(0);
