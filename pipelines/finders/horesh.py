#!/usr/bin/env python
import sys
import os

import pyfits;
import pywcs;
import numpy;

import logging;
import re;
import math as m;

from sltools.image import header_funcs;
from sltools.pipelines import finders;
from sltools.catalog import fits_data;

finderspath = os.path.abspath(finders.__path__[0]);

def calc_magzpt(hdr):
    """Compute Magnitude zero-point from header keywords"""
    
    # HST images
    phzpt = hdr.get('PHOTZPT',default=0);
    pflam = hdr.get('PHOTFLAM',default=0);
    etime = hdr.get('EXPTIME',default=0);
    logging.debug("Header PHOTZPT,PHOTFLAM,EXPTIME: %s,%s,%s",phzpt,pflam,etime);
    if not etime:
        logging.error("No exposure time were found.");
        return False;
    
    zeropt = -2.5 * m.log(pflam,10) + phzpt;
    
    ncomb = hdr.get('NCOMBINE',default=1);
    ttime = hdr.get('TEXPTIME',default=0);
    atime = ttime/ncomb;
    logging.info("Image Zero-point,Avg-time: %s,%s",zeropt,atime);    

    # Mag zero-point for N combined (average or median) images
    magzpt = zeropt + 2.5 * m.log(atime,10);
    
    return magzpt;
    
    
def init(imagefile):

    # ------------------------------
    # List of images for processing
    # ------------------------------
    hdr = pyfits.getheader(imagefile);
    y_size = hdr['naxis1'];
    x_size = hdr['naxis2'];
    exptime = hdr.get('texptime',hdr.get('exptime',0));
    if not exptime:
        logging.error("Exposure time was not found on given image header.");
        return False;

    imagelistfile = 'images.lst';
    fp = open(imagelistfile,'w')
    fp.write("%s %d %d %d\n" % (imagefile,x_size,y_size,exptime))
    fp.close()
    # ------------------------------


    # ----------------------
    # Input parameters file
    # ----------------------
    inputfile = 'input.in';
    fp = open(inputfile,'w');

    # List-of-images filename
    fp.write(str(imagelistfile)+"\n");
    
    # Output filename, with the number of detected files
    fp.write("nArcs.out\n");

    # Pixel Scale
    dimpix = header_funcs.read_pixelscale(hdr);
    fp.write(str(dimpix)+"\n");

    # Lenght-to-width threshold value
    LW = 5;
    fp.write(str(LW)+"\n");

    # Magnitude Zeropoint
    magzpt = hdr.get('MAGZP',default=None)
    if not magzpt:
        magzpt = calc_magzpt(hdr);
    fp.write(str(magzpt)+"\n");

    # Phot flam param (HST)
    photflam = hdr.get('PHOTFLAM',0);
    fp.write(str(photflam)+"\n");

    # GAMBIARRA!!!   Path to SE files
    os.system('cp -R %s/horesh_sex_files .' % (finderspath))
    fp.write("./horesh_sex_files\n");

    fp.close();
    
    return inputfile;
    
# --
def process_output(rootname,clean=False):

    # Clean final but unwanted files:
    os.system('rm -f %s %s' % (rootname+'.cat_old',rootname+'new.cat'));
    os.system('rm -f %s %s' % (rootname+'new_final2.fits',rootname+'.cat'));
    # Clean tmp files:
    if clean:
        os.system('rm -rf images.lst input.in horesh_sex_files');
    # Rename final image:
    os.system('mv %s %s' % (rootname+'new_final.fits',rootname+'.final.fits'));    
    finalimg = rootname+'.final.fits';
    finalcat = rootname+'.final.cat';
    
    return finalcat;

# --
def run(imgfile='img.fits', inputfile=None, clean=False):

    # Initialize Horesh: read out necessary values from header
    # and write it down to the input file.
    if not inputfile:
        inputfile = init(imgfile);

    os.system('horesh < %s' % (inputfile));

    # Process outputs
    rootname = imgfile[:-5];
    finalcat = process_output(rootname,clean);
    
    return finalcat;
    
# ---
if __name__ == '__main__':

    import optparse;

    usage="\n  %prog [options] <image.fits>"
    parser = optparse.OptionParser(usage=usage);

    parser.add_option('-p',
                        dest='input', default=None,
                        help="Horesh' arcfinder input parameters file");
                        
    (opts,args) = parser.parse_args();

    if len(args) == 0:
        parser.error("No input images were given.");

    inputfile = opts.input;

    imgfile = args[0];
    
    # Run Horesh's arcfinder
    finalcat = run(imgfile,inputfile)

    FITSout = True;
    
    # To finish the process, convert the ASCII table into FITS if asked for.
    if FITSout:
#        finalfit = finalcat[:-4]+".fit";
        finalfit = 'Catalog_arcs_horesh.out.fit';
        cat = numpy.loadtxt(finalcat).T;
        if cat.ndim == 1:
            cat.shape = (len(cat),1);
        dcat = {'ID':cat[0], 'X':cat[1], 'Y':cat[2], 'STmag':cat[3], 'Length':cat[4], 'Width':cat[5], 'L/W':cat[6], 'Area':cat[7]};
        tbhdu = fits_data.dict_to_tbHDU(dcat,tbname="Horesh_arcsfound");
        tbhdu.writeto(finalfit,clobber=True);
#        os.system('rm -f %s' % (finalcat));

    sys.exit(0);
