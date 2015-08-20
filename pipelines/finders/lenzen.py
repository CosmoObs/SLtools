#!/usr/bin/env python
import sys;
import os;
import logging

"""Python interface fo Lenzen's arc detector"""

import pyfits
import ConfigParser
import re

from sltools.catalog import ascii_data
from sltools.catalog import fits_data
from sltools.io import config_parser as cp

# ---

class Config(dict):
    """Store input.par section contents"""
    
    def __init__(self,imagefile=None,configfile=None):
        self.default()
        self.sections = [('INPUT',['Infile','Resize','x1','y1','x2','y2','Origin_x','Origin_y']),
                        ('SCALE_INPUT',['Mode','Min','Max']),
                        ('SMOOTHING',['Steps','Tau','K','Sigma','Rho','Outputfile','Outputformat']),
                        ('DETECTION',['Minlevel','Mindiff','Minecc','Maxthick','Use_center_of_gc','Center_x','Center_y']),
                        ('OUTPUT',['Extra_output','Outputfile','Data_directory','Catalogfile','Catalogtype'])]
        self.head_lines = ["Temp_dir:   temp\n","Data_dir:   data\n"]
        self.begin_ini = 0
        if configfile:
            self.read_config(configfile)
        if imagefile:
            input_sec = self['INPUT']
            input_sec['Infile'] = imagefile
        
    
    def default(self):
        """Default Arc Detector config params"""
        self['INPUT'] = {'Infile':'<?>',
                            'Resize':'no',
                            'x1':'',
                            'y1':'',
                            'x2':'',
                            'y2':'',
                            'Origin_x':'0',
                            'Origin_y':'0'
                            }
        self['SCALE_INPUT'] = {'Mode':'4',
                                'Min':'auto',
                                'Max':'99%'
                                }
        self['SMOOTHING'] = {'Steps':'1',
                                'Tau':'10',
                                'K':'0.00001',
                                'Sigma':'1',
                                'Rho':'10',
                                'Outputfile':'smooth',
                                'Outputformat':'fits'
                            }
        self['DETECTION'] = {'Minlevel':'0.1',
                                'Mindiff':'0.1',
                                'Minecc':'0.6',
                                'Maxthick':'50',
                                'Use_center_of_gc':'no',
                                'Center_x':'',
                                'Center_y':''
                            }
        self['OUTPUT'] = {'Extra_output':'none',
                            'Outputfile':'output_img.fits',
                            'Data_directory':'data/',
                            'Catalogfile':'output_cat.txt',
                            'Catalogtype':'small'
                            }

    def get(self,section,key):
        """ = """
        sec = self[section]
        return sec[key]
        
    def set(self,section,key,value):
        """ = """
        sec = self[section]
        if not sec.has_key(key):
            logging.error("Key %s not found in section %s" % (key,section))
            return
        sec[key] = str(value)
        
    def __head_lines__(self,cfg):
        head_lines = [];
        while True:
            pos = cfg.tell();
            line = cfg.readline();
            if (line.find('[') >= 0) and (line.find('#') != 0): break;
            if (line.rstrip('\n') != ''): head_lines.append(line);
        self.head_lines = head_lines
        self.begin_ini = pos
        
    def read_config(self,configfile):
        """ Returns a dictionary containing config sec->key:val structure """
        parser = ConfigParser.ConfigParser()
        parser.optionxform = str
        cfg = open(configfile);
        self.__head_lines__(cfg)
        cfg.seek(self.begin_ini)
        parser.readfp(cfg)
        config = {}
        for section in parser.sections():
            sec = self[section.upper()]
            for k,v in parser.items(section):
                sec[k.capitalize()] = re.sub('\s+','',re.sub('#.*','',v))
        cfg.close()
    
    def write_config(self,filename='input.par'):
        """ Function to write an INI file: [section] key:value """
        
        fp = open(filename,'w')
        [ fp.write(line) for line in self.head_lines ]
        fp.close()
        
        fp = open(filename,'a')
        for section,options in self.sections:
            fp.write("\n")
            _sec = self[section]
            fp.write("["+section+"]\n")
            for key in options:
                fp.write("%s: %s\n" % (key.capitalize(),_sec[key]))
        fp.close()
        
"""
Temp_dir:      temp
Data_dir:      data

[INPUT]
Infile: <?>                 # FITS or PGM file name
Resize: no                  # 'yes' or 'no' 
x1:                         # Rectangle to use
y1:
x2:
y2:
Origin_x: 0
Origin_y: 0

[SCALE_INPUT]
Mode:  5                    # Scaling mode,recommended: 4 or 5   
Min:   <?>                  # 'auto' or value
Max:   <?>
# Scaling mode:
# 0 - no scaling
# 1 - I(x)/Max, Intensities above Max are set to 1
# 2 - intensities below Min are set to Min, above Max are set to Max
# 3 - I(x)/Max ,intensities below Min are set to Min, above Max are set to 1
# 4 - (I(x)-Min)/(Max-Min),intensities below Min are set to 0, above Max are set to 1
# 5 -  sqrt((I(x)-Min)/(Max-Min)),intensities below Min are set to 0, above Max are set to 1

[SMOOTHING]
Steps:         1            # Numsteps*Tau define duration of diffusion process
#                           # After each step the result can be stored in an outputfile (see below)
Tau:           <?>          # typical range[1-15]
K:             0.0001       # typical range[0.001-.00001]
Sigma:         2            # typical range [1-5]
Rho:           <?>          # typical value for arclets ~7-10, arcs ~10-20
Outputfile:    none         # 'none' or filename without ending
Outputformat:  pgm          # 'pgm' or 'fits'

[DETECTION]
Minlevel:      <?>          # [0,1] local maxima below are ignored
Mindiff:                    # [0,1] minimal aquired difference between object and background
#                           # default: same as Minlevel
Minecc:        .7           # [0,1] minimal aquired eccenticity
Maxthick:      10           # maximal thickness 
Use_center_of_gc: no        # 'yes' or 'no'
Center_x:                   # due to original image size
Center_y: 

[OUTPUT]
Extra_output:   none        # 'none' or file name
#                           # ending is automatically attachend (PGM/PPM)
Outputfile:     <?>         # file name with fits- or PPM-ending
Data_directory: DATA/       # where to put the tempfiles, e.g. ./
Catalogfile:    none        # 'none' or file name
Catalogtype:    small       # 'small' or 'large'
"""

# --
def run(imagefile=None, inputfile=None, clean=False):
    """
    Function to run Lenzen arc detector pipeline
    
    Input:
     - imagefile   str : FITS image filename to segment
     - inputfile   str : Lenzen's input filename (input.par)
     - clean
    
    Output:
     - Lenzen's code outputs. Stay calm, "ls" your dir ;P
    
    ---
    """
    
    if (inputfile == None) and (imagefile == None):
        logging.error("Neither an image file nor an input.par file given. Nothing to do.");
        return False;
    
    config = Config(imagefile=imagefile,configfile=inputfile)
    
    imageinput = config.get('INPUT','Infile')
    outputimage = re.sub(".fits","_out_img.fits",imageinput)
    config.set('OUTPUT','Outputfile',outputimage)
    outputcat = re.sub(".fits","_out_cat.txt",imageinput)
    config.set('OUTPUT','Catalogfile',outputcat)
    smoothfile = re.sub(".fits","_out_smooth",imageinput) 
    config.set('SMOOTHING','Outputfile',smoothfile)
    
    lenzeninput = "run_lenzen.par"
    config.write_config(lenzeninput)
    
    os.system('mkdir temp data &> /dev/null')
    os.system('ln -sf ../%s data/. &> /dev/null' % (imageinput))
    os.system('lenzen %s' % (lenzeninput))
    
    process_output(outputimage,outputcat,clean)
    
    return

# --
def process_output(imgfile,catfile,clean=False):

    from scipy import ndimage
    import numpy as np
    
    # ----------
    # Catalog
    catalog_file = catfile
    os.system("sed -e 's/^[ ]*//' %s | tr -s ' ' > %s" % ('data/'+catfile,catalog_file))
    columns = ['X','Y','ecc','lambda1','lambda2','A','B','theta','size']
    Dcatin = ascii_data.dict_from_csv(catalog_file,columns,header_lines=10,delimiter=' ')

    # Lenzen outputs the 'Y' positions inverted. Lets fix that:
    hdr = pyfits.getheader('data/'+imgfile)
    y_size = hdr['naxis1']
    x_size = hdr['naxis2']
    Y_inv = Dcatin['Y']
    
    finalimg = re.sub("_out_img.fits","_out_result.fits",imgfile)
    if len(Y_inv):
        Dcatin['Y'] = [ str(y_size-int(i)) for i in Y_inv ]
        Catin = fits_data.dict_to_tbHDU(Dcatin)

        # Select entries to output arc data
        Catout = fits_data.sample_entries(Catin,ecc=0.7)
        Catout.name = "Lenzen_arcsfound"
        finalcat = catalog_file[:-4]+'.fits'
        Catout.writeto(finalcat,clobber=True)
        
        # Image
        img = ndimage.imread('data/result.ppm')
        comb = np.maximum(img[...,0],img[...,1])
        diff = img[...,2] - comb
        pyfits.writeto(finalimg,diff[::-1],clobber=True)
        
    else:
        blank_array = np.zeros((y_size,x_size))
        pyfits.writeto(finalimg,blank_array,clobber=True)
    
    if clean:
        os.system('rm -rf temp/ %s %s' % (imgfile,catfile));

    return

# ---
if __name__ == '__main__':

    import optparse;

    usage="\n  %prog [options] <image.fits>"
    parser = optparse.OptionParser(usage=usage);

    parser.add_option('-p',
                        dest='input', default=None,
                        help="Lenzen' arcfinder input parameters file (-create)");
                        
    (opts,args) = parser.parse_args();

    if len(args) == 0:
        imgfile = None
    else:
        imgfile = args[0];
    
    inputfile = opts.input;
    
    # Run Lenzen's arcfinder
    run(imgfile,inputfile)

    sys.exit(0);
