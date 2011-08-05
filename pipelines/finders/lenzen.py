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
        self.sections = ['INPUT','SCALE_INPUT','SMOOTHING','DETECTION','OUTPUT']
        self.head_lines = []
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
                                'Mindiff':'0.3',
                                'Minecc':'0.7',
                                'Maxthick':'30',
                                'Use_center_of_gc':'no',
                                'Center_x':'',
                                'Center_y':''
                            }
        self['OUTPUT'] = {'Extra_output':'none',
                            'Outputfile':'output_img.fits',
                            'Data_directory':'data/',
                            'Catalogfile':'output_cat.fit',
                            'Catalogtype':'large'
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
    
    def __head_lines__(self,cfg):
        head_lines = [];
        while True:
            pos = cfg.tell();
            line = cfg.readline();
            if (line.find('[') >= 0) and (line.find('#') != 0): break;
            if (line.rstrip('\n') != ''): head_lines.append(line);
        self.head_lines = head_lines
        self.begin_ini = pos
        
    def write_config(self,filename='input.par'):
        """ Function to write an INI file: [section] key:value """
        fp = open(filename,'w')
        [ fp.write(line+'\n') for line in self.head_lines ]
        fp.close()
        for section in self.sections:
            fp = open(filename,'a')
            _sec = self[section]
            config.add_section(section);
            for key,value in _sec.items():
                config.set(section, key, value);
        config.write(fp);
        fp.seek(0);
        

    
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
    Dcatin = ascii_data.dict_from_csv(catfile[:-4]+'.txt',columns,header_lines=10,delimiter=' ');

    # Lenzen outputs the 'Y' positions inverted. Lets fix that:
    y_size = pyfits.getheader('data/'+imgfile)['naxis1'];
    Y_inv = Dcatin['Y'];

    Dcatin['Y'] = [ str(y_size-int(i)) for i in Y_inv ];
    Catin = fits_data.dict_to_tbHDU(Dcatin);
    print Dcatin.keys();
    print Catin.data
    # Select entries to output arc data
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
    pyfits.writeto(finalimg,diff[::-1],clobber=True)
    # ----------
    
    if clean:
        os.system('rm -rf temp/ %s %s' % (imgfile,catfile));

    return finalcat;

# --

def run(imgfile=None, inputfile=None, clean=False):

    if inputfile == None:
    
        if imgfile == None:
            logging.error("Neither an image file nor an input.par file given. Nothing to do.");
            return False;
        else:
            inputfile = init(imgfile,ecc,outputimage,outputcat);

    else:
        pass;
        
    outputimage = 'lenzen_out.fits';
    outputcat = 'lenzen_out.cat';
    ecc = 0.7;
    
    # Initialize Horesh: read out necessary values from header
    # and write it down to the input file.
    if inputfile is None:
        pass;
        
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
