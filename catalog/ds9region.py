""" Module to deal with ds9 region files"""

##@package ds9region

import re;
import string;


def write_cat(centroids,size=20,marker='circle',color='red',imagefile='null.fits',outputfile='ds9.reg'):
    """ Function to write a ds9 region file given a set of centroids
    
    It works only with a circular 'marker' with fixed
    radius for all (x,y) - 'centroids' - given.
    
    Input:
     - centroids : [(x0,y0),]
     - size : int
     - marker : str
     - outputfile : str

    Output:
     <bool>

    """

    output = open(outputfile,'w');
    
    # DS9 region file header
    output.write("# Region file format: DS9 version 4.1 \n");
    output.write("# Filename: %s \n" % (imagefile));
    output.write("global color="+color+" dashlist=8 3 width=1 font=\"helvetica 10 normal\" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n");
    output.write("image\n");

    # Circle entities
    for o_o in centroids:
        output.write("%s(%d,%d,%d) \n" % (marker,o_o[0],o_o[1],size));
        
    output.close();
    
    print "New file: ds9.reg"
    return (True);
    
# ---

def read_cat(regionfile):
    """ Function to read ds9 region file

    Only regions marked with a 'circle' or 'box' are read.
    'color' used for region marks (circle/box) are given as
    output together with 'x','y','dx','dy' as list in a 
    dictionary. The key 'image' in the output (<dict>) gives
    the filename in the 'regionfile'.
    
    Input:
     - regionfile   :   ASCII (ds9 format) file
     
    Output:
     -> {'image':str,'color':[],'x':[],'y:[]','dx':[],'dy':[]}
    
    """
    
    out = {'image':'', 'color':[], 'x':[], 'y':[], 'dx':[], 'dy':[]};

    fp = open(regionfile,'r');

    for line in fp.readlines():

        if (re.search("^#",line)):
            if (re.search("Filename",line)):
                imagename = string.split(line,"/")[-1];
                out['image'] = imagename.rstrip('\n');
            continue;

        else:
            try:
                _cl = re.search('(?<=color\=).*',line).group();
                color = string.split(_cl)[0];
            except AttributeError:
                pass;

            try:
                _fg = re.sub("\)","",re.search('(?<=box\().*\)',line).group());
                x,y,dx,dy = string.split(_fg,sep=",")[:4];
                out['x'].append(x);
                out['y'].append(y);
                out['dx'].append(dx);
                out['dy'].append(dy);
                out['color'].append(color);
                continue;
            except AttributeError:
                pass;

            try:
                _fg = re.sub("\)","",re.search('(?<=circle\().*\)',line).group());
                x,y,R = string.split(_fg,sep=",")[:3];
                out['x'].append(x);
                out['y'].append(y);
                out['dx'].append(R);
                out['dy'].append(R);
                out['color'].append(color);
                continue;
            except AttributeError:
                pass;

    fp.close();
    return out;

if __name__ == '__main__' :
    import sys;
    import string;
    import numpy;
    import optparse;
    
    parser = optparse.OptionParser();
    
    parser.add_option('-r',
                    dest='regionfile',default=None,
                    help='DS9 Region file');
    parser.add_option('-c',
                    dest='block_data', default=None,
                    help='ASCII file with values in columns whitespaced (block data). Header should be commented (#)');
    parser.add_option('--cols_xy',
                    dest='xy_columns',default='1,2',
                    help='Column id. So far, just where "x" and "y" are concerns so far');
    parser.add_option('--color',
                    dest='mark_color',default='green',
                    help='Color to use for region file markers');
    
    (opts,args) = parser.parse_args();
    
    reg = opts.regionfile;
    cat = opts.block_data;
    cols = opts.xy_columns;
    color = opts.mark_color;
    
    if not (reg or cat) :
        parser.print_help();
        parser.exit(msg='No (ds9) region file nor block data given.\n');
        
    xcol,ycol = string.split(cols,sep=',');
    xcol = int(xcol);
    ycol = int(ycol);

    if cat:
        data = numpy.loadtxt(cat,comments='#').T;
        x,y = data[xcol-1],data[ycol-1];
        x = x.tolist();
        y = y.tolist();
        write_cat(zip(x,y),imagefile='F814W_sci.fits',color=color);
        
    if reg:
        ret = read_cat(reg);
        print zip(ret['x'],ret['y'],ret['dx'],ret['dy']);
