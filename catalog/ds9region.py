""" Module to deal with ds9 region files"""

##@package ds9region

import re;
import string;

def read_cat(regionfile):
    """ Function to read ds9 region file

    Only regions marked with a 'circle' or 'box' are read.
    Color used for region marks (circle/box) are also read
    with the purpose of being used as a 'tag'.
    'x','y','dx','dy' are given in a dictionary as output.
    
    Input:
     - regionfile   :   ASCII (ds9 format) file
     
    Output:
     -> {'image':str,'tag':[],'x':[],'y:[]','dx':[],'dy':[]}
    
    """
    
    out = {'image':'', 'tag':[], 'x':[], 'y':[], 'dx':[], 'dy':[]};

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
                out['tag'].append(color);
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
                out['tag'].append(color);
                continue;
            except AttributeError:
                pass;

    fp.close();
    return out;

if __name__ == '__main__' :
    import sys;
    ret = read_cat(sys.argv[1]);
    print zip(ret['x'],ret['y'],ret['dx'],ret['dy']);
