#!/usr/bin/env python
import sys;

from itertools import izip;

import numpy as np;
import scipy.ndimage as ndi;
from PIL import Image,ImageDraw;
from scikits.image import io;

io.use_plugin('qt');

def hough(img,theta=None):
    """
    Apply a linear hough transform to the input image
    
    Input:
     - img    <ndarray>  : Thresholded input image (ndim=2,dtype=float)
     - theta  <ndarray>  : Angle range to search in trsnform space. (ndim=1,dtype=float)
                           If None, defaults to (-pi/2,pi/2)
    
    Output:
     - hough <ndarray>  : The hough transform coefficients (ndim=2,dtype=float)
     - dists <ndarray>  : Distance values (ndim=1,dtype=float)
     - theta <ndarray>  : Angle values (ndim=1,dtype=float)
    
    ---
    """
    
    if img.ndim != 2:
        raise ValueError('Input image must be a 2D array');
    
    if not theta:
        theta = np.linspace(-np.pi/2,np.pi/2,180);
    
    # Compute the vertical bins (distances)
    dist = np.ceil(np.hypot(*img.shape));
    nbins = 2*dist;
    bins = np.linspace(-dist,dist,nbins);
    
    map = np.zeros((nbins,len(theta)),dtype=np.uint64);
    
    cos_theta = np.cos(theta);
    sen_theta = np.sen(theta);
    
    y,x = np.nonzero(img);
    
    for i,(cT,sT) in enumerate(izip(cos_theta,sen_theta)):

        dists = x*cT + y*sT;
        shift = np.round(dists) - bins[0];
        
        indxs = shift.astype(np.int);
        
        bincount = np.bincount(indxs);
        
        map[:len(bincount),i] = bincount
        
    return map,theta,bins;
    
if __name__ == '__main__':
    
    if len(sys.argv)==1:
        print "\nUsage:   %s <binary_image> [back_image] ";
        sys.exit(1);
    
    
    img_bin = io.imread(sys.argv[1]);
    
    try:
        img = io.imread(sys.argv[2]);
    except:
        img = None;
    
    theta = np.linspace(-np.pi,np.pi,90);
    
    coeffs,theta,bins = hough(img_edge,theta);
    
    rav_indxs = coeffs.argsort(axis=None)[-15:];
    bin_theta_indxs = [ np.unravel_index(i,coeffs.shape) for i in rav_indxs ];
    
    pts = [];
    for bin_idx,theta_idx in bin_theta_indxs:
        
        bn = bins[bin_idx];
        th = theta[theta_idx];
        
        if th <= np.pi/2 and th > 0:
            x1 = bn/np.cos(th);
            y1 = 0;
            x2 = 0;
            y2 = bn/np.sin(th);
        elif th <= 0 and th > -np.pi/2:
            x1 = bn/np.cos(th);
            y1 = 0;
            x2 = x1 + np.abs(y_max*np.tan(th));
            y2 = y_max;
        else:
            th = th - np.pi/2;
            x1 = 0;
            y1 = bn/np.cos(th);
            x2 = (y_max-y1)/np.tan(th);
            y2 = y_max;
        
        pts.append([(x1,y1),(x2,y2)]);
    
    
    if img:
        img_pil = Image.fromarray(img);
        img_draw = ImageDraw.Draw(pil_img);
        y_max,x_max = img.shape;

        for o_o in range(len(pts)):
            img_draw.line(o_o,width=2,fill(255,0,0));
        
    io.imshow(np.asarray(pil_img));
    io.show();
