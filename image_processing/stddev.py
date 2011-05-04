import sys;

import numpy as np;
import scipy.ndimage as ndi;

def mean(img, size=3):
    """
    Simple mean value filter

    Input:
     - img <ndarray>
     - size <int> : size of the filter window

    Output:
     <ndarray>
    
    ---
    """
    size = int(size);
    kernel = np.ones((size,size)) / float(size**2);
    return ndi.convolve(img,kernel)

# --

def median(img, size=3):
    """
    Simple median value filter
    
    Input:
     - img <ndarray>
     - size <int> : size of the filter window
    
    Output:
     <ndarray>
    
    ---
    """
    
    return ndi.median_filter(img,size);

# --

def gaussian(img, sigma=[3,3]):
    """
    Simple gaussian filter
    
    Input:
     - img <ndarray>
     - sigma <[int,int]> : sigma window
    
    Output:
     - <ndarray>
    
    ---
    """

    return ndi.gaussian_filter(img,sigma);
    
# --

def stddev(img, size=3):
    """
    Standard deviation transform
    
    Input:
     - img <ndarray>
     - size <int> : size of the window to compute StdDev
    
    Output:
     ndarray
    
    ---
    """
    mean_img = mean(img,size);
    diff_img = img.copy();
    diff_img[1:-1,1:-1] = np.power(mean_img[1:-1,1:-1]-img[0:-2,0:-2],2) \
                            +np.power(mean_img[1:-1,1:-1]-img[2:,2:],2) \
                            +np.power(mean_img[1:-1,1:-1]-img[0:-2,2:],2) \
                            +np.power(mean_img[1:-1,1:-1]-img[2:,0:-2],2) \
                            +np.power(mean_img[1:-1,1:-1]-img[1:-1,1:-1],2) \
                            +np.power(mean_img[1:-1,1:-1]-img[1:-1,0:-2],2) \
                            +np.power(mean_img[1:-1,1:-1]-img[1:-1,2:],2) \
                            +np.power(mean_img[1:-1,1:-1]-img[0:-2,1:-1],2) \
                            +np.power(mean_img[1:-1,1:-1]-img[2:,1:-1],2);
    return diff_img;
    
