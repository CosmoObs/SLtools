from scipy import weave;
import numpy as np;

def tanh(img,slope=100,cut=None):
    """
    Modifies image intensities scale
    
    Input:
     - img <ndarray>
     - slope <float> : tanh('slope'*image)
     - cut <float> : If 'cut' is not given, histogram's maximum is used
    
    Output:
     - ndarray
    
    ---
    """
    
    new_img = np.zeros(img.shape,dtype=img.dtype);
    
    min = img.min();
    max = img.max();
    
    if cut == None:
        steps = 10000;
        delta_h = float(max-min)/(steps-1);
        hist,bins = np.histogram(img.flatten(),bins=steps);
        cut = min + (hist.argmax()+1) * delta_h;
        print "Hist max:",cut
    
    new_img = np.tanh(slope*(img-cut)) + 1.;
    
    return new_img;

# --
def histogram(image, percent=99):
    """
    Modifiies image's histogram

    Input:
     - image <ndarray>     : Image array (ndim=2,dtype=float)
     - percent <float>     : Range (percent) of image's intensity to maintain

    Output:
     - new_image <ndarray> : Modified image

    ---
    """
    
    steps = 10000;
    
    width,height = image.shape;
    lowercut, uppercut = 0,0;
    
    size = image.size;
    min = image.min();
    max = image.max();
    
    if min==max:
        return False;
    
    delta_h = float(max-min)/(steps-1);
    hist,bins = np.histogram(image.ravel(),bins=steps);
    
    #histmax_val = hist.max();
    histmax_ind = hist.argmax();
    
    lowercut = min + (histmax_ind+1) * delta_h;

    hist_low = hist[:histmax_ind+1];
    hist_high = hist[histmax_ind+1:];
    hist_cum = np.cumsum(hist_high);
    hist_normcum = hist_cum/float(image.size-hist_low.sum());

    hist_cutlen = np.where(hist_normcum < percent/100.)[0].size;
    hist_cutlen += hist_low.size;
    uppercut = min + (hist_cutlen) * delta_h;
    
    image_new = image.copy();
    image_new[np.where(image<lowercut)] = lowercut;
    image_new[np.where(image>uppercut)] = uppercut;
    
    return image_new;
