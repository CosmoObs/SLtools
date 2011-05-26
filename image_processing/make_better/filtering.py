import random
import numpy as np
import scipy as sc
import scipy.ndimage as ndi

def average(img, size=3):
  """ Simple mean value filter

  Input:
   - img <ndarray>
   - size <int> : size of the filter window

  Output:
   <ndarray>
  """
  
  size = int(size)
  kernel = np.ones((size,size)) / float(size**2)

  return ndi.convolve(img,kernel)

# ---
median_filter = ndi.median_filter

# ---
def sharp_laplace1(img):
  """ Sharpen the image using laplacian operator

  Input:
   - img <ndarray>

  Output:
   <ndarray>
  """

  # Shapening the image with laplacian involves adding the image concolved
  # with the laplacian back to the original image. Since laplace operator
  # can generate negative values we need to use a int type image
  img = np.asarray(img, dtype=np.int)

  # Perform the operation
  sharp = img - ndi.laplace(img)

  # Clip, cast and return the result
  return np.asarray(np.clip(sharp, 0, 255), dtype=np.uint8)

# ---
def sharp_laplace2(img):
  """ Sharpen the image using laplacian operator.
  This version does shapening in the x,y directions,
  taking into account diagonal effects. Produce more
  sharpening than shap_laplace1.

  Input:
   - img <ndarray>

  Output:
   <ndarray>
  """

  # the laplacian kernel
  laplace_kern = np.array([[-1., -1., -1.],
                           [-1., 0., -1.],
                           [-1., -1., -1.]])

  # See laplace 1 above for the following lines
  img = np.asarray(img, dtype=np.int)
  sharp = img + ndi.convolve(img, laplace_kern)

  return np.asarray(np.clip(sharp, 0, 255), dtype=np.uint8)

# ---
def unsharp_mask(img, size=3):
  """ Sharpen image using the unsharp masking.

  Input:
   - img <ndarray>
   - size <int> : size of the filter window

  Output:
   <ndarray>
  """

  # apply the averaging filter
  avg = average(img, size)

  # subtract the average from the image, for a "diference" mask
  int_img = np.asarray(img, np.int)
  diff_mask = int_img - avg

  # Finally add the mask to the original image
  sharp = int_img + diff_mask

  return np.asarray(np.clip(sharp, 0, 255), dtype=np.uint8)
