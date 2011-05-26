import random
import numpy as np
import scipy as sc
import scipy.ndimage as ndi

# ---
def gen_snp_noise(img, perc=10):
  """ Generate salt-and-pepper noise in an image. Salt-and-Pepper noise is
  defined as randomly dispersed values of 0 and 255.

  Input:
   - img
   - perc
   
  Output:
   - <ndarray>
  """

  # Create a flat copy of the image
  flat = img.ravel().copy

  # The total number of pixels
  total = len(flat)

  # The number of pixels we need to modify
  nmod = int(total * (perc/100.))

  # Random indices to modify
  indices = np.random.random_integers(0, total-1, (nmod,))

  # Set the first half of the pixels to 0
  flat[indices[:nmod/2]] = 0

  # Set the second half of the pixels to 255
  flat[indices[nmod/2:]] = 255

  return flat.reshape(img.shape)

# ---
def gen_gaussian_noise(img, stdev=3.):
  """ Generate zero-mean additive gaussian noise

  Input:
   - img <ndarray>
   - stdev <float>

  Output:
   <ndarray>
  """

  noise = np.random.normal(0., stdev, img.shape)
  noisy_img = (1. * img) + noise

  return np.asarray(np.clip(noisy_img, 0, 255), dtype=np.uint8)

