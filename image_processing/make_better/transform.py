import numpy as np

std_rng = np.arange(256, dtype=np.uint8)
flt_rng = np.linspace(0., 1., 256.)
pi2_rng = np.linspace(-np.pi/2, np.pi/2, 256)

#------------------------------------
# Gamma correction
#
def power_transform(img, gamma, const=1.0):
  """Apply a power law transformation (gamma correction) to an image.
  Output is clipped to the range [0:255]

  Input:
   - img   : ndarray, ndims=2, dtype=np.uint8
   - gamma : scalar
   - const : multiplicative factor

  Output:
   - <ndarray>
  """

  # Convert to a floating point image for exponentiation:
  img = img/255.
  powimg = const * img**gamma

  # Scale and clip to the proper range:
  clipped = np.clip(255 * powimg, 0, 255)

  return np.asarray(clipped, dtype=np.uint8)

# ---
def power_transform_lut(img, gamma, const=1.0):
  """A faster version of power_transform, which uses a lookup table for assignement

  Input:
   - img   : ndarray, ndims=2, dtype=np.uint8
   - gamma : scalar
   - const : multiplicative factor

  Output:
   - <ndarray>
  """

  # Apply the transform to the floating point base range
  pow_lut = (255 * const) * flt_rng**gamma

  # Clip the lut to the range (0,255) with the proper dtype
  clipped_lut = np.asarray(np.clip(pow_lut, 0, 255), dtype=np.uint8)

  # Use the lut to get the new image values
  return clipped_lut[img]

# -

# Atan transformation
#
def _atan_lut(factor):
  """ Compute a LUT based on the arc-tangent.
  """

  # First apply the atan computation on the base LUT
  atan = np.arctan(factor * pi2_rng)

  # The scaling is done by shifting the lut so the min is zero,
  # then linearly scaling to the range [0:255]
  shifted = atan - np.min(atan)
  scaled = (255./np.max(shifted)) * shifted

  # Finally cast back to uint8 and return the final LUT
  return np.asarray(scaled, dtype=np.uint8)

def atan_transform(img, factor):
  """ Apply an arc-tangent transformation to an image.
  Useful for contrast stretching.

  Input:
   - img    :  ndarray, dtype=uint8
   - factor : scalar

  Output:
   <ndarray>
  """

  lut = _atan_lut(factor)

  return lut[img]

# Histogram Equalization
#
def equalize_hist(img):
  """ Equalize the histogram of a grayscale image.

  Input:
   - img

  Output:
   - ndarray
   """

  # First, compute the normalized histogram of the image using numpy's histogram routine.
  hist_value, bins = np.histogram(img, bins=range(257), normed=True)

  # Now, compute the cumulative distribution of the probability distribution function.
  # This can be easily with numpy's cum sum function.
  cum_dist = np.cumsum(hist_value)

  # The cumulative distribution function becomes the transformation function once we
  # scale it back to the propery range.
  transform = np.asarray(255 * cum_dist, dtype=img.dtype)

  # Transform in an array of 256 values. Think of this array as a function where the
  # indices are the domain and the values are the range. Creating a new image from this
  # function is easily done by indexing the transform array with the original image.
  equalized_img = transform[img]

  return equalized_img
