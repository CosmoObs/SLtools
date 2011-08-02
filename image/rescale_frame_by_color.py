#!/usr/bin/env python
# ==================================
# Authors:
# Pedro Ferreira - pferreira@dfte.ufrn.br
# ==================================

""" Rescalonates a frame to simulate other band """

##@package rescale_frame_by_color
#
#
# Rescale counts of a frame using magzpt of both reference image and desired image and the color



# =================================================================================================
def create_color_img(frame_data, ref_magzpt, color_magzpt, color):
    """
    Rescales a given frame to simulate other band of an object in the image.

    This function uses the magnitude definition \f$ m = m_0 - 2.5log(F/F_0) \f$ to calculate the
    factor that multiplies every pixel count of the original image:
    \f$ factor = 10^{ (2/5)( magzpt_{color} - magzpt_{ref} + color ) }\f$,
    where color_magzpt and ref_magzpt are the zero point magnitudes of the output and input 
    images, while color is the magnitude difference mag_ref - mag_color (see input below).
    Multiplying an image by this factor simulates an image in the desired band for the special case 
    all objects present in the image derive from a single object.

    Input:
     - frame_data <ndarray> : numpy.ndarray(ndim=2,dtype=float). Input data of the reference frame
     - ref_magzpt   <float> : the zero point magnitude of the input image
     - color_magzpt <float> : the zero point magnitude of the image in the desired band
     - color        <float> : difference of magnitudes mag_ref - mag_color(magnitude of an object in the 
                              reference band minus the magnitude of the same object in the output image/band

    Output:
     - <ndarray> : numpy.ndarray(ndim=2,dtype=float). Color rescaled data.

    """

    factor = 10**( (2./5)*( color_magzpt - ref_magzpt + color ) ) # color_magzpt - ref_magzpt - mag_color + mag_ref

    return frame_data*factor


