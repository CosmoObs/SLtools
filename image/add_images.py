"""Add images"""

##@package add_images

import string
import logging

def add_images(image_a, image_b, x, y):

    """ Given a common reference position (x,y) add two imagens 
    
    Input:
     
     - image_a : image array
     - image_b : image array, restriction: image_a > image_b
     - x, y : common reference position
    
    Output:
     
     - image_a """

    logging.debug('Reference position for adding images: (%d, %d)' % (x,y))

    DY2 = int(len(image_b)/2)
    if ( len(image_b)%2 ):
        minDY = y-DY2-1
    else:
        minDY = y-DY2
    
    maxDY = y+DY2

    DX2 = int(len(image_b[0])/2)
    if ( len(image_b[0])%2 ):
        minDX = x-DX2-1
    else:
        minDX = x-DX2
    maxDX = x+DX2
    
    image_a[minDY:maxDY,minDX:maxDX] = image_a[minDY:maxDY,minDX:maxDX] + image_b

    return (image_a)
