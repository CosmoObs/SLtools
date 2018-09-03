#!/usr/bin/env python

import sys
import os
import argparse
import astropy.io.fits as fits
import numpy as np

parser=argparse.ArgumentParser(description='...')
parser.add_argument('src', help='Source image')
parser.add_argument('mask', help='Mask image')
parser.add_argument('value', type=float, help='Value to be set')
args = parser.parse_args()

def fatal(msg):
    sys.stderr.write('error: ' + msg + '\n')
    sys.exit(1)

def getimage(path):
    try:
        hdus = fits.open(path, memmap=True)
    except IOError as e:
        fatal(str(e))

    return hdus[0].data, hdus[0].header

src, header = getimage(args.src)
mask, _ = getimage(args.mask)

if np.shape(mask) != np.shape(src):
    fatal('mask dimensions different than source dimensions')
     
indices = np.where(mask != 0)
src[indices] = args.value

out = fits.HDUList(fits.PrimaryHDU(data=src, header=header))
out.writeto(sys.stdout, checksum=True)
