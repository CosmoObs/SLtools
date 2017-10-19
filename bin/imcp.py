#!/usr/bin/env python

import sys
import os
import pyfits
import string
import re

from optparse import OptionParser
from sltools.io import log
from sltools.image import imcp

usage = "%prog [options] image.fits [x y]"
parser = OptionParser(usage = usage)

parser.add_option('--debug',
    action = 'store_true', dest = 'debug', default = False,
    help = 'Enable debug logging')

parser.add_option('-v',
    action = 'store_true', dest = 'verbose', default = False,
    help = 'Enable verbose logging')

parser.add_option('--coord-degrees',
    action = 'store_true', dest = 'is_coord_deg', default = False,
    help = "If given, 'x,y' arguments (i.e, central coordinates) are in degrees. If not given (default), coordinates are expected to be in pixels")

parser.add_option('--shape-degrees',
    action = 'store_true', dest = 'is_shape_deg', default = False,
    help = "If given, cutout 'shape' (see option '-s') are in degrees. If not given (default), cutout shape is expected to be in pixels.")

parser.add_option('-s',
    dest = 'shape', default = '250,250',
    help = "Output image shape. Default is 250,250")

# XXX is this option really used?
parser.add_option('--increase',
    dest = 'incr', default = 0,
    help = "Amount of border (enlarge) of object poststamp. It is an additive or multiplicative factor. See 'relative_increase' flag")

# XXX is this option really used?
parser.add_option('--relative_increase',
    action = 'store_true', dest = 'rel_incr', default = False,
    help = "Turn on the multiplicative factor for increase image, otherwise, 'increase' will be treated as an absolute value.")

parser.add_option('--prefix',
    dest = 'prefix', default = '',
    help = "Output filename prefix")

parser.add_option('--no-header',
    action = 'store_false', dest = 'use_header', default = True,
    help = "Do not use input image header.")

parser.add_option('--segimg',
    dest = 'segimg', default = None,
    help = "Segmented image (e.g, SE's SEGMENTATION output)")

parser.add_option('--objIDs',
    dest = 'IDs', default = '',
    help = "Comma-separated list of object IDs found on 'segimg' to extract")

parser.add_option('--hdu-n', default = 0, help = "HDU position")

(opts, args) = parser.parse_args()

logging = log.init(debug = opts.debug, verbose = opts.verbose)

nrequired = opts.segimg and 1 or 3
if len(args) < nrequired:
    parser.print_help()
    sys.exit(1)

IDs = opts.IDs
prefix = opts.prefix

infits = args[0]

coord_unit = opts.is_coord_deg and 'degrees' or 'pixel'
size_unit = opts.is_shape_deg and 'degrees' or 'pixel'

hdulist = pyfits.open(infits, memmap = True)

if opts.use_header:
    image = hdulist[int(opts.hdu_n)].data
    hdr = hdulist[int(opts.hdu_n)].header
else:
    image = hdulist[int(opts.hdu_n)].data
    hdr = None

if opts.segimg:
    IDs = [ int(ID) for ID in re.split(',', opts.IDs) ]

    seghdu = pyfits.open(opts.segimg, memmap = True)
    seg = seghdu[0].data
    hdr = seghdu[0].header
    print(opts.segimg)

    if prefix == '':
        prefix = 'pstamp'

    for i in range(len(IDs)):
      outpath = '%s_%s.fits' % (prefix, IDs[i])
      print opts.segimg, infits

      image, header = imcp.segstamp(seg, IDs[i], objimg = image,
          hdr = hdr, increase = opts.incr, relative_increase = opts.rel_incr)

      try:
          pyfits.writeto(outpath, image, header)
      except IOError as e:
          sys.stderr.write('error: %s\n' % e)
          sys.exit(1)

else:
    if len(args) < 3:
        parser.print_help()
        sys.exit(1)

    x = args[1]
    y = args[2]
    dx, dy = string.split(opts.shape, sep = ',')

    image, header = imcp.cutout(image, hdr = hdr, coord_unit = coord_unit,
        xc = x, yc = y, size_unit = size_unit, x_size = dx, y_size = dy)
    if prefix == '':
        prefix = 'cut'
    outpath = '%s.fits' % prefix

    try:
        pyfits.writeto(outpath, image, header)
    except IOError as e:
        sys.stderr.write('error: %s\n' % e)
        sys.exit(0)

sys.stderr.write('Output image: %s\n' % outpath)
