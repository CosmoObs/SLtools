#!/usr/bin/env python

import sys
import sltools
import pyfits
import sltools.image.imcp as imcp

from optparse import OptionParser

usage = "%prog [-o output] objects segmentation object-id"
parser = OptionParser(usage = usage)
parser.add_option('-o', dest = 'output', help = "Output filename")
(opts, args) = parser.parse_args()

if len(args) < 3:
    parser.print_help()
    sys.exit(1)

def getimage(path):
    try:
        hdulist = pyfits.open(path, memmap = True)
    except IOError as e:
        sys.stderr.write(str(e) + '\n')
        sys.exit(1)

    return hdulist[0].data, hdulist[0].header


obj, hdr = getimage(args[0])
seg, hdr = getimage(args[1])
id = int(args[2])

pstamp, hdr = imcp.segstamp(seg, id, hdr = hdr, objimg = obj)

outfile = opts.output or ('ps-%d.fits' % id)
pyfits.writeto(outfile, pstamp, hdr)

print 'Output file: %s' % outfile
