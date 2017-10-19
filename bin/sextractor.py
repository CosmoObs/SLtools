#!/usr/bin/env python

import sys
import os
import pyfits

from optparse import OptionParser
from sltools.io import config_parser as CP
from sltools.image import sextractor

usage = '%prog [options] image.fits'
parser = OptionParser(usage = usage)

parser.add_option('-p',
    dest = 'params_file', default = None,
    help = 'Sextractor params file (e.g. default.param)')

parser.add_option('-c',
    dest = 'config_file', default = None,
    help = 'Sextractor config file (e.g. default.sex)')

parser.add_option('--segment',
    action = 'store_true', dest = 'run_seg', default = False,
    help = 'Run sextractor in SEGMENTATION mode')

parser.add_option('--list-presets',
    action = 'store_true', dest = 'list_presets', default = False,
    help = "List available instruments' presets parameters")

parser.add_option('--preset',
    dest = 'preset', default = '',
    help = "Run sextractor with given instruments' preset parameters")

parser.add_option('-o',
    dest = 'temp_dir', default = '',
    help = 'Temporary directory for run-time files.')

parser.add_option('-q',
    action = 'store_true', dest = 'quiet', default = False,
    help = 'Make SExtractor be quiet about diagnostic output.')
    
(opts, args) = parser.parse_args()

if opts.list_presets:
    inst = ''
    for i in sextractor.instrument_presets:
        inst = '%s%s ' % (inst, i)
    print 'Available instruments: %s' % inst
    sys.exit(0)

if len(args) != 1:
    parser.print_usage()
    print "Try '--help' for more info."
    sys.exit(0)

infits = args[0]
config = opts.config_file and read_config(opts.config_file) or {}
params = opts.params_file and [ line.rstrip('\n') for line in open(opts.params_file).readlines() ] or []

if opts.run_seg:
    out = sextractor.run_segobj(infits, params = params, args = config,
               preset = opts.preset, temp_dir = opts.temp_dir, 
               quiet = opts.quiet)
    if not out:
        sys.exit(1)

    print 'Output files'
    print 'CATALOG: %s' % (out['CATALOG'])
    print 'OBJECTS: %s' % (out['OBJECTS'])
    print 'SEGMENT: %s' % (out['SEGMENTATION'])
else:
    if not sextractor.run(infits, params = params, args = config, preset = opts.preset,
              temp_dir = opts.temp_dir, quiet = opts.quiet):
        sys.exit(1)
