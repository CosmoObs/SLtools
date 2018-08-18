#!/usr/bin/env python

import sys
import os
import re
import datetime
import argparse
from multiprocessing import Process
from sltools.image import sextractor
from sltools.io.rundir import Rundir

#
# path is a file in the following format:
#
# KEY1  VAL11, VAL12  # comment 1
# KEY2  VAL21, VAL22  # comment 2
#
def dict_read(file):
    d = {}
    for line in file:
        m = re.match('^(\w+)\s+([^#]+)', line)
        if not m: continue

        d.update({m.group(1) : m.group(2)})
    return d

def dict_format(d):
    s = ''
    for key in sorted(d):
        s = '%s%s %s\n' % (s, key, d[key])
    return s

def list_format(L):
    s = ''
    for element in sorted(L):
        s = '%s%s\n' % (s, element)
    return s

def save(data, path, formatter):
    file = open(path, 'w')
    file.write(formatter(data))

def getconfig(path):
    try:
        file = open(path)
    except:
        file = []
    return dict_read(file)

def parse_overrides(args):
    infiles = []
    overrides = {}

    i = 0
    while i < len(args):
        arg = args[i]
        if arg[0] == '-':
            key = args[i][1:]
            value = args[i+1]
            print key, value
            overrides[key] = value
            i = i + 1
        else:
            infiles.append(arg)
        i = i + 1

    return infiles, overrides

def parse_params(path):
    return path and [ line.rstrip('\n') for line in open(path).readlines() ] or []

def override_config(config, overrides):
    for key in overrides:
        config[key] = overrides[key]
    return config


def setup(config_path, params_path, images, cmdline_config,
        psf=None):
    rd = Rundir('SErun')

    config = getconfig(config_path)
    config = override_config(config, cmdline_config)
    params = parse_params(params_path)

    newimages = []
    for image in images:
        newpath = rd.link(image)
        newimages.append(newpath)

    rd.copy('default.conv')
    psf = psf and rd.copy(psf)

    os.chdir(rd.path)

    save(config, os.path.basename('default.sex'), dict_format)
    save(params, os.path.basename('default.param'), list_format)

    return rd.path, config, params, newimages, psf

def process_file(infile, func = None, params = None, args = None,
    preset = None, quiet = False):

    return func(infile, params = params, args = args,
            preset = preset, quiet = quiet)

parser=argparse.ArgumentParser(description='A wrapper for sextractor')
parser.add_argument('files', nargs= '+', help='Input files')
parser.add_argument('-p', dest='params_path', default='default.param',
    help='SExtractor params file')
parser.add_argument('-c', dest='config_path', default='default.sex',
    help='SExtractor config file')
parser.add_argument('-psf', dest='psf',
    help='PSF model from PSFex')
parser.add_argument('-nproc', type=int, help='Number of processors')
parser.add_argument('-segobj', action='store_true',
    help='Run sextractor in SEGMENTATION mode')
parser.add_argument('-list-presets', action='store_true',
    help="List available instruments' presets parameters")
parser.add_argument('-preset', default='',
    help="Run sextractor with given instruments' preset parameters")
parser.add_argument('-q', action='store_true', dest='quiet',
    help='Make SExtractor be quiet about diagnostic output.')
args=parser.parse_args()

if args.list_presets:
    inst = ''
    for i in sextractor.instrument_presets:
        inst = '%s%s ' % (inst, i)
    print 'Available instruments: %s' % inst
    sys.exit(0)

infiles, overrides = parse_overrides(args.files)
rundir, config, params, infiles, psf = setup(args.config_path, args.params_path,
    infiles, overrides, psf=args.psf)

if args.segobj:
    run = sextractor.run_segobj
else:
    run = sextractor.run

nproc = args.nproc or len(infiles)
nblock = len(infiles) / nproc
if (len(infiles) % nproc) > 0:
    nblock = nblock + 1

for i in range(0, nblock):
    start = i * nproc
    end = (i+1) * nproc
    procs = []

    for infile in infiles[start:end]:
        lastdot = infile.rfind('.fits')
        catalog_name = infile[:lastdot] + '_cat.' + infile[lastdot+1:]
        config.update({'CATALOG_TYPE': 'FITS_LDAC', 'CATALOG_NAME': catalog_name})

        if psf:
            config.update({'PSF_NAME': psf})

        p = Process(target=process_file, args=(infile, run, params, config,
                args.preset, args.quiet))
        procs.append(p)
        p.start()

    for p in procs:
        p.join()

print 'Output directory: ' + rundir
