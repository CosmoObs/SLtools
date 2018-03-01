#!/usr/bin/env python

import sys
import os
import re
import datetime
import shutil

from optparse import OptionParser
from sltools.image import sextractor

#
# path is a file in the following format:
#
# KEY1  VAL11, VAL12  # comment 1
# KEY2  VAL21, VAL22  # comment 2
#
def dat2dict(file):
    d = {}
    for line in file:
        m = re.match('^(\w+)\s+([^#]+)', line)
        if not m: continue

        d.update({m.group(1) : m.group(2)})
    return d

def dict2datstring(d):
    s = ''
    for key in d:
        s = '%s\n%s %s' % (s, key, d[key])
    s = s + '\n'
    return s

def list2datstring(L):
    s = ''
    for element in L:
        s = '%s\n%s' % (s, element)
    s = s + '\n'
    return s

def dictsave(d, path):
    file = open(path, 'w')
    file.write(dict2datstring(d))

def listsave(L, path):
    file = open(path, 'w')
    file.write(list2datstring(L))

def getconfig(path):
    try:
        file = open(path)
    except:
        file = []
    return dat2dict(file)

def override_config(config, overrides):
    for i in range(0, len(overrides) - 1):
        key = overrides[i][1:]    # skip leading -
        value = overrides[i+1]
        config[key] = value
    return config

def create_rundir():
    now = datetime.datetime.now().isoformat().replace(':', '')
    path = 'SErun-%s' % now

    os.mkdir(path)
    return path, now

def create_status(now, config_path, params_path, images_paths):
    status = open('status.txt', 'w')
    status.write('# date,config file,params file,input files\n')

    for image_path in images_paths:
        images = '%s,' % image_path

    status.write(now + ',' + config_path + ',' + params_path + ',' + images + '\n')
    status.close()

def absolute_link(src, dst):
    os.symlink(os.path.abspath(src), dst)

def setup(config_path, params_path, images_paths, cmdline_config):
    rundir, now = create_rundir()

    config = getconfig(config_path)
    config = override_config(config, cmdline_config)

    params = params_path and [ line.rstrip('\n') for line in open(params_path).readlines() ] or []

    new_images_paths = []
    for image_path in images_paths:
        link_path = os.path.basename(image_path)
        absolute_link(image_path, '%s/%s' % (rundir, link_path))
        new_images_paths.append(link_path)

    absolute_link('default.conv', '%s/default.conv' % rundir)

    os.chdir(rundir)

    dictsave(config, os.path.basename(config_path))
    listsave(params, os.path.basename(params_path))

    create_status(now, config_path, params_path, new_images_paths)
    return rundir, config, params, new_images_paths

usage = '%prog [options] image.fits [-- -OPTION1 VALUE1 ... -OPTIONK VALUEK]'
parser = OptionParser(usage = usage)

parser.add_option('-p',
    dest = 'params_path', default = 'default.param',
    help = 'Sextractor params file (e.g. default.param)')

parser.add_option('-c',
    dest = 'config_path', default = 'default.sex',
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

if len(args) < 1:
    parser.print_usage()
    print "Try '--help' for more info."
    sys.exit(0)


infiles = []
overrides = []

i = 0
while i < len(args):
    arg = args[i]
    if arg[0] == '-':
        overrides.append(args[i])
        overrides.append(args[i+1])
        i = i + 1
    else:
        infiles.append(arg)
    i = i + 1

print infiles
print overrides
rundir, config, params, inpaths = setup(opts.config_path, opts.params_path, infiles, overrides)

for infits in inpaths:
    if opts.run_seg:
        out = sextractor.run_segobj(infits, params = params, args = config,
                   preset = opts.preset, quiet = opts.quiet)
        if not out:
            sys.exit(1)

        print 'Output files'
        print 'CATALOG: %s' % (out['CATALOG'])
        print 'OBJECTS: %s' % (out['OBJECTS'])
        print 'SEGMENT: %s' % (out['SEGMENTATION'])
    else:
        if not sextractor.run(infits, params = params, args = config, preset = opts.preset,
            quiet = opts.quiet):
            sys.exit(1)

print 'Output directory: ' + rundir
