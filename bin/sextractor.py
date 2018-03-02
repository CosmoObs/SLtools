#!/usr/bin/env python

import sys
import os
import re
import datetime
import argparse

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
    for key in sorted(d):
        s = '%s\n%s %s' % (s, key, d[key])
    s = s + '\n'
    return s

def list2datstring(L):
    s = ''
    for element in sorted(L):
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

def parse_rest(args):
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

def override_config(config, overrides):
    for key in overrides:
        config[key] = overrides[key]
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


parser = argparse.ArgumentParser(description = 'A wrapper for running sextractor')
parser.add_argument('files', nargs = '+', help = 'Input files')
parser.add_argument('-p', dest = 'params_path', default = 'default.param',
    help = 'SExtractor params file')
parser.add_argument('-c', dest = 'config_path', default = 'default.sex',
    help = 'SExtractor config file')
parser.add_argument('-segobj', action = 'store_true',
    help = 'Run sextractor in SEGMENTATION mode')
parser.add_argument('-list-presets', action = 'store_true',
    help = "List available instruments' presets parameters")
parser.add_argument('-preset', default = '',
    help = "Run sextractor with given instruments' preset parameters")
parser.add_argument('-q', action = 'store_true', dest = 'quiet',
    help = 'Make SExtractor be quiet about diagnostic output.')
args = parser.parse_args()

if args.list_presets:
    inst = ''
    for i in sextractor.instrument_presets:
        inst = '%s%s ' % (inst, i)
    print 'Available instruments: %s' % inst
    sys.exit(0)

if args.segobj:
    run = sextractor.run_segobj
else:
    run = sextractor.run


infiles, overrides = parse_rest(args.files)
rundir, config, params, infiles = setup(args.config_path, args.params_path, infiles, overrides)

for infile in infiles:
    if not run(infile, params = params, args = config,
            preset = args.preset, quiet = args.quiet):
        sys.exit(1)

print 'Output directory: ' + rundir
