#!/usr/bin/env python

import sys
import os
import argparse
import subprocess
from datetime import datetime

parser=argparse.ArgumentParser(description='A wrapper for PSFex')
parser.add_argument('catalog', help='Input catalog')
parser.add_argument('-c', dest='config_path', help='SExtractor config file')
parser.add_argument('-nproc', default=1,
    help='Number of cores to use in processing')
args = parser.parse_args()

def mkrundir():
    now = datetime.now().isoformat().replace(':', '')
    path = 'modelfit-%s' % now

    os.mkdir(path)
    return path, now

def mklink(src, dst):
    os.symlink(os.path.abspath(src), dst)

def setup(args):
    rundir, now = mkrundir()

    link = os.path.basename(args.catalog)
    mklink(args.catalog, '%s/%s' % (rundir, link))

    os.chdir(rundir)
    return rundir

def run(args):
    L = ['psfex']

    if args.config_path:
        L.append('-c')
        L.append(args.config_path)
    L.append(args.catalog)

    try:
        subprocess.check_output(L)
    except subprocess.CalledProcessError as e:
        print 'Process error: ', e
        return

rundir = setup(args)
run(args)
print 'Output directory: ' + rundir
