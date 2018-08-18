#!/usr/bin/env python

import sys
import os
import argparse
import subprocess
from datetime import datetime
from sltools.io.rundir import Rundir

parser=argparse.ArgumentParser(description='A wrapper for PSFex')
parser.add_argument('catalog', help='Input catalog')
parser.add_argument('-c', dest='config', help='PSFex config file')
parser.add_argument('-nproc', default=1,
    help='Number of cores to use in processing')
args = parser.parse_args()

def setup(args):
    rd = Rundir('psfex')

    catalog = rd.link(args.catalog)
    config = args.config and rd.copy(args.config) or None

    os.chdir(rd.path)
    return rd.path, catalog, config

def run(catalog, config=None):
    L = ['psfex']

    if config:
        L.append('-c')
        L.append(config)
    L.append(catalog)

    try:
        subprocess.check_output(L)
    except subprocess.CalledProcessError as e:
        print 'Process error: ', e
        return

rdir, catalog, config = setup(args)
run(catalog, config)
print 'Rundir: ' + rdir
