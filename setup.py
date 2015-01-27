#!/usr/bin/env python2.7
from __future__ import print_function

import glob
import os
from setuptools import setup, find_packages

SETUP_PTH = os.path.dirname(os.path.abspath(__file__))

setup_args = dict(
    name             = 'dftscripts',
    version          = '0.1',
    description      = 'Scripts on top of pymatgen',
    author           = 'Liam Damewood',
    author_email     = 'damewood@physics.ucdavis.edu',
    package_dir      = {'dftscripts' : 'dftscripts'},
    packages         = ['dftscripts'],
    scripts          = glob.glob(os.path.join(SETUP_PTH, "scripts", "*"))
)

if __name__ == '__main__':
    setup(**setup_args)
