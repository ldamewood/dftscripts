#!/usr/bin/env python2.7
from __future__ import print_function

import os
from setuptools import setup, find_packages

setup_args = dict(
    name             = 'myscripts',
    version          = '0.1',
    description      = 'Scripts I use for my work',
    author           = 'Liam Damewood',
    author_email     = 'damewood@physics.ucdavis.edu',
    package_dir      = {'myscripts' : 'myscripts'},
    packages         = ['myscripts'],
)

if __name__ == '__main__':
    setup(**setup_args)
