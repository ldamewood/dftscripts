#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function, division

def main(filename = 'OSZICAR'):
    with open(filename, 'r') as oszicar:
        lines = oszicar.readlines()
    
    moment = float(lines[-1].split()[-1])
    
    for i in range(len(lines) - 2,1,-1):
        if len(lines[i].split()) == 7:
            error = float(lines[i].split()[-1])
            break
    
    print('(%f±%f)μʙ' % (moment, 2. * error))
