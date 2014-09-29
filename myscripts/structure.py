#!/usr/bin/env python

import os
import re
import numpy

from abipy.core.structure import Structure, Lattice

__all__ = [
    'HalfHeusler',
    'ZincBlende',
    'Cu2Sb',
]

def equation_to_elements(equation):
    return re.findall('[A-Z][^A-Z]*', equation)

fcc_lattice = [[.0,.5,.5],[.5,.0,.5],[.5,.5,.0]]

def wyckoff(spgrp, site, x = 0., y = 0., z = 0., coordinate = 0):
    coords = 230*[{}]
    coords[129] = {
    	'2a' : [[.0,.0,.0],[.5,.5,.0]],
        '2b' : [[.0,.0,.5],[.5,.5,.5]],
        '2c' : [[.0,.5,numpy.mod(z,1.0)],[.5,.0,numpy.mod(-z,1.0)]],
    }
    coords[216] = {
        '4a' : [[.00,.00,.00]],
        '4b' : [[.50,.50,.50]],
        '4c' : [[.25,.25,.25]],
        '4d' : [[.75,.75,.75]],
    }
    return coords[spgrp][site][coordinate]

def ZincBlende(elements, a):
    positions = [
        wyckoff(216,'4a'),
        wyckoff(216,'4c')
    ]
    if hasattr(a,'to') and callable(getattr(a,'to')):
        a = a.to('bohr')
    lattice = Lattice(float(a) * numpy.array(fcc_lattice))
    return Structure(lattice, elements, positions).get_sorted_structure()

def HalfHeusler(elements, phase, a):
    positions = {
        'alpha' : [ wyckoff(216,'4c'),
                    wyckoff(216,'4b'),
                    wyckoff(216,'4a')],
        'beta'  : [ wyckoff(216,'4b'),
                    wyckoff(216,'4a'),
                    wyckoff(216,'4c')],
        'gamma' : [ wyckoff(216,'4a'),
                    wyckoff(216,'4c'),
                    wyckoff(216,'4b')],
    }
    if hasattr(a,'to') and callable(getattr(a,'to')):
        a = a.to('bohr')
    lattice = Lattice(float(a) * numpy.array(fcc_lattice))
    return Structure(lattice, elements, positions[phase]).get_sorted_structure()
    
def Cu2Sb(elements, a, c, z1 = 1./3, z2 = 1./4):
    positions = [
        wyckoff(129, '2a', coordinate = 0),
        wyckoff(129, '2a', coordinate = 1),
        wyckoff(129, '2c', coordinate = 0, z = z1),
        wyckoff(129, '2c', coordinate = 1, z = z1),
        wyckoff(129, '2c', coordinate = 0, z = z2),
        wyckoff(129, '2c', coordinate = 1, z = z2),
    ]
    if hasattr(a,'to') and callable(getattr(a,'to')):
        a = a.to('ang')
        
    if hasattr(c,'to') and callable(getattr(c,'to')):
        c = c.to('ang')
        
    lattice = Lattice([[float(a),0,0],[0,float(a),0],[0,0,float(c)]])
    return Structure(lattice, elements, positions).get_sorted_structure()