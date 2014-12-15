#!/usr/bin/env python

import os
import re
import numpy

from pymatgen.core.structure import Structure, Lattice
from pymatgen.core.units import ArrayWithUnit

__all__ = [
    'HalfHeusler',
    'ZincBlende',
    'Cu2Sb',
    'structure_to_abivars'
]

def structure_to_abivars(self):
    """Returns a dictionary with the abinit variables."""
    types_of_specie = self.types_of_specie
    natom = self.num_sites

    znucl_type = [specie.number for specie in types_of_specie]

    znucl_atoms = self.atomic_numbers

    typat = numpy.zeros(natom, numpy.int)
    for (atm_idx, site) in enumerate(self):
        typat[atm_idx] = types_of_specie.index(site.specie) + 1

    rprim = ArrayWithUnit(self.lattice.matrix, "ang").to("bohr")
    xred = numpy.reshape([site.frac_coords for site in self], (-1,3))

    # Set small values to zero. This usually happens when the CIF file
    # does not give structure parameters with enough digits.
    #rprim = np.where(np.abs(rprim) > 1e-8, rprim, 0.0)
    #xred = np.where(np.abs(xred) > 1e-8, xred, 0.0)

    d = dict(
        natom=natom,
        ntypat=len(types_of_specie),
        typat=typat,
        xred=xred,
        znucl=znucl_type)

    d.update(dict(
        acell=3 * [1.0],
        rprim=rprim))

    #d.update(dict(
    #    acell=3 * [1.0],
    #    angdeg))

    return d

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