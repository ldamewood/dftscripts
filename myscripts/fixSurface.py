#!/usr/bin/env python

import numpy

from pymatgen.io.vaspio import Poscar
from pymatgen.core import Element

si = Poscar.from_file('/Users/liam/Downloads/Si.POSCAR')
s = si.structure

s.remove_species("H")

for element in s:
    if element.z < 11.4 and element.specie.symbol == 'Si':
        coords = element.coords + numpy.array([0,1.19207,-0.842944])
        s.append('H',coords,coords_are_cartesian=True)
        coords = element.coords + numpy.array([0,-1.19207,-0.842944])
        s.append('H',coords,coords_are_cartesian=True)
