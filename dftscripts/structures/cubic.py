import numpy

from pymatgen.core.structure import Structure, Lattice
from .wyckoff import wyckoff

fcc_lattice = [[.0,.5,.5],[.5,.0,.5],[.5,.5,.0]]

def ZincBlende(elements, a):
    positions = [
        wyckoff(216,'4a'),
        wyckoff(216,'4c')
    ]
    if hasattr(a,'to') and callable(getattr(a,'to')):
        a = a.to('bohr')
    lattice = Lattice(float(a) * numpy.array(fcc_lattice))
    return Structure(lattice, elements, positions).get_sorted_structure()

def Heusler(elements, phase, a):
    positions = {
        'alpha' : [ wyckoff(216,'4c'),
                    wyckoff(216,'4b'),
                    wyckoff(216,'4a'),
                    wyckoff(216,'4d')],
        'beta'  : [ wyckoff(216,'4b'),
                    wyckoff(216,'4a'),
                    wyckoff(216,'4c'),
                    wyckoff(216,'4d')],
        'gamma' : [ wyckoff(216,'4a'),
                    wyckoff(216,'4c'),
                    wyckoff(216,'4b'),
                    wyckoff(216,'4d')],
    }
    if hasattr(a,'to') and callable(getattr(a,'to')):
        a = a.to('bohr')
    lattice = Lattice(float(a) * numpy.array(fcc_lattice))
    return Structure(lattice, elements, positions[phase]).get_sorted_structure()

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