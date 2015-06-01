import numpy

from pymatgen.core.periodic_table import Element
from abipy.abilab import Structure
from .wyckoff import wyckoff

fcc_lattice = [[.0, .5, .5], [.5, .0, .5], [.5, .5, .0]]


def ZincBlende(elements, a):
    positions = numpy.array([
        wyckoff(216, '4a'),
        wyckoff(216, '4c')
    ]).flatten()
    if hasattr(a, 'to') and callable(getattr(a, 'to')):
        a = a.to('bohr')
    return Structure.from_abivars(dict(
        natom=2,
        ntypat=len(numpy.unique(elements)),
        typat=range(len(numpy.unique(elements))),
        znucl=[Element(e).Z for e in elements],
        acell=3*[a],
        rprim=fcc_lattice,
        xred=positions,
    ))


def Heusler(elements, phase, a):
    positions = {
        'alpha': [wyckoff(216, '4c'),
                  wyckoff(216, '4b'),
                  wyckoff(216, '4a'),
                  wyckoff(216, '4d')],
        'beta': [wyckoff(216, '4b'),
                 wyckoff(216, '4a'),
                 wyckoff(216, '4c'),
                 wyckoff(216, '4d')],
        'gamma': [wyckoff(216, '4a'),
                  wyckoff(216, '4c'),
                  wyckoff(216, '4b'),
                  wyckoff(216, '4d')],
    }
    if hasattr(a, 'to') and callable(getattr(a, 'to')):
        a = a.to('bohr')
    return Structure.from_abivars(dict(
        natom=3,
        ntypat=len(numpy.unique(elements)),
        typat=range(len(numpy.unique(elements))),
        znucl=[Element(e).Z for e in elements],
        acell=3*[a],
        rprim=fcc_lattice,
        xred=positions[phase],
    ))


def HalfHeusler(elements, phase, a):
    positions = {
        'alpha': [wyckoff(216, '4c'),
                  wyckoff(216, '4b'),
                  wyckoff(216, '4a')],
        'beta': [wyckoff(216, '4b'),
                 wyckoff(216, '4a'),
                 wyckoff(216, '4c')],
        'gamma': [wyckoff(216, '4a'),
                  wyckoff(216, '4c'),
                  wyckoff(216, '4b')],
    }
    if hasattr(a, 'to') and callable(getattr(a, 'to')):
        a = a.to('bohr')
    return Structure.from_abivars(dict(
        natom=4,
        ntypat=len(numpy.unique(elements)),
        typat=range(len(numpy.unique(elements))),
        znucl=[Element(e).Z for e in elements],
        acell=3*[a],
        rprim=fcc_lattice,
        xred=positions[phase],
    ))
