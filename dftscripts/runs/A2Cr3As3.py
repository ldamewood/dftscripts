import numpy

from pymatgen.core.periodic_table import Element
from pymatgen.core.lattice import Lattice
from pymatgen.symmetry.groups import SpaceGroup

from dftscripts.structure.wyckoff import wyckoff

a = 9.98329
b = 4.23044
lattice = Lattice.hexagonal(1.,1.*b/a)

positions = [
  wyckoff(187,'1c'),                            # Li (K2)
  wyckoff(187,'3k',x = 0.5387, coordinate = 0), # Li (K1)
  wyckoff(187,'3k',x = 0.5387, coordinate = 1), # Li (K1)
  wyckoff(187,'3k',x = 0.5387, coordinate = 2), # Li (K1)
  wyckoff(187,'3j',x = 0.0898, coordinate = 0), # Cr2
  wyckoff(187,'3j',x = 0.0898, coordinate = 1), # Cr2
  wyckoff(187,'3j',x = 0.0898, coordinate = 2), # Cr2
  wyckoff(187,'3k',x = 0.9127, coordinate = 0), # Cr1
  wyckoff(187,'3k',x = 0.9127, coordinate = 1), # Cr1
  wyckoff(187,'3k',x = 0.9127, coordinate = 2), # Cr1
  wyckoff(187,'3k',x = 0.1675, coordinate = 0), # As2
  wyckoff(187,'3k',x = 0.1675, coordinate = 1), # As2
  wyckoff(187,'3k',x = 0.1675, coordinate = 2), # As2
  wyckoff(187,'3j',x = 0.8339, coordinate = 0), # As1
  wyckoff(187,'3j',x = 0.8339, coordinate = 1), # As1
  wyckoff(187,'3j',x = 0.8339, coordinate = 2), # As1
]

with open('/Users/liam/Work/POSCAR_Li2Cr3As3','w') as f:
    f.write('(Li2Cr3As3)2\n')
    f.write(str(a))
    f.write('\n')
    f.write(str(lattice))
    f.write('\n')
    f.write('Li Cr As\n')
    f.write('4 6 6\n')
    f.write('Direct\n')
    for pos in positions:
            f.write('%+f %+f %+f\n' % tuple(numpy.array(pos)%1.0))