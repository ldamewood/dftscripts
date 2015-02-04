from .wyckoff import wyckoff
from pymatgen.core.structure import Structure, Lattice

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