import sys
sys.path.append("/Applications/VisIt.app/Contents/Resources/2.8.0/darwin-x86_64/lib/site-packages")

from pymatgen.core.structure import Structure, Lattice
from pymatgen.core.units import Length
from pymatgen.core.periodic_table import Element
from pymatgen.io.vaspio import vasp_output
from pymatgen.transformations.standard_transformations import RotationTransformation
from math import sqrt
import numpy as np
import visit
import os
import tempfile

def make_Si001_surface(Si_layers = 4, vacuum_layers = 4):
    total_layers = Si_layers + vacuum_layers
    a = Length(5.431, "ang")
    fcc_lattice = np.array([[.0,.5,.5],[.5,.0,.5],[.5,.5,.0]])
    lattice = Lattice(fcc_lattice * a)
    surface = Structure(lattice, ['Si','Si'], [[0.00,0.00,0.00],[0.25,0.25,0.25]])
    surface.make_supercell([[0,0,1],[1,-1,0],[1,1,-1]])
    surface.make_supercell([[2,0,0],[0,2,0],[0,0,total_layers]])
    to_remove = []
    for idx,site in enumerate(surface):
        if (site.frac_coords[2] > 1.0*Si_layers/total_layers):
            to_remove.append(idx)
    surface.remove_sites(to_remove)
    print to_remove
    selective_dynamics = len(surface) * [[True,True,True]]
    for idx,site in enumerate(surface):
        if (site.frac_coords[2] < 0.001):
            selective_dynamics[idx] = [False,False,False]
            surface.replace(idx,Element('H'))
    return surface.get_sorted_structure(),np.array(selective_dynamics)

def make_InAs001_surface(inas_layers = 4, licras_layers = 4, vacuum_layers = 8):
    total_layers = inas_layers + licras_layers + vacuum_layers
    err = 1./total_layers/10
    a = Length(6.0583, "ang")
    fcc_lattice = np.array([[.0,.5,.5],[.5,.0,.5],[.5,.5,.0]])
    lattice = Lattice(fcc_lattice * a)
    surface = Structure(lattice, ['Li', 'Cr', 'As'], [[0.50,0.50,0.50],[0.00,0.00,0.00],[0.25,0.25,0.25]])
    surface.make_supercell([[0,0,1],[1,-1,0],[1,1,-1]])
    surface.make_supercell([[1,0,0],[0,1,0],[0,0,total_layers]])
    to_remove = []
    print (1.*inas_layers / total_layers)
    print 1.*(inas_layers+licras_layers)/total_layers
    print 1. - 1. / vacuum_layers
    print err
    for idx,site in enumerate(surface):
        if (site.frac_coords[2] - 1.*inas_layers / total_layers) < -err:
            if site.specie.symbol == 'Li': to_remove.append(idx)
            if site.specie.symbol == 'Cr': surface.replace(idx,Element('In'))
        if (site.frac_coords[2] - 1. + 1./vacuum_layers/8.) > -err:
            surface.replace(idx,Element('H'))
        elif (1.*(inas_layers+licras_layers)/total_layers) - site.frac_coords[2] < err:
            to_remove.append(idx)
    surface.remove_sites(to_remove)
    selective_dynamics = len(surface) * [[True,True,True]]
    for idx,site in enumerate(surface):
        if (np.linalg.norm(site.frac_coords) < err):
            selective_dynamics[idx] = [False,False,False]
    return surface.get_sorted_structure(),np.array(selective_dynamics)

def make_InAs001_slab(inas_layers = 4, licras_layers = 4):
    total_layers = inas_layers + licras_layers
    err = 1./total_layers/10
    a = Length(6.0583, "ang")
    fcc_lattice = np.array([[.0,.5,.5],[.5,.0,.5],[.5,.5,.0]])
    lattice = Lattice(fcc_lattice * a)
    surface = Structure(lattice, ['Li', 'Cr', 'As'], [[0.50,0.50,0.50],[0.00,0.00,0.00],[0.25,0.25,0.25]])
    surface.make_supercell([[0,0,1],[1,-1,0],[1,1,-1]])
    surface.make_supercell([[1,0,0],[0,1,0],[0,0,total_layers]])
    to_remove = []
    print (1.*inas_layers / total_layers)
    print 1.*(inas_layers+licras_layers)/total_layers
    print err
    for idx,site in enumerate(surface):
        if (site.frac_coords[2] - 1.*inas_layers / total_layers) < -err:
            if site.specie.symbol == 'Li': to_remove.append(idx)
            if site.specie.symbol == 'Cr': surface.replace(idx,Element('In'))
    surface.remove_sites(to_remove)
    selective_dynamics = len(surface) * [[True,True,True]]
    for idx,site in enumerate(surface):
        if (np.linalg.norm(site.frac_coords) < err):
            selective_dynamics[idx] = [False,False,False]
    return surface.get_sorted_structure(),np.array(selective_dynamics)

def make_trilayer(xyscale = 2, zscale = 4):
    a = Length(5.431,"ang")
    fcc_lattice = np.array([[.0,.5,.5],[.5,.0,.5],[.5,.5,.0]])
    lattice = Lattice(fcc_lattice * a)
    trilayer = Structure(lattice, ['Si','Si'], [[0.00,0.00,0.00],[0.25,0.25,0.25]])
    trilayer.make_supercell([[0,0,1],[1,-1,0],[1,1,-1]])
    trilayer.make_supercell([[xyscale,0,0],[0,xyscale,0],[0,0,zscale]])
    # Dope the Ga sites
    for idx,site in enumerate(trilayer):
        if np.linalg.norm(site.frac_coords - np.array([+.0,+.5,.0])) < 1e-10:
            trilayer.replace(idx,Element('Ga'))
        if np.linalg.norm(site.frac_coords - np.array([+.5,+.0,.0])) < 1e-10:
            trilayer.replace(idx,Element('Ga'))
    # Insert the Mn sites
    trilayer.append('Mn', [0.0,0.0,zscale*a-a/2], coords_are_cartesian=True)
    trilayer.append('Mn', [0,a,zscale*a-a/2], coords_are_cartesian=True)
    trilayer.append('Mn', [0.0,0.0,a/2], coords_are_cartesian=True)
    trilayer.append('Mn', [0,a,a/2], coords_are_cartesian=True)
    
    return trilayer.get_sorted_structure()

def visit_structure(structure, vdir="/Applications/VisIt.app/Contents/Resources/bin/"):
    fi, path = tempfile.mkstemp()
    pos = vasp_output.Poscar(structure)
    poscarfile = "%s.POSCAR" % path
    pos.write_file(poscarfile)
    visit.OpenDatabase(poscarfile, 0)
    visit.AddPlot("Molecule", "element", 1, 0)
    MoleculeAtts = visit.MoleculeAttributes()
    MoleculeAtts.atomSphereQuality = MoleculeAtts.Super
    MoleculeAtts.bondCylinderQuality = MoleculeAtts.Super
    visit.SetPlotOptions(MoleculeAtts)
    visit.DrawPlots()

surface,selective_dynamics = make_InAs001_slab(inas_layers = 2, licras_layers = 2)
pos = vasp_output.Poscar(surface, selective_dynamics = selective_dynamics)
pos.write_file('POSCAR')
print surface

#for i in [4, 6, 8, 10]:
#    trilayer = make_trilayer(zscale = i)
#    pos = vasp_output.Poscar(trilayer)
#    pos.write_file('POSCAR_Trilayer_z%d' % i)