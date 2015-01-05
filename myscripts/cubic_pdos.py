#!/usr/bin/env python

import numpy, pandas, sys
from pymatgen.io.vaspio.vasp_output import Vasprun
from pymatgen.electronic_structure.core import Spin, Orbital
from pymatgen.core.periodic_table import Element

if len(sys.argv) != 2:
	sys.exit('Please specify the vasprun.xml file.')
run = Vasprun(sys.argv[1])
totalDOS = run.tdos.densities
partialDOS = run.pdos
energies = run.tdos.energies - run.tdos.efermi
structure = run.structures[0]

dos = pandas.DataFrame(
  {
    'total.up': run.complete_dos.densities[Spin.up],
    'total.dn': run.complete_dos.densities[Spin.down],
  }
)
dos.index = ['%+8.4f' % e for e in energies]
for element in run.structures[0].composition.elements:
	for orbital in ['S','P','D']:
		for spin in Spin.all_spins:
			col = '%s.%s.%s' % (element, orbital.lower(), 'up' if spin == Spin.up else 'dn')
			dos[col] = run.complete_dos.get_element_spd_dos(element)[orbital].densities[spin]
for site in run.structures[0]:
	for orbital in ['t2g','e_g']:
		for spin in Spin.all_spins:
			col = '%s.%s.%s' % (site.specie, orbital, 'up' if spin == Spin.up else 'dn')
			dos[col] = run.complete_dos.get_site_t2g_eg_resolved_dos(site)[orbital].densities[spin]

print dos.to_csv(float_format = '%8.3e', index_label = 'energy')
