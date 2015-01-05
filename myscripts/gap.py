#!/usr/bin/env python

import numpy, pandas
from pymatgen.io.vaspio.vasp_output import Vasprun
from pymatgen.electronic_structure.core import Spin, Orbital
from pymatgen.core.periodic_table import Element

run = Vasprun('DOS/vasprun.xml')
dos = run.tdos

try:
	gap = dos.get_interpolated_gap(spin = -1)[0]
	gap = dos.get_interpolated_gap(spin = 1)[0]
except:
	pass
finally:
	print gap
