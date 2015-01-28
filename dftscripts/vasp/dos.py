#!/usr/bin/env python
from __future__ import print_function

import numpy
from pymatgen.io.vaspio.vasp_output import Vasprun
from pymatgen.electronic_structure.core import Spin

def main(vasprun = 'vasprun.xml', outfile = 'dos.csv'):
    print(u'Reading "%s"' % vasprun)
    run = Vasprun(vasprun)

    has_gamma = (numpy.linalg.norm(run.kpoints.kpts,axis=1) < 0.00001).any()
    if not has_gamma:
        print('WARNING: KPOINTS should include the Gamma point for improved accuracy.')
        print('         You may fix this by using a zone centered KPOINTS file, or')
        print('         if you are using the KGAMMA flag in your INCAR, set it to .TRUE.')
    
    dos = run.complete_dos
    nedos = len(run.complete_dos.energies)

    print(u'Fermi energy is %f. Shifting the energies so Efermi = 0.' % run.efermi)
    print(u'Gap might be: %f' % dos.get_interpolated_gap()[0])
    print(u' -> You should verify this manually!')
    dos_array = numpy.zeros([nedos, run.is_spin + 2])
    dos_array[:,0] = dos.energies[:] - dos.efermi
    dos_array[:,1] = dos.densities[Spin.up]
    if run.is_spin:
        dos_array[:,2] = dos.densities[Spin.down]
    
    header = u' energy, spin up, spin down' if run.is_spin else u' energy, dos'
    with open(outfile, 'w+') as f:
        numpy.savetxt(f, dos_array, delimiter = u',', header = header)
    print(u'Written file as "%s"' % outfile)