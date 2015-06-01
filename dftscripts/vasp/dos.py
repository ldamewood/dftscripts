#!/usr/bin/env python
from __future__ import print_function

import numpy
from pymatgen.io.vaspio.vasp_output import Vasprun
from pymatgen.electronic_structure.core import Spin


def main(vasprun='vasprun.xml', outfile='dos.csv'):
    print('NOTE: This program will only convert the total or spin polarized DOS.')
    print('Reading "%s"' % vasprun)
    run = Vasprun(vasprun)
    if run.is_spin:
        print('Extracting spin polarized DOS')
    else:
        print('Extracting total DOS')

    has_gamma = (numpy.linalg.norm(run.kpoints.kpts, axis=1) < 0.00001).any()
    if not has_gamma:
        print('WARNING: KPOINTS should include the Gamma point for improved accuracy.')
        print('         You may fix this by using a zone centered KPOINTS file, or')
        print('         if you are using the KGAMMA flag in your INCAR, set it to .TRUE.')
    
    dos = run.complete_dos
    nedos = len(run.complete_dos.energies)

    print('Fermi energy is %f. Shifting the energies so Efermi = 0.' % run.efermi)
    print('Gap might be: %f' % dos.get_interpolated_gap()[0])
    print(' -> You should verify this manually!')
    dos_array = numpy.zeros([nedos, run.is_spin + 2])
    dos_array[:, 0] = dos.energies[:] - dos.efermi
    dos_array[:, 1] = dos.densities[Spin.up]
    if run.is_spin:
        dos_array[:, 2] = dos.densities[Spin.down]
    
    header = ' energy, spin up, spin down' if run.is_spin else u' energy, dos'
    with open(outfile, 'w+') as f:
        numpy.savetxt(f, dos_array, delimiter=',', header=header)
    print('Written file as "%s"' % outfile)
