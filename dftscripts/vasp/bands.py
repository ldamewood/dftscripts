#!/usr/bin/env python
from __future__ import print_function

from pymatgen.io.vaspio.vasp_output import Vasprun
from pymatgen.io.vaspio.vasp_input import Incar, Kpoints
from pymatgen.symmetry.bandstructure import HighSymmKpath
from shutil import copy2 as cp, rmtree
from subprocess import check_call, CalledProcessError
from tempfile import mkdtemp
import os
import numpy
import gzip
import argparse

def _check(filename):
    """ This routine raises an exception when the file cannot be found. """
    if not os.path.isfile(filename):
        raise Exception("%s not found" % filename)

def _automatic_kpoints(density, ibz):
    """ Automatic line-mode with specified density and BZ """
    kpoints = list()
    labels = list()
    for path in ibz.kpath["path"]:
        kpoints.append(ibz.kpath["kpoints"][path[0]])
        labels.append(path[0])
        for i in range(1, len(path) - 1):
            kpoints.append(ibz.kpath["kpoints"][path[i]])
            labels.append(path[i])
            kpoints.append(ibz.kpath["kpoints"][path[i]])
            labels.append(path[i])

        kpoints.append(ibz.kpath["kpoints"][path[-1]])
        labels.append(path[-1])

    return Kpoints("Line_mode KPOINTS file",
                    style=Kpoints.supported_modes.Line_mode,
                    coord_type="Reciprocal",
                    kpts=kpoints,
                    labels=labels,
                    num_kpts=int(density))

def _clean_exit(original_path, temp_path, err = 1):
    """ Cleanly quit in case of error """
    os.chdir(original_path)
    rmtree(temp_path)
    exit(err)

def main():
    """ Main routine """
    parser = argparse.ArgumentParser(description='Calculate the band structure of a vasp calculation.')
    parser.add_argument('-d', '--density', nargs='?', default=10, type=int,
                        help = 'k-point density of the bands (default: 10)')
    parser.add_argument('-s', '--sigma', nargs='?', default=0.04, type=float,
                        help = 'SIGMA value in eV (default: 0.04)')
    parser.add_argument('-l', '--linemode', nargs='?', default=None, type=str,
                        help='linemode file')
    args = parser.parse_args()
    
    line_density = args.density
    line_file = args.linemode

    # Check that required files exist
    _check(u'vasprun.xml')
    _check(u'INCAR')
    _check(u'POSCAR')
    _check(u'POTCAR')
    _check(u'CHGCAR')

    # Get IBZ from vasprun
    vasprun = Vasprun(u'vasprun.xml')
    ibz = HighSymmKpath(vasprun.final_structure)

    # Create a temp directory
    print('Create temporary directory ... ', end='')
    tempdir = mkdtemp()
    print(tempdir)

    # Edit the inputs
    print('Saving new inputs in temporary directory ... ')
    incar = Incar.from_file(u'INCAR')
    print('Making the following changes to the INCAR:')
    print('  ICHARG = 11')
    print('  ISMEAR = 0')
    print('  SIGMA = %f' % args.sigma)
    incar[u'ICHARG'] = 11  # Constant density
    incar[u'ISMEAR'] = 0 # Gaussian Smearing
    incar[u'SIGMA'] = args.sigma # Smearing temperature
    incar.write_file(os.path.join(tempdir, u'INCAR'))
    # Generate line-mode kpoint file

    if line_file is None:
        print('Creating a new KPOINTS file:')
        kpoints = _automatic_kpoints(line_density, ibz)
        kpoints.write_file(os.path.join(tempdir, u'KPOINTS'))
        print('### BEGIN KPOINTS')
        print(kpoints)
        print('### END KPOINTS')
    else:
        cp(line_file, os.path.join(tempdir, u'KPOINTS'))
    
    # Copy other files (May take some time...)
    print('Copying POSCAR, POTCAR and CHGCAR to the temporary directory.')
    cp(u'POSCAR', os.path.join(tempdir, u'POSCAR'))
    cp(u'POTCAR', os.path.join(tempdir, u'POTCAR'))
    cp(u'CHGCAR', os.path.join(tempdir, u'CHGCAR'))

    # cd to temp directory and run vasp
    path = os.getcwd()
    os.chdir(tempdir)
    print('Running VASP in the temporary directory ...')
    try:
        check_call(u'vasp')
    except CalledProcessError:
        print('There was an error running VASP')
        _clean_exit(path, tempdir)

    # Read output
    vasprun = Vasprun('vasprun.xml')
    ibz = HighSymmKpath(vasprun.final_structure)
    kpath, labels = ibz.get_kpoints(line_density)
    bands = vasprun.get_band_structure()
    print('Success! Efermi = %f' % bands.efermi)
    
    # Backup vasprun.xml only
    print('Making a gzip backup of vasprun.xml called bands_vasprun.xml.gz')
    zfile = os.path.join(path, 'bands_vasprun.xml.gz')
    try:
        with gzip.open(zfile, 'wb') as gz, open('vasprun.xml', 'rb') as vr:
            gz.writelines(vr)
    except Exception:
        print('There was an error with gzip')
        _clean_exit(path, tempdir)

    # Return to original path
    os.chdir(path)

    # Write band structure
    
    # There may be multiple bands due to spin or noncollinear.
    for key, item in bands.bands.items():
        print('Preparing bands_%s.csv' % key)
        
        kpts = []

        # Get list of kpoints
        for kpt in bands.kpoints:
            kpts.append(kpt.cart_coords)
        
        # Subtract fermi energy
        print('Shifting energies so Efermi = 0.')
        band = numpy.array(item) - bands.efermi

        # Prepend kpoint vector as label (TODO: get real label Gamma, etc.)
        out = numpy.hstack([kpts, band.T])
        
        final = []
        lastrow = float('inf') * numpy.ones(3)
        for row in out:
            if numpy.linalg.norm(row[:3] - lastrow) > 1.e-12:
                final.append(row)
            lastrow = row[:3]

        # Write bands to csv file.
        print('Writing bands_%s.csv to disk' % key)
        with open('bands_%s.csv' % key, 'w+') as f:
            numpy.savetxt(f, numpy.array(final), delimiter = ',', header = ' kptx, kpty, kptz, band1, band2, ... ')

    # Delete temporary directory
    _clean_exit(path, tempdir, 0)
