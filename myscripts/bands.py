#!/usr/bin/env python

from pymatgen.io.vaspio.vasp_output import Vasprun
import numpy

run = Vasprun('BAND/vasprun.xml')
dos = Vasprun('DOS/vasprun.xml')
allbands = run.get_band_structure('BAND/IBZKPTS', dos.efermi, line_mode = False)
prevkpt = numpy.array([5,5,5])
skip = [False] * len(allbands.kpoints)
for i,kpt in enumerate(allbands.kpoints):
	if numpy.linalg.norm(kpt.frac_coords - prevkpt) < 1.e-5:
		skip[i] = True
	prevkpt = kpt.frac_coords
# Spin up
bands = numpy.array(allbands.bands[1])
for band in bands:
	ikpts = list(range(len(band) - sum(skip)))
	band = numpy.array(band) - dos.efermi
        band = band[numpy.array(skip) == False]
	str = ' '.join(['(%d,%5.3f)' % (i,j) for (i,j) in numpy.array([ikpts,band]).T])
	print '\\addplot[spin up] coordinates {%s};' % str
# Spin dn
bands = numpy.array(allbands.bands[-1])
for band in bands:
	ikpts = list(range(len(band) - sum(skip)))
	band = numpy.array(band) - dos.efermi
        band = band[numpy.array(skip) == False]
	str = ' '.join(['(%d,%5.3f)' % (i,j) for (i,j) in numpy.array([ikpts,band]).T])
	print '\\addplot[spin dn] coordinates {%s};' % str
