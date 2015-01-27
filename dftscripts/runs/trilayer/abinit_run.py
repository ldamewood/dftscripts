#!/usr/bin/env python

import os
import logging
import numpy

from pymatgen.io.abinitio.works import RelaxWork
from pymatgen.io.abinitio.tasks import TaskManager
from pymatgen.io.abinitio.flows import Flow
from pymatgen.io.abinitio.strategies import RelaxStrategy
from pymatgen.io.abinitio.abiobjects import KSampling, RelaxationMethod, AbiStructure

from pymatgen.core import Structure

from myscripts.pseudos import all_pseudos

scratchdir = '/p/lscratchd/damewood'
basename = '2014_Trilayer/abinit_6z'
workdir = os.path.join(scratchdir,basename)
logging.basicConfig()
structure = Structure.from_file('trilayer_6z.json')
manager = TaskManager.from_user_config()
ksampling = KSampling(mode='monkhorst',kpts=((6,6,2),), kpt_shifts=((0.5,0.5,0.5),(0.5,0.0,0.0),(0.0,0.5,0.0),(0.0,0.0,0.5)))
relax_ion = RelaxationMethod(ionmov = 2, optcell = 0)
relax_ioncell = RelaxationMethod(ionmov = 2, optcell = 1)

pseudos = all_pseudos()
flow = Flow(manager = manager, workdir = os.path.join(workdir, 'trilayer_6z'))
spins = numpy.zeros([len(structure),3])
spins[:4,2] = 3.
ion_input = RelaxStrategy(structure, pseudos, ksampling,
    relax_ion, accuracy="high", smearing = "fermi_dirac:0.025 eV",
    ecut = 40., pawecutdg = 80., chkprim = 0, tolmxf = 5.e-6, spinat = spins,
    restartxf = -2, nband = 60, nstep = 100)
ioncell_input = ion_input.copy()
ioncell_input.relax_algo = relax_ioncell
work = RelaxWork(ion_input, ioncell_input, manager = manager)
flow.register_work(work, workdir = 'relax')
flow = flow.allocate()