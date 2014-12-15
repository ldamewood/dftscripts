#!/usr/bin/env python

import numpy
import os
import logging
import re
import itertools
import matplotlib.pyplot as plt

from pymatgen.io.abinitio.works import RelaxWork
from pymatgen.io.abinitio.tasks import TaskManager
from pymatgen.io.abinitio.flows import Flow
from pymatgen.io.abinitio.strategies import RelaxStrategy
from pymatgen.io.abinitio.abiobjects import KSampling, RelaxationMethod, AbiStructure
from pymatgen.core.units import Energy

from myscripts.pseudos import get_psp
from myscripts.structure import HalfHeusler

scratchdir = '/p/lscratchd/damewood'
basename = 'LixMn4Z4'
workdir = os.path.join(scratchdir,basename)
logging.basicConfig()
manager = TaskManager.from_user_config()

acell_opt = {
    'N': {
            'alpha': 4.961,
            'beta' : 4.912,
            'gamma': 5.139,
         },
    'P': {
            'alpha': 5.600,
            'beta' : 5.717,
            'gamma': 5.715,
         },
    'Si': {
            'alpha': 5.629,
            'beta' : 5.778,
            'gamma': 5.788,
         },
}

ksampling = KSampling(mode='monkhorst',kpts=((11,11,11),), kpt_shifts=((0.5,0.5,0.5),(0.5,0.0,0.0),(0.0,0.5,0.0),(0.0,0.0,0.5)))
relax_ion = RelaxationMethod(ionmov = 3, optcell = 0)
relax_ioncell = RelaxationMethod(ionmov = 3, optcell = 1, tolmxf = 5e-06)

flows = []
for (Z, i) in itertools.product(["P","N","Si"],range(4)):
    x = 4 - i
    name = 'Li%dMn4%s4' % (x,Z)
    print(name)
    flow = Flow(manager = manager, workdir = os.path.join(workdir, name))
    for (phase) in ["alpha","beta","gamma"]:
        structure = HalfHeusler(['Li','Mn',Z], phase, acell_opt[Z][phase])
        structure.make_supercell([[-1,1,1],[1,-1,1],[1,1,-1]])
        to_remove = i
        while(to_remove>0):
            for atom in structure:
                if atom.specie.symbol == u'Li':
                    structure.remove(atom)
                    to_remove -= 1
        ion_input = RelaxStrategy(AbiStructure(structure), get_psp(structure), ksampling,
            relax_ion, accuracy="high", smearing = "fermi_dirac:0.025 eV",
            ecut = Energy(40., "eV"), pawecutdg = Energy(80., "eV"))
        ioncell_input = ion_input.copy()
        ioncell_input.relax_algo = relax_ioncell
        work = RelaxWork(ion_input, ioncell_input, manager = manager)
        flow.register_work(work, workdir = phase)
    flow = flow.allocate()
    flows.append(flow)