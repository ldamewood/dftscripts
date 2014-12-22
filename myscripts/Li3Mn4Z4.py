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

from myscripts.pseudos import all_pseudos
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

ksampling = KSampling(mode='monkhorst',kpts=((12,12,12),), kpt_shifts=((0.5,0.5,0.5),(0.5,0.0,0.0),(0.0,0.5,0.0),(0.0,0.0,0.5)))
relax_ion = RelaxationMethod(ionmov = 2, optcell = 0)
relax_ioncell = RelaxationMethod(ionmov = 2, optcell = 1)

pseudos = all_pseudos()
flows = []
for (Z, i) in itertools.product(['P','N','Si'],range(4)):
    x = 4 - i
    name = u'Li%d Mn4 %s4' % (x,Z)
    print(name)
    flow = Flow(manager = manager, workdir = os.path.join(workdir, name))
    for (phase) in ["alpha","beta","gamma"]:
        structure = HalfHeusler(['Li','Mn',Z], phase, acell_opt[Z][phase])
        structure.make_supercell([[-1,1,1],[1,-1,1],[1,1,-1]])
        structure.remove_sites(list(range(i)))
        structure.sort(key=lambda j: j.specie.Z)
        assert name == structure.formula
        spins = numpy.zeros([len(structure),3])
        for j,atom in enumerate(structure):
            if atom.specie.symbol == 'Li':
                spins[j,2] = 1.
            if atom.specie.symbol == 'Mn':
                spins[j,2] = 3.
        ion_input = RelaxStrategy(AbiStructure(structure), pseudos, ksampling,
            relax_ion, accuracy="high", smearing = "fermi_dirac:0.025 eV",
            ecut = 40., pawecutdg = 80., chkprim = 0, tolmxf = 5.e-6,
            spinat = spins)
        ioncell_input = ion_input.copy()
        ioncell_input.relax_algo = relax_ioncell
        work = RelaxWork(ion_input, ioncell_input, manager = manager)
        flow.register_work(work, workdir = phase)
    flow = flow.allocate()
    flows.append(flow)

def build_and_pickle_dump():
    for flow in flows:
        flow.build_and_pickle_dump()

def rapidfire():
    for flow in flows:
        flow.rapidfire()

#build_and_pickle_dump()