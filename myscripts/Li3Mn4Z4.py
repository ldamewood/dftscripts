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
            'Li1' : {
                'alpha': 4.961,
                'beta' : 4.404,
                'gamma': 4.996,
            },
            'Li2' : {
                'alpha': 4.500,
                'beta' : 4.500,
                'gamma': 4.500,
            },
            'Li3' : {
                'alpha': 4.552,
                'beta' : 4.552,
                'gamma': 4.552,
            },
         },
    'P': {
            'Li1' : {
                'alpha': 5.300,
                'beta' : 5.422,
                'gamma': 5.260,
            },
            'Li2' : {
                'alpha': 5.220,
                'beta' : 5.220,
                'gamma': 5.220,
            },
            'Li3' : {
                'alpha': 5.502,
                'beta' : 5.502,
                'gamma': 5.502,
            },
         },
    'Si': {
            'Li1' : {
                'alpha': 4.894,
                'beta' : 5.527,
                'gamma': 5.517,
            },
            'Li2' : {
                'alpha': 4.894,
                'beta' : 5.527,
                'gamma': 5.517,
            },
            'Li3' : {
                'alpha': 5.629,
                'beta' : 5.778,
                'gamma': 5.788,
            },
         },
}

ksampling = KSampling(mode='monkhorst',kpts=((12,12,12),), kpt_shifts=((0.5,0.5,0.5),(0.5,0.0,0.0),(0.0,0.5,0.0),(0.0,0.0,0.5)))
relax_ion = RelaxationMethod(ionmov = 2, optcell = 0)
relax_ioncell = RelaxationMethod(ionmov = 2, optcell = 1)

pseudos = all_pseudos()
flows = []
for (Z, i) in itertools.product(['P','N','Si'],range(1,4)):
    x = 4 - i
    name = u'Li%dMn4%s4' % (x,Z)
    print(name)
    flow = Flow(manager = manager, workdir = os.path.join(workdir, name))
    li = u'Li%d' % x
    for (phase) in ["alpha","beta","gamma"]:
        structure = HalfHeusler(['Li','Mn',Z], phase, acell_opt[Z][li][phase])
        structure.make_supercell([[-1,1,1],[1,-1,1],[1,1,-1]])
        structure.remove_sites(list(range(i)))
        structure.sort(key=lambda j: j.specie.Z)
        assert name == structure.formula.replace(' ','')
        spins = numpy.zeros([len(structure),3])
        for j,atom in enumerate(structure):
            if atom.specie.symbol == 'Li':
                spins[j,2] = 1.
            if atom.specie.symbol == 'Mn':
                spins[j,2] = 3.
        ion_input = RelaxStrategy(AbiStructure(structure), pseudos, ksampling,
            relax_ion, accuracy="high", smearing = "fermi_dirac:0.025 eV",
            ecut = 40., pawecutdg = 80., chkprim = 0, tolmxf = 5.e-6,
            spinat = spins, restartxf = -2, nband = 60, nstep = 100)
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

def get_status():
    for flow in flows:
        for work in flow:
            for task in work:
                try:
                    task.check_status()
                except TypeError:
                    pass

#build_and_pickle_dump()