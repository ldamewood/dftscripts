#!/usr/bin/env python

import numpy
import os
import logging
import re
import itertools
import matplotlib.pyplot as plt

from pymatgen.io.abinitio.workflows import Workflow
from pymatgen.io.abinitio.tasks import TaskManager, ScfTask, NscfTask
from pymatgen.io.abinitio.flows import AbinitFlow
from pymatgen.symmetry.bandstructure import HighSymmKpath

from abipy.htc.input import AbiInput
from abipy.core.constants import Energy

from myscripts.pseudos import get_psp
from myscripts.structure import HalfHeusler

scratchdir = '/p/lscratchd/damewood'
basename = os.path.dirname(os.path.abspath(__file__)).split(os.path.sep)[-1]
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
    'C': {
            'alpha': 4.961,
            'beta' : 4.912,
            'gamma': 5.139,
        },
}

ksampling = dict(
    kptopt = 1,
    nshiftk = 4,
    shiftk = [
        [.5,.5,.5],
        [.5,.0,.0],
        [.0,.5,.0],
        [.0,.0,.5],
    ],
    ngkpt = [16,16,16],
)

electrons = dict(
    ecut   = Energy(40.,'Ha'),
    occopt = 3,
    nsppol = 2,
    tsmear = Energy(1.e-3,'eV'),
    nband = 16,
    tolvrs = 1.e-10,
)

def valence(Z):
    return sum([int(re.findall(r'[\d]+',shell)[-1]) for shell in re.sub('<[^<]+?>','',Z.electronic_structure).split('.')[1:]])

def spin(structure, anti = False):
    spins = []
    up = True
    for element in structure:
        if element.specie.symbol == 'Li':
            spins.append([0.,0.,1.])
        elif element.specie.symbol == 'Mn':
            if up == True:
                spins.append([0.,0.,7.])
                if anti: up = False
            else:
                spins.append([0.,0.,-7.])
                up = True
        else:
            spins.append([0.,0.,-valence(element.specie)])
    return dict(spinat = spins)

def get_input(structure):
    pseudos = get_psp(structure)
    inp = AbiInput(pseudos)
    inp.set_variables(**structure.to_abivars())
    inp.set_variables(**ksampling)
    inp.set_variables(**electrons)
    inp.set_variables(**spin(structure))
    return inp

def dos_input(inp):
    inp = inp.deepcopy()
    inp.set_variables(**dict(
        iscf = -3,
        prtdos = 2,
        shiftk = [[.0,.0,.0],
                  [.0,.5,.5],
                  [.5,.0,.5],
                  [.5,.5,.0]],
        nshiftk = 4,
        tolvrs = 0,
        tolwfr = 1.e-15,
    ))
    return inp

def bands_input(inp):
    inp = inp.deepcopy()
    a = numpy.linspace(0,0.5,26)
    kpts = numpy.array([a,a,a-a]).T
    inp.set_variables(**dict(
        iscf = -2,
        kptopt = 0,
        nkpt = len(kpts),
        kpt = numpy.array(kpts),
        tolvrs = 0,
        tolwfr = 1.e-15,
    ))
    return inp

flow = AbinitFlow(manager = manager, workdir = workdir)
tasks = 9*[{}]
for i,Z in enumerate(['N','P']):
    phase = 'beta'
    structure_opt = HalfHeusler(['Li','Mn',Z], phase, acell_opt[Z][phase])
    structure_opt.make_supercell([[0,0,1],[1,-1,0],[1,1,-1]])
    
    inp_opt = get_input(structure_opt)
    
    work_Z = Workflow()
    scf_opt_task = work_Z.register(inp_opt, manager = manager, task_class=ScfTask)
    flow.register_work(work_Z)
    tasks[i] = {
        "opt": {
            "input" : inp_opt,
            "task": scf_opt_task,
        },
    }

for i,Z in enumerate(['N','P']):
    phase = 'beta'
    work_Nscf = Workflow()
    work_Nscf.register(dos_input(tasks[i]["opt"]["input"]),   deps = {tasks[i]["opt"]["task"]: "DEN"}, manager = manager, task_class=NscfTask)
    flow.register_work(work_Nscf)
        
flow = flow.allocate()
