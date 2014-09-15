#!/usr/bin/env python

import numpy
import os
import logging
import re
import itertools
import matplotlib.pyplot as plt

from pymatgen.io.abinitio.workflows import Workflow
from pymatgen.io.abinitio.tasks import TaskManager, ScfTask, NscfTask, RelaxTask
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
}

acell_hm = {
    'N': {
            'alpha': 5.491,
            'beta' : 5.641,
            'gamma': 5.571,
         },
    'P': {
            'alpha': 6.525,
            'beta' : 6.809,
            'gamma': 7.248,
         },
    'Si': {
            'alpha': 6.274,
            'beta' : 6.590,
            'gamma': 6.749,
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
    ngkpt = [6,6,6],
)

electrons = dict(
    ecut   = Energy(30.,'Ha'),
    occopt = 3,
    nsppol = 2,
    tsmear = Energy(1.e-3,'eV'),
    nband = 48,
    tolvrs = 0,
)

def valence(Z):
    return sum([int(re.findall(r'[\d]+',shell)[-1]) for shell in re.sub('<[^<]+?>','',Z.electronic_structure).split('.')[1:]])

def spin(structure):
    spins = []
    for element in structure:
        if element.specie.symbol == 'Li':
            spins.append([0.,0.,1.])
        elif element.specie.symbol == 'Mn':
            spins.append([0.,0.,7.])
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

Z = 'Si'
phase = 'beta'
flow = AbinitFlow(manager = manager, workdir = workdir)
work = Workflow()
for acell in [5.3, 5.4, 5.5, 5.6, 5.7, 5.8]:
    structure = HalfHeusler(['Li','Mn',Z], phase, acell)
    structure.make_supercell([[-1,1,1],[1,-1,1],[1,1,-1]])
    structure.remove_sites([0])
    inp = get_input(structure)
    inp.set_variables(
        **dict(
            ionmov = 7,
            ntime = 80,
            optcell = 0,
            tolmxf = 5e-05,
            toldff = 5e-06,
            amu = '*12',
        )
    )
    work.register(inp, manager = manager, task_class=RelaxTask)
flow.register_work(work)
flow = flow.allocate()

#flow = AbinitFlow(manager = manager, workdir = workdir)
#tasks = 9*[{}]
#for i,(Z,phase) in enumerate(itertools.product(['N','P','Si'],['alpha','beta','gamma'])):
#    structure_opt = HalfHeusler(['Li','Mn',Z], phase, acell_opt[Z][phase])
#    structure_hm  = HalfHeusler(['Li','Mn',Z], phase, acell_hm[Z][phase] )
#    
#    inp_opt = get_input(structure_opt)
#    inp_hm  = get_input(structure_hm)
#    
#    work_Z = Workflow()
#    scf_opt_task = work_Z.register(inp_opt, manager = manager, task_class=ScfTask)
#    scf_hm_task  = work_Z.register(inp_hm , manager = manager, task_class=ScfTask)
#    flow.register_work(work_Z)
#    tasks[i] = {
#        "opt": {
#            "input" : inp_opt,
#            "task": scf_opt_task,
#        },
#        "hm" : {
#            "input" : inp_hm,
#            "task": scf_hm_task,
#        }
#    }
#
#for i,(Z,phase) in enumerate(itertools.product(['N','P','Si'],['alpha','beta','gamma'])):
#    work_Nscf = Workflow()
#    work_Nscf.register(dos_input(tasks[i]["opt"]["input"]),   deps = {tasks[i]["opt"]["task"]: "DEN"}, manager = manager, task_class=NscfTask)
#    work_Nscf.register(bands_input(tasks[i]["opt"]["input"]), deps = {tasks[i]["opt"]["task"]: "DEN"}, manager = manager, task_class=NscfTask)
#    work_Nscf.register(dos_input(tasks[i]["hm"]["input"]),   deps = {tasks[i]["hm"]["task"]: "DEN"}, manager = manager, task_class=NscfTask)
#    work_Nscf.register(bands_input(tasks[i]["hm"]["input"]), deps = {tasks[i]["hm"]["task"]: "DEN"}, manager = manager, task_class=NscfTask)
#    flow.register_work(work_Nscf)
#        
#flow = flow.allocate()
