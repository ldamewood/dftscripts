#!/usr/bin/env python

import numpy
import os
import logging

from pymatgen.io.abinitio.workflows import Workflow
from pymatgen.io.abinitio.tasks import TaskManager, RelaxTask
from pymatgen.io.abinitio.flows import AbinitFlow

from abipy.core.structure import Structure, Lattice
from abipy.htc.input import AbiInput
from abipy.core.constants import Length, Energy

from myscripts.pseudos import get_paw

scratchdir = '/p/lscratchd/damewood'
basename = os.path.dirname(os.path.abspath(__file__)).split(os.path.sep)[-1]
workdir = os.path.join(scratchdir,basename)
logging.basicConfig()
manager = TaskManager.from_user_config()

acell = {
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

ksampling = dict(
    kptopt = 1,
    nshiftk = 1,
    shiftk = [
        [.5,.5,.5],
        [.5,.0,.0],
        [.0,.5,.0],
        [.0,.0,.5],
    ],
    ngkpt = [12,12,12],
)

electrons = dict(
    ecut   = Energy(40.,'Ha'),
    occopt = 3,
    nsppol = 2,
    tsmear = Energy(1.e-3,'eV'),
    nband = 24,
)

def LiMnZ_structure(Z, phase, acell):
    wk = {
        '4a' : numpy.array([.00,.00,.00]),
        '4b' : numpy.array([.50,.50,.50]),
        '4c' : numpy.array([.25,.25,.25]),
        '4d' : numpy.array([.75,.75,.75]),
    }
    positions = {
        'alpha' : [wk['4c'],wk['4b'],wk['4a']],
        'beta'  : [wk['4b'],wk['4a'],wk['4c']],
        'gamma' : [wk['4a'],wk['4c'],wk['4b']],
        'zincb' : [wk['4a'],wk['4c']],
    }
    if hasattr(acell,'to') and callable(getattr(acell,'to')):
        acell = acell.to('bohr')
    lattice = Lattice(float(acell) * numpy.array([[.5,.5,.0],[.5,.0,.5],[.0,.5,.5]]))
    if phase in ['alpha','beta','gamma']:
        elements = ['Li','Mn',str(Z)]
    if phase in 'zincb':
        elements = ['Mn',str(Z)]
    return Structure(lattice, elements, positions[phase])

def relax_input(structure):
    pseudos = [get_paw(atom.specie.symbol) for atom in structure]
    inp = AbiInput(pseudos = pseudos)
    inp.set_variables(**structure.to_abivars())
    inp.set_variables(**ksampling)
    inp.set_variables(**electrons)
    inp.set_variables(**dict(
        ionmov = 0,
        optcell = 1,
    ))
    return inp

flow = AbinitFlow(manager = manager, workdir = workdir)
for Z in ['N','P','Si']:
    work_Z = Workflow()
    for phase in ['alpha','beta','gamma']:
        structure = LiMnZ_structure(Z, phase, acell[Z][phase])
        work_Z.register(relax_input(structure), manager = manager, task_class=RelaxTask)
    flow.register_work(work_Z)
flow = flow.allocate()