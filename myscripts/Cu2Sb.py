#!/usr/bin/env python
from __future__ import division, print_function

import os
import logging
import numpy

from pymatgen.core.units import Length, Energy

from pymatgen.io.abinitio.tasks import TaskManager, RelaxTask, NscfTask, ScfTask
from pymatgen.io.abinitio.flows import AbinitFlow
from pymatgen.io.abinitio.workflows import Workflow
from pymatgen.symmetry.bandstructure import HighSymmKpath

from abipy.htc.input import AbiInput

from myscripts.pseudos import get_psp
from myscripts.structure import Cu2Sb

scratchdir = '/p/lscratchd/damewood'
basename = os.path.dirname(os.path.abspath(__file__)).split(os.path.sep)[-1]
workdir = os.path.join(scratchdir,basename)
logging.basicConfig()

structure = Cu2Sb(['Mn','Li','Li','Li','Sb','Sb'], Length(4.326,'ang'), Length(1.58 * 4.326,'ang'))

ksampling = dict(
    kptopt = 1,
    nshiftk = 1,
    shiftk = [
        [.5,.5,.5],
    ],
    ngkpt = [15,15,11],
)

electrons = dict(
    ecut   = Energy(1.,'Ha'),
    occopt = 3,
    nsppol = 2,
    tsmear = Energy(1.e-3,'eV'),
    nband = 24,
)

def spin(structure):
    spins = []
    for element in structure:
        if element.specie.symbol == 'Li':
            spins.append([0.,0.,1.])
        elif element.specie.symbol == 'Mn':
            spins.append([0.,0.,7.])
        else:
            spins.append([0.,0.,-3])
    return dict(spinat = spins)

def get_input(structure, method = 'paw'):
    pseudos = get_psp(structure,method=method)
    inp = AbiInput(pseudos)
    inp.set_variables(**structure.to_abivars())
    inp.set_variables(**ksampling)
    inp.set_variables(**electrons)
    inp.set_variables(**spin(structure))
    return inp

def relax_input(inp, optcell = 1):
    if optcell == 0:
        tolmxf = 1.e-10
        toldfe = 0
    else:
        toldfe = 1.e-10
        tolmxf = 0
    inp = inp.deepcopy()
    inp.set_variables(**dict(
        ntime = 100,
        ionmov = 3,
        optcell = optcell,
        ecutsm = Energy(0.5, 'Ha'),
        dilatmx = 1.2,
        tolmxf = tolmxf,
        toldfe = toldfe,
    ))
    return inp
    
def dos_input(inp):
    inp = inp.deepcopy()
    inp.set_variables(**dict(
        iscf = -3,
        prtdos = 2,
        shiftk = [.0,.0,.0],
        nshiftk = 1,
    ))
    return inp

def bands_input(inp):
    inp = inp.deepcopy()
    hs = HighSymmKpath(inp.structure)
    kpts = hs.get_kpoints()[0]
    inp.set_variables(**dict(
        iscf = -2,
        kptopt = 0,
        nkpt = len(kpts),
        kpt = numpy.array(kpts),
    ))
    return inp

def update_structure(fromtask, totasks):
    structure = fromtask.read_final_structure()
    for task in totasks:
        task.strategy.abinit_input.set_structure(structure)
        task.build()

manager = TaskManager.from_user_config()
flow = AbinitFlow(manager = manager, workdir=workdir)
inp = get_input(structure)

relax_work1 = Workflow()
relax_task1 = relax_work1.register(relax_input(inp, optcell = 0), manager = manager, task_class=RelaxTask)

relax_work2 = Workflow()
relax_task2 = relax_work2.register(relax_input(inp, optcell = 2), manager = manager, task_class=RelaxTask)

relax_task1.on_ok = lambda: update_structure(relax_task1, [relax_task2])

scf_work = Workflow()
scf_task = scf_work.register(dos_input(inp), manager = manager, task_class=ScfTask)

nscf_work = Workflow()
dos_task = nscf_work.register(dos_input(inp), deps = {scf_task: ["DEN","WFK"]}, manager = manager, task_class=NscfTask)
bands_task = nscf_work.register(bands_input(inp), deps = {scf_task: ["DEN","WFK"]}, manager = manager, task_class=NscfTask)

relax_task2.on_ok = lambda: update_structure(relax_task2, [scf_task, dos_task, bands_task])

flow.register_work(relax_work1)
flow.register_work(relax_work2)
flow.register_work(scf_work)
flow.register_work(nscf_work)
flow = flow.allocate()