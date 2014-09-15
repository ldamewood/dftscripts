#!/usr/bin/env python
from __future__ import division, print_function

import os
import logging
import numpy

from pymatgen.core.units import Length, Energy

from pymatgen.io.abinitio.tasks import TaskManager, RelaxTask, NscfTask, ScfTask
from pymatgen.io.abinitio.flows import AbinitFlow
from pymatgen.io.abinitio.workflows import Workflow
from pymatgen.io.abinitio.abiobjects import RelaxationMethod, Electrons, KSampling
from pymatgen.symmetry.bandstructure import HighSymmKpath

from abipy.htc.input import AbiInput

from myscripts.pseudos import get_psp
from myscripts.structure import Cu2Sb

scratchdir = '/p/lscratchd/damewood'
workdir = os.path.join(scratchdir,'Cu2Sb/abinit')
logging.basicConfig()

relax = dict(
    ionmov = 3,
    ntime = 80,
    optcell = 0,
    tolmxf = 1e-06,
    toldff = 1e-07,
)

ksampling = dict(
    kptopt = 1,
    nshiftk = 1,
    shiftk = [
        [.5,.5,.5],
    ],
    ngkpt = [15,15,11],
)

electrons = dict(
    ecut   = Energy(1000.,'eV'),
    occopt = 3,
    nsppol = 2,
    tsmear = Energy(1.e-3,'eV'),
    nband = 24,
)

def spin(structure):
    spins = []
    for element in structure:
        if element.specie.symbol == 'Li':
            spins.append([0.,0.,0.])
        elif element.specie.symbol == 'Mn':
            spins.append([0.,0.,7.])
        else:
            spins.append([0.,0.,0.])
    return dict(spinat = spins)

def get_input(structure, method = 'paw'):
    pseudos = get_psp(structure,method=method)
    inp = AbiInput(pseudos)
    inp.set_variables(**structure.to_abivars())
    inp.set_variables(**ksampling)
    inp.set_variables(**electrons)
    inp.set_variables(**spin(structure))
    return inp

def relax_input(inp):
    inp = inp.deepcopy()
    inp.set_variables(**relax)
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

manager = TaskManager.from_user_config()
flow = AbinitFlow(manager = manager, workdir=workdir)
relax_work = Workflow()

equation = ['Mn','Li','Li','Li','As','As']

structure = Cu2Sb(equation, Length(4.119,'ang'), Length(4.119*1.5978,'ang'), z1 = 0.35428, z2 = 0.23058)
inp = get_input(structure)
relax_work.register(relax_input(inp), manager = manager, task_class=RelaxTask)

flow.register_work(relax_work)

#try:
#    structure = relax_task.read_final_structure()
#except:
#    pass
#else:
#    inp = get_input(structure)
#    scf_work = Workflow()
#    scf_task = scf_work.register(dos_input(inp), manager = manager, task_class=ScfTask)
#    nscf_work = Workflow()
#    dos_task = nscf_work.register(dos_input(inp), deps = {scf_task: ["DEN","WFK"]}, manager = manager, task_class=NscfTask)
#    bands_task = nscf_work.register(bands_input(inp), deps = {scf_task: ["DEN","WFK"]}, manager = manager, task_class=NscfTask)    
#    flow.register_work(scf_work)
#    flow.register_work(nscf_work)

flow = flow.allocate()