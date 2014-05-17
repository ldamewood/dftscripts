#!/usr/bin/env python
from __future__ import division, print_function

import os
import logging
import pickle

from pymatgen.core.units import Length, Energy

from pymatgen.io.abinitio.tasks import TaskManager

from abipy.htc.input import AbiInput

from myscripts.flows import PhononFlow
from myscripts.pseudos import get_psp
from myscripts.structure import HalfHeusler

scratchdir = '/p/lscratchd/damewood'
basename = os.path.dirname(os.path.abspath(__file__)).split(os.path.sep)[-1]
workdir = os.path.join(scratchdir,basename)
logging.basicConfig()

structure = HalfHeusler(['Li', 'Mn', 'P'], 'beta', Length(5.717,'ang'))

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
    ecut   = Energy(1.,'Ha'),
    occopt = 3,
    nsppol = 2,
    tsmear = Energy(1.e-3,'eV'),
    nband = 24,
)

def get_input(structure):
    pseudos = list(set([get_psp(atom.specie.symbol) for atom in structure]))
    inp = AbiInput(pseudos)
    inp.set_variables(**structure.to_abivars())
    inp.set_variables(**ksampling)
    inp.set_variables(**electrons)
    inp.set_variables(
        spinat    = [[0., 0., 0.,],
                    [0., 0., +5.],
                    [0., 0., -3.]],
        istwfk    = '*1',       # do not use time reversal symmetry
        nstep     = 2000,
        tolvrs    = 1.e-2,
    )
    return inp

def load_flow(inp, runpath):
    pickle_file = os.path.join(runpath,'__AbinitFlow__.pickle')
    try:
        flow = pickle.Unpickler(open(pickle_file)).load()
    except IOError as e:
        manager = TaskManager.from_user_config()
        flow = PhononFlow(runpath, manager, inp, ngqpt=[2,2,2], do_nscf = False)
        flow.ph_tolvrs = 1.e-2
        flow.build_tasks()
        flow = flow.allocate()
    return flow

flow = load_flow(get_input(structure), workdir)
