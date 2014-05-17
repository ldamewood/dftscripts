#!/usr/bin/env python

import numpy
import os
import logging
import re
import matplotlib.pyplot as plt

from pymatgen.io.abinitio.workflows import Workflow
from pymatgen.io.abinitio.tasks import TaskManager, RelaxTask, ScfTask, NscfTask
from pymatgen.io.abinitio.flows import AbinitFlow

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
    ngkpt = [12,12,12],
)

electrons = dict(
    ecut   = Energy(40.,'Ha'),
    occopt = 3,
    nsppol = 2,
    tsmear = Energy(1.e-3,'eV'),
    nband = 16,
    tolvrs = 1.e-10,
)

electrons_ZB = dict(
    ecut   = Energy(40.,'Ha'),
    occopt = 3,
    nsppol = 2,
    tsmear = Energy(1.e-3,'eV'),
    nband = 16,
    tolvrs = 1.e-10,
    spinat = [
        [0.0,0.0,+6.0],
        [0.0,0.0,-3.0],
    ],
    prtdos = 2,
    nstep = 2000,
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

def relax_input(structure):
    pseudos = [get_psp(atom.specie.symbol) for atom in structure]
    inp = AbiInput(pseudos)
    inp.set_variables(**structure.to_abivars())
    inp.set_variables(**ksampling)
    inp.set_variables(**electrons)
    inp.set_variables(**dict(
        ntime = 100,
        ionmov = 2,
        optcell = 1,
        ecutsm = Energy(0.5, 'Ha'),
        dilatmx = 1.2,
        tolmxf = 1.e-7,
    ))
    inp.set_variables(**spin(structure))
    return inp

def bands_input(structure):
    pseudos = [get_paw(atom.specie.symbol) for atom in structure]
    inp = AbiInput(pseudos)
    inp.set_variables(**structure.to_abivars())
    inp.set_variables(**ksampling)
    inp.set_variables(**electrons)
    inp.set_variables(**dict(
        prtdos = 2,
        kptopt = 1,
        nshiftk = 4,
        shiftk = [
            [.0,.0,.0],
            [.0,.5,.5],
            [.5,.0,.5],
            [.5,.5,.0],
        ],
        ngkpt = [12,12,12],
    ))
    inp.set_variables(**spin(structure))
    return inp

def fix_input(structure):
    pseudos = [get_paw(atom.specie.symbol) for atom in structure]
    inp = AbiInput(pseudos)
    inp.set_variables(**structure.to_abivars())
    inp.set_variables(**ksampling)
    inp.set_variables(**electrons)
    inp.set_variables(**spin(structure))
    return inp

def loaddos(dosfilename):
    lines = open(dosfilename).readlines()[:15]
    dimline = lines[3].split(',')
    try:
        nsppol = int(dimline[0].split()[-1])
    except:
        nsppol = 1
    #nkpt = int(dimline[1].split()[-1])
    #nband = int(dimline[2].split()[-1])
    efermi = float(lines[7].split()[-1])
    #nene = int(lines[10].split()[2])
    #elo = float(lines[11].split()[2])
    #ehi = float(lines[11].split()[4])
    dos = numpy.loadtxt(dosfilename)
    dos[:,0] = dos[:,0] - efermi
    n = dos.shape[0]
    dos = dos.reshape(nsppol,n/nsppol,3)
    return dos

def plotdos(dosfilename):
    dos = loaddos(dosfilename)
    dosup,dosdn = dos[:,numpy.logical_and(dos[1,:,0] > -0.1,dos[1,:,0] < 0.1),:2]
    plt.plot(dosup[:,0],+dosup[:,1])
    plt.plot(dosdn[:,0],-dosdn[:,1])
    plt.show()

flow = AbinitFlow(manager = manager, workdir = workdir)
for Z in ['N','P','Si']:
    for phase in ['alpha','beta','gamma']:
        work_Z = Workflow()
        
        structure = HalfHeusler(['Li','Mn',Z], phase, acell_opt[Z][phase])
        relax_task = work_Z.register(relax_input(structure), manager = manager, task_class=RelaxTask)
	nscf_task = work_Z.register(bands_input(structure), deps = {relax_task: "DEN"}, manager = manager, task_class=NscfTask)
        
        structure = HalfHeusler(['Li','Mn',Z], phase, acell_hm[Z][phase])
        hm_task = work_Z.register(fix_input(structure), manager = manager, task_class=ScfTask)
        nscf_task = work_Z.register(bands_input(structure), deps = {hm_task: "DEN"}, manager = manager, task_class=NscfTask)
        
        flow.register_work(work_Z)
flow = flow.allocate()
