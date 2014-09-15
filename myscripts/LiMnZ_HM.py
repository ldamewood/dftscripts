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
from pymatgen.io.abinitio.netcdf import NetcdfReader
from pymatgen.io.abinitio.abiinspect import GroundStateScfCycle

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

acell_low = {
    'N': {
            'alpha': 5.25,
            'beta' : 0,
            'gamma': 0,
         },
    'P': {
            'alpha': 6.15,
            'beta' : 6.40,
            'gamma': 0,
         },
    'Si': {
            'alpha': 0,
            'beta' : 6.10,
            'gamma': 6.20,
         },
}

acell_hm = {
    'N': {
            'alpha': 5.341,
            'beta' : 0,
            'gamma': 0,
         },
    'P': {
            'alpha': 6.375,
            'beta' : 6.659,
            'gamma': 0,
         },
    'Si': {
            'alpha': 0,
            'beta' : 6.274,
            'gamma': 6.396,
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

def dos_from_file(filename, limits = [-0.1,0.1]):
    lines = open(filename).readlines()[:15]
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
    dos = numpy.loadtxt(filename)
    n = dos.shape[0]
    dos = dos.reshape(nsppol,n/nsppol,3)
    x = dos[0,:,0]
    return x,dos[:,:,1],efermi

def plot_opt():
    work_task = [9,10,11,12,13,14,15,16,17]
    work_name = [
        "LiMnN alpha",  "LiMnN beta",  "LiMnN gamma",
        "LiMnP alpha",  "LiMnP beta",  "LiMnP gamma",
        "LiMnSi alpha", "LiMnSi beta", "LiMnSi gamma",
    ]
    
    for i in range(9):
        filename = 'work_%s/task_2/outdata/out_DOS' % str(work_task[i])
        print(work_name[i])
        print(is_hm(filename))
        #x,dos,efermi = dos_from_file(filename)
        #x = x - efermi
        #idx = numpy.logical_and(x >= -0.1,x <= 0.1)
        #plt.plot(x[idx],dos[0,idx],x[idx],-dos[1,idx])
        #plt.ylim([-150,150])
        #plt.show()

def check_hm(flow, n = 9):
    tasks = []
    i = 0
    for fl in flow[:n]:
        for task in fl:
            tasks.append(task)
    for i,task in enumerate(flow[n]):
        dos = task.outdir.has_abiext('DOS')
        if dos:
            ncdf = NetcdfReader(task.outdir.has_abiext('GSR'))
            acell = 2*ncdf.read_variable('primitive_vectors')[0][1]
            cycle = GroundStateScfCycle.from_file(tasks[i].output_file.path)
            moment = cycle.last_iteration['magn']
            etotal = cycle.last_etotal
            gapup,gapdn = is_hm(dos)
            name = tasks[i]._name.split('_')
            alloy = name[0]
            phase = name[1]
            print '%s\t%s\t%f\t%f\t%f\t%f\t%f' % (alloy,phase,acell,etotal,moment,gapup,gapdn)
            

def is_hm(dosfile):
    x,dos,efermi = dos_from_file(dosfile)
    x = x - efermi
    izero = numpy.abs(x).argmin()
    
    gapup = 0.
    gapdn = 0.
    istart = izero
    iend = izero
    while dos[0,istart] < 1e-5:
        istart = istart - 1
    while dos[0,iend] < 1e-5:
        iend = iend + 1
    if iend > istart:
        gapup = x[iend]-x[istart]
    
    istart = izero
    iend = izero
    while dos[1,istart] < 1e-5:
        istart = istart - 1
    while dos[1,iend] < 1e-5:
        iend = iend + 1
    if iend > istart:
        gapdn = x[iend]-x[istart]
    return (gapup,gapdn)

flow = AbinitFlow(manager = manager, workdir = workdir)
tasks = []
inpts = []
for Z,phase in itertools.product(['N','P','Si'],['alpha','beta','gamma']):
    work = Workflow()
    if acell_low[Z][phase] < 0.05:
        continue
    for acell in numpy.arange(acell_low[Z][phase],acell_hm[Z][phase]+0.05,0.05):
        structure = HalfHeusler(['Li','Mn',Z], phase, acell)
        inp = get_input(structure)
        task = work.register(inp, manager = manager, task_class=ScfTask)
        task._name = 'LiMn%s_%s_%f' % (Z,phase,acell)
        print(task._name)
        tasks.append(task)
        inpts.append(inp)
    flow.register_work(work)
work = Workflow()
for task,inp in zip(tasks,inpts):
    work.register(dos_input(inp), deps={task:"DEN"},
                  manager=manager, task_class=NscfTask)
flow.register_work(work)
flow = flow.allocate()
