#!/usr/bin/env python
from __future__ import division, print_function

import os
import logging
from abipy import abilab
from myscripts.structure import Cu2Sb

logging.basicConfig()

pseudo_dir = os.path.expanduser("~/.abinit/paw")

def make_input():
    inp = abilab.AbiInput(
        pseudos = [
            os.path.join(pseudo_dir, "Li.GGA_PBE-JTH-paw.xml"),
            os.path.join(pseudo_dir, "Mn.GGA_PBE-JTH-paw.xml"),
            os.path.join(pseudo_dir, "As.GGA_PBE-JTH-paw.xml")], 
            ndtset = 4)
    equation = ['Mn','Li','Li','Li','As','As']
    structure = Cu2Sb(equation, 4.119, 6.5813382, z1 = 0.35428, z2 = 0.23058)
    inp.set_structure(structure)
    inp.set_kmesh(
        shiftk = [
            [.5,.5,.5],
        ],
        ngkpt = [15,15,11],
    )
    inp.set_variables(
        ecut   = "1000 eV",
        occopt = 3,
        nsppol = 2,
        tsmear = "0.001 eV",
        nband = 24,
    )
    spins = []
    for element in structure:
        if element.specie.symbol == 'Li':
            spins.append([0.,0.,0.])
        elif element.specie.symbol == 'Mn':
            spins.append([0.,0.,7.])
        else:
            spins.append([0.,0.,0.])
    inp.set_variables(spinat = spins)
    
    # Dataset 1 (SCF)
    inp[1].set_kmesh(shiftk = [[.5,.5,.5]],ngkpt = [15,15,11])
    inp[1].set_variables(tolvrs = 1e-8)
    # Dataset 2 (WFK)
    inp[2].set_kmesh(shiftk = [[.0,.0,.0]],ngkpt = [2, 2, 2])
    inp[2].set_variables(tolwfr=1e-12)
    # Dataset 3 (SCR)
    inp[3].set_kmesh(shiftk = [[.0,.0,.0]],ngkpt = [2, 2, 2])
    inp[3].set_variables(ecuteps = "500 eV", nband = 100)
    # Dataset 4 (SIG)
    inp[4].set_kmesh(shiftk = [[.0,.0,.0]],ngkpt = [2, 2, 2])
    inp[4].set_variables(ecutsigx = "500 eV", nband = 100)
    
    return inp.split_datasets()
        
def gw_flow(workdir):
    inps = make_input();
    manager = abilab.TaskManager.from_user_config()
    flow = abilab.g0w0_flow(workdir, manager, inps[0], inps[1], inps[2], [inps[3]])
    return flow

flow = gw_flow('/tmp/flows')

#def dos_input(inp):
#    inp = inp.deepcopy()
#    inp.set_variables(**dict(
#        iscf = -3,
#        prtdos = 2,
#        shiftk = [.0,.0,.0],
#        nshiftk = 1,
#    ))
#    return inp
#
#def bands_input(inp):
#    inp = inp.deepcopy()
#    hs = HighSymmKpath(inp.structure)
#    kpts = hs.get_kpoints()[0]
#    inp.set_variables(**dict(
#        iscf = -2,
#        kptopt = 0,
#        nkpt = len(kpts),
#        kpt = numpy.array(kpts),
#    ))
#    return inp
#
#manager = TaskManager.from_user_config()
#flow = AbinitFlow(manager = manager, workdir=workdir)
#relax_work = Workflow()
#
#equation = ['Mn','Li','Li','Li','As','As']
#
#structure = Cu2Sb(equation, Length(4.119,'ang'), Length(4.119*1.5978,'ang'), z1 = 0.35428, z2 = 0.23058)
#inp = get_input(structure)
#relax_work.register(relax_input(inp), manager = manager, task_class=RelaxTask)
#
#flow.register_work(relax_work)
#
##try:
##    structure = relax_task.read_final_structure()
##except:
##    pass
##else:
##    inp = get_input(structure)
##    scf_work = Workflow()
##    scf_task = scf_work.register(dos_input(inp), manager = manager, task_class=ScfTask)
##    nscf_work = Workflow()
##    dos_task = nscf_work.register(dos_input(inp), deps = {scf_task: ["DEN","WFK"]}, manager = manager, task_class=NscfTask)
##    bands_task = nscf_work.register(bands_input(inp), deps = {scf_task: ["DEN","WFK"]}, manager = manager, task_class=NscfTask)    
##    flow.register_work(scf_work)
##    flow.register_work(nscf_work)
#
#flow = flow.allocate()