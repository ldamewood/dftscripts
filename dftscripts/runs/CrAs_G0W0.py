#!/usr/bin/env python
from __future__ import division, print_function

import os
import logging
from abipy import abilab
from myscripts.structure import ZincBlende

logging.basicConfig()

pseudo_dir = os.path.expanduser("~/.abinit/paw")

def make_input():
    inp = abilab.AbiInput(
        pseudos = [
            os.path.join(pseudo_dir, "Cr.GGA_PBE-JTH-paw.xml"),
            os.path.join(pseudo_dir, "As.GGA_PBE-JTH-paw.xml")], 
            ndtset = 4)
    equation = ['Cr', 'As']
    structure = ZincBlende(equation, a = 5.66)
    inp.set_structure(structure)
    inp.set_kmesh(
        shiftk = [
            [.5,.5,.5],
        ],
        ngkpt = [15,15,11],
    )
    inp.set_variables(
        ecut   = "1600 eV",
        pawecutdg = "2000 eV",
        occopt = 3,
        nsppol = 2,
        tsmear = "0.001 eV",
        nband = 24,
    )
    spins = []
    for element in structure:
        if element.specie.symbol == 'As':
            spins.append([0.,0.,3.])
        elif element.specie.symbol == 'Cr':
            spins.append([0.,0.,-5.])
        else:
            spins.append([0.,0.,0.])
    inp.set_variables(spinat = spins)
    
    # Dataset 1 (SCF)
    inp[1].set_kmesh(shiftk = [[.5,.5,.5],[.5,.0,.0],[.0,.5,.0],[.0,.0,.5]],ngkpt = [20,20,20])
    inp[1].set_variables(tolvrs = 1e-12)
    # Dataset 2 (WFK)
    inp[2].set_kmesh(shiftk = [[.5,.5,.5],[.5,.0,.0],[.0,.5,.0],[.0,.0,.5]],ngkpt = [2,2,2])
    inp[2].set_variables(tolwfr=1e-20, istwfk="*1", nband = 1000)
    # Dataset 3 (SCR)
    inp[3].set_kmesh(shiftk = [[.5,.5,.5],[.5,.0,.0],[.0,.5,.0],[.0,.0,.5]],ngkpt = [2,2,2])
    inp[3].set_variables(
        optdriver=3, istwfk="*1",
        ecuteps = "1400 eV", ecutwfn = "1600 eV", nband = 1000
    )
    # Dataset 4 (SIG)
    inp[4].set_kmesh(shiftk = [[.5,.5,.5],[.5,.0,.0],[.0,.5,.0],[.0,.0,.5]],ngkpt = [2,2,2])
    inp[4].set_variables(
        optdriver=4, istwfk="*1",
        ecutsigx = "1400 eV", ecutwfn = "1600 eV", nband = 1000
    )
    
    return inp.split_datasets()
        
def gw_flow(workdir):
    inps = make_input();
    manager = abilab.TaskManager.from_user_config()
    flow = abilab.g0w0_flow(workdir, manager, inps[0], inps[1], inps[2], [inps[3]])
    return flow