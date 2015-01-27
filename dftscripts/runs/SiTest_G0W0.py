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
            os.path.join(pseudo_dir, "Si.GGA_PBE-JTH-paw.xml")], 
            ndtset = 4)
    equation = ['Si', 'Si']
    structure = ZincBlende(equation, a = 5.431)
    inp.set_structure(structure)
    inp.set_kmesh(
        shiftk = [
            [.5,.5,.5],
            [.5,.0,.0],
            [.0,.5,.0],
            [.0,.0,.5],
        ],
        ngkpt = [4,4,4],
    )
    inp.set_variables(
        ecut   = "300 eV",
        pawecutdg = "600 eV",
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
    inp[1].set_kmesh(shiftk = [[.5,.5,.5],[.5,.0,.0],[.0,.5,.0],[.0,.0,.5]],ngkpt = [4,4,4])
    inp[1].set_variables(tolvrs = 1e-12)
    # Dataset 2 (WFK)
    inp[2].set_kmesh(shiftk = [[.5,.5,.5],[.5,.0,.0],[.0,.5,.0],[.0,.0,.5]],ngkpt = [2,2,2])
    inp[2].set_variables(tolwfr=1e-12, istwfk="*1", nband = 100)
    # Dataset 3 (SCR)
    inp[3].set_kmesh(shiftk = [[.5,.5,.5],[.5,.0,.0],[.0,.5,.0],[.0,.0,.5]],ngkpt = [2,2,2])
    inp[3].set_variables(
        optdriver=3, istwfk="*1",
        ecuteps = "300 eV", nband = 100
    )
    # Dataset 4 (SIG)
    inp[4].set_kmesh(shiftk = [[.5,.5,.5],[.5,.0,.0],[.0,.5,.0],[.0,.0,.5]],ngkpt = [2,2,2])
    inp[4].set_variables(
        optdriver=4, istwfk="*1",
        ecutsigx = "300 eV", nband = 100
    )
    
    return inp.split_datasets()
        
def gw_flow(workdir):
    inps = make_input();
    manager = abilab.TaskManager.from_user_config()
    flow = abilab.g0w0_flow(workdir, manager, inps[0], inps[1], inps[2], [inps[3]])
    return flow
