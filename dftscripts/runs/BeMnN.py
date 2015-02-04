#!/usr/bin/env python
from __future__ import division, print_function

import os
import logging
from abipy import abilab
from dftscripts.structure.cubic import HalfHeusler

logging.basicConfig()

pseudo_dir = os.path.expanduser("~/.abinit/paw")

def make_input():
    inp = abilab.AbiInput(
        pseudos = [
            os.path.join(pseudo_dir, "Be.GGA_PBE-JTH.xml"),
            os.path.join(pseudo_dir, "Mn.GGA_PBE-JTH.xml"),
            os.path.join(pseudo_dir, "N.GGA_PBE-JTH.xml")], 
            ndtset = 1)
    equation = ['Be', 'Mn', 'N']
    structure = HalfHeusler(equation, 'gamma', a = 12.283218926)
    inp.set_structure(structure)
    inp.set_kmesh(
        shiftk = [
            [.5,.5,.5],
            [.5,.0,.0],
            [.0,.5,.0],
            [.0,.0,.5],
        ],
        ngkpt = [12,12,12],
    )
    inp.set_vars(
        ecut   = 40,
        pawecutdg = 60,
        occopt = 3,
        nsppol = 2,
        tsmear = "0.04 eV",
        nband = 24,
    )
    spins = []
    for element in structure:
        if element.specie.symbol == 'Mn':
            spins.append([0.,0.,3.])
        else:
            spins.append([0.,0.,0.])
    inp.set_vars(spinat = spins)
        
    return inp
    
from dftscripts.flows.phonon import auto_phonon_flow
manager = abilab.TaskManager.from_user_config()
flow = auto_phonon_flow('.', manager, make_input(), with_nscf = True)