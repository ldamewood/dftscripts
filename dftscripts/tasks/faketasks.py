import os

import sys

from pymatgen.io.abinitio.works import Work
from pymatgen.io.abinitio.tasks import TaskManager
from pymatgen.io.abinitio.netcdf import NetcdfReader
from pymatgen.io.abinitio.abiinspect import yaml_read_irred_perts

__all__ = [
    'get_memory',
    'get_all_kpoints',
    'get_qpts',
    'get_irred_perts',
]        

def _setup_faketask(name, manager = None, workdir = "."):
    tmp_dir = os.path.join(workdir, "__" + str(name) + "_run__")
    if manager is None:
        manager = TaskManager.from_user_config()
    manager = manager.to_shell_manager()
    return Work(workdir=tmp_dir, manager=manager)

def get_memory(inp, manager = None, workdir="."):

    vars = dict(
        prtvol = -1,
    )

    w = _setup_faketask('mem', manager = manager, workdir = workdir)
    fake_input = inp.deepcopy()
    fake_input.set_variables(**vars)
    fake_task = w.register(fake_input)
    w.allocate()
    w.start(wait=True)

    mem = 0
    for line in fake_task.log_file.readlines():
        if "P This job should need less than" in line:
            mem = float(line.split()[7])
            break

    w.rmtree()
    return mem

def get_all_kpoints(inp, manager = None, workdir = "."):

    vars = dict(
        prtvol = -1,
        kptopt = 3,
    )

    w = _setup_faketask('kpt', manager = manager, workdir = workdir)
    fake_input = inp.deepcopy()
    fake_input.set_variables(**vars)
    fake_task = w.register(fake_input)
    w.allocate()
    w.start(wait=True)

    nc = NetcdfReader(fake_task.opath_from_ext('OUT'))
    kpts = nc.read_variable('kpt')[:].reshape((-1,3))
    mem = 0
    for line in fake_task.log_file.readlines():
        if "P This job should need less than" in line:
            mem = float(line.split()[7])
            break

    w.rmtree()
    return kpts, mem

def get_qpts(inp, ngqpt, manager = None, workdir = "."):

    vars = dict(
        ngkpt = ngqpt,
        shiftk = [0.,0.,0.],
        nshiftk = 1,
        prtvol = -2,
    )

    w = _setup_faketask('qpt', manager = manager, workdir = workdir)
    fake_input = inp.deepcopy()
    fake_input.set_variables(**vars)
    fake_task = w.register(fake_input)
    w.allocate()
    w.start()

    nc = NetcdfReader(fake_task.opath_from_ext('OUT'))
    qpts = nc.read_variable('kpt')[:].reshape((-1,3))
    mem = 0
    for line in fake_task.log_file.readlines():
        if "P This job should need less than" in line:
            mem = float(line.split()[7])
            break

    w.rmtree()
    return qpts, mem

def get_irred_perts(qpt, inp, manager = None, workdir = "."):

    vars = dict(
        rfphon=1,
        nqpt=1,
        qpt=qpt,
        paral_rf=-1,
        rfatpol=[1, len(inp.structure)],
        rfdir=[1, 1, 1],
    )

    w = _setup_faketask('ph', manager = manager, workdir = workdir)
    fake_input = inp.deepcopy()
    fake_input.set_variables(**vars)
    fake_task = w.register(fake_input)
    w.allocate()
    w.start(wait=True)

    irred_perts = yaml_read_irred_perts(fake_task.log_file.path)
    mem = 0
    for line in fake_task.log_file.readlines():
        if "P This job should need less than" in line:
            mem = float(line.split()[7])
            break

    w.rmtree()
    return irred_perts, mem
