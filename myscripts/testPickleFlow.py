import os

from pymatgen.io.abinitio.flows import Flow
from pymatgen.io.abinitio.tasks import TaskManager

workdir = '/p/'
name = 'test'

flow1 = Flow(manager = TaskManager.from_user_config(), workdir = os.path.join(workdir, name))
flow1.build_and_pickle_dump()

flow2 = Flow.pickle_load(os.path.join(workdir, name))