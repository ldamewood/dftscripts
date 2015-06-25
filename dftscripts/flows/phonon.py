#

from pymatgen.io.abinitio.flows import phonon_flow

from ..tasks.faketasks import get_qpts

__all__ = ['auto_phonon_flow']


def auto_phonon_flow(workdir, manager, inp, ngqpt=(4, 4, 4), with_nscf=False):
    qpts, mem = get_qpts(inp.deepcopy(), ngqpt, workdir=workdir)
    ph_inputs = []
    for qpt in qpts:
        ph_input = inp.deepcopy()
        ph_input.set_vars(dict(qpt=qpt))
        ph_inputs.append(ph_input)
    return phonon_flow(workdir, manager, inp, ph_inputs, with_nscf=with_nscf)
