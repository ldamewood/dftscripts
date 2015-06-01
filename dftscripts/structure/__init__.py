import numpy

from pymatgen.core.units import ArrayWithUnit


def structure_to_abivars(self):
    """Returns a dictionary with the abinit variables."""
    types_of_specie = self.types_of_specie
    natom = self.num_sites

    znucl_type = [specie.number for specie in types_of_specie]

    # znucl_atoms = self.atomic_numbers

    typat = numpy.zeros(natom, numpy.int)
    for (atm_idx, site) in enumerate(self):
        typat[atm_idx] = types_of_specie.index(site.specie) + 1

    rprim = ArrayWithUnit(self.lattice.matrix, "ang").to("bohr")
    xred = numpy.reshape([site.frac_coords for site in self], (-1, 3))

    # Set small values to zero. This usually happens when the CIF file
    # does not give structure parameters with enough digits.
    # rprim = np.where(np.abs(rprim) > 1e-8, rprim, 0.0)
    # xred = np.where(np.abs(xred) > 1e-8, xred, 0.0)

    d = dict(
        natom=natom,
        ntypat=len(types_of_specie),
        typat=typat,
        xred=xred,
        znucl=znucl_type)

    d.update(dict(
        acell=3 * [1.0],
        rprim=rprim))

    # d.update(dict(
    #     acell=3 * [1.0],
    #     angdeg))

    return d
