import numpy
import os

from pymatgen.core.units import Energy, EnergyArray

from pymatgen.io.abinitio.netcdf import NetcdfReader

from abipy.electrons.edos import ElectronDOS
from abipy.electrons.ebands import ElectronBands, ElectronsReader

__all__ = ['dos_from_file']

def dos_from_file(filename, limits = [-0.1,0.1]):
    lines = open(filename).readlines()[:15]
    dimline = lines[3].split(',')
    try:
        nsppol = int(dimline[0].split()[-1])
    except:
        nsppol = 1
    #nkpt = int(dimline[1].split()[-1])
    #nband = int(dimline[2].split()[-1])
    efermi = Energy(float(lines[7].split()[-1]),'Ha')
    #nene = int(lines[10].split()[2])
    #elo = float(lines[11].split()[2])
    #ehi = float(lines[11].split()[4])
    dos = numpy.loadtxt(filename)
    n = dos.shape[0]
    dos = dos.reshape(nsppol,n/nsppol,3)
    x = EnergyArray(dos[0,:,0],'Ha')
    return x,dos[:,:,1],efermi

def bands_from_file(filename, limit = [-0.1,0.1]):
    
    if not os.path.exists(filename):
        raise Exception()
    
    reader = ElectronsReader(filename)
    bands = ElectronBands(
	structure=reader.read_structure(),
	kpoints=reader.read_kpoints(),
        eigens=reader.read_eigenvalues(),
        fermie=reader.read_fermie(),
        occfacts=reader.read_occupations(),
        nelect=reader.read_nelect(),
        nband_sk=reader.read_nband_sk(),
        smearing=reader.read_smearing(),
    )
    bands._fix_fermie()
    return bands
    

#def plotdos(dosfilename):
#    dos = loaddos(dosfilename)
#    dosup,dosdn = dos[:,numpy.logical_and(dos[1,:,0] > -0.1,dos[1,:,0] < 0.1),:2]
#    plt.plot(dosup[:,0],+dosup[:,1])
#    plt.plot(dosdn[:,0],-dosdn[:,1])
#    plt.show()