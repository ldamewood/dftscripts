import numpy

from matplotlib import pyplot as plt

from pymatgen.io.abinitio.netcdf import NetcdfReader
from pymatgen.io.abinitio.abiinspect import GroundStateScfCycle

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

def check_hm(flow):
    tasks = []
    i = 0
    for fl in flow[:9]:
        for task in fl:
            tasks.append(task)
    for i,task in enumerate(flow[9]):
        dos = task.outdir.has_abiext('DOS')
        if dos:
            print(tasks[i]._name)
            ncdf = NetcdfReader(task.outdir.has_abiext('GSR'))
            acell = 2*ncdf.read_variable('primitive_vectors')[0][1]
            cycle = GroundStateScfCycle.from_file(tasks[i].output_file.path)
            moment = cycle.last_iteration['magn']
            etotal = cycle.last_etotal
            print('  acell:' + str(acell*0.529177249))
            print('  momen:' + str(moment))
            print('  etotl:' + str(etotal))
            is_hm(dos)
            

def is_hm(dosfile):
    x,dos,efermi = dos_from_file(dosfile)
    x = x - efermi
    izero = numpy.abs(x).argmin()
    
    istart = izero
    iend = izero
    while dos[0,istart] < 1e-5:
        istart = istart - 1
    while dos[0,iend] < 1e-5:
        iend = iend + 1
    if iend > istart:
        gapup = x[iend]-x[istart]
        print('    gap dn: '+str(gapup*27.3))
    
    istart = izero
    iend = izero
    while dos[1,istart] < 1e-5:
        istart = istart - 1
    while dos[1,iend] < 1e-5:
        iend = iend + 1
    if iend > istart:
        gapdn = x[iend]-x[istart]
        print('    gap up: '+str(gapdn*27.3))