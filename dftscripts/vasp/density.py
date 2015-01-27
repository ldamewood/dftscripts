#!/usr/bin/env python
# plot_density.py

import argparse
import numpy
from pymatgen.io.vaspio.vasp_output import Chgcar
from pymatgen.electronic_structure.core import Spin

from .util.interpolate import interpolate_plane, plane_to_cs

def _build_parser():
    # Setup parser
    parser = argparse.ArgumentParser(prog='plot_chgcar',description='2d projection of CHG')

    # Add options
    parser.add_argument('filename',type=str,nargs='?',default='CHG',
        help = 'The CHG file (default: CHG)')
    parser.add_argument('-s','--spin-channel',choices=['total','polarization','up','down'],default='total',
        help = 'For spin polarized calculations, the component of the density (default: total)')

    # Set default commands
    # Coordinate system options
    group_coord = parser.add_argument_group('coordinate system')
    group_coord.add_argument('-ca','--center-atom',nargs='?', type=int,
        help = 'Atom at the center (default: None)')
    group_coord.add_argument('-cc','--center-cartesian',nargs=3, type=float,
        help = 'Cartesian coordinate for the center of the plot (default: None)')
    group_coord.add_argument('-cr','--center-direct',nargs=3, type=float, default=[0.5,0.5,0.5],
        help = 'Direct (fractional) coordinate for the center (default: 0.5 0.5 0.5)')
    group_coord.add_argument('-xa','--x-axis',nargs=3, type=float, default=[1.,0.,0.],
        help = 'Direction of the x-axis (default: 1 0 0)')
    group_coord.add_argument('-ya','--y-axis',nargs=3, type=float, default=[0.,0.,1.],
        help = 'Direction of the y-axis (default: 0 0 1)')

    # Interpolation options
    group_interpolate = parser.add_argument_group('interpolation')
    group_interpolate.add_argument('-xrep','--x-repeat',nargs='?', type=int, default=1, 
        help = 'Approximate number of repetitions of the cell in the x-direction (default: 1)')
    group_interpolate.add_argument('-yrep','--y-repeat',nargs='?', type=int, default=1,
        help = 'Approximate number of repetitions of the cell in the y-direction (default: 1)')
    group_interpolate.add_argument('-xres','--x-res',nargs='?', type=int, default=100,
        help = 'Resolution of the x-axis (default: 100)')
    group_interpolate.add_argument('-yres','--y-res',nargs='?', type=int, default=100,
        help = 'Resolution of the y-axis (default: 100)')

    # Output parameters
    group_output = parser.add_argument_group('output')
    group_output.add_argument('-fmt','--format',choices=['xyz','matrix'],default='xyz',
        help = 'Output format: xyz = columns of x,y,z data, matrix = raw density data (default: xyz)')
    group_output.add_argument('-gz','--gzip', action='store_true',
        help='Enables gzip compression of the output (default: disabled)')
    return parser

def main():
    # Set default variables
    options = dict(
        spin_channel = 'total',
        center_cartesian = [0.,0.,0.],
        x_axis = [1.,0.,0.],
        y_axis = [0.,0.,1.],
        x_repeat = 1, # approximate number of unit cells to interpolate
        y_repeat = 1, # approximate number of unit cells to interpolate
        x_res = 500, # points to interpolate
        y_res = 500, # points to interpolate
        format = 'xyz',
    )
    options.update(**vars(_build_parser().parse_args()))

    chg = Chgcar.from_file(options['filename'])
    structure = chg.structure
    rprim = chg.structure.lattice.matrix
    # Decide where to obtain center from
    # Cartisian option set
    center = options['center_direct']
    
    if options['center_cartesian'] is not None:
        # Convert from cartesian to reduced coordinates
        center = numpy.linalg.solve(rprim,options['center_cartesian']).tolist()
    
    # Atomic option is set
    if options['center_atom'] is not None:
        # Convert from atom number to reduced coordinates
        center = structure[options['center_atom']].frac_coords.tolist()
    
    plane = numpy.array([options['x_axis'],options['y_axis']],dtype=float)
    repeat = [options['x_repeat'],options['y_repeat']]
    res = [options['x_res'],options['y_res']]
    if options['spin_channel'] == 'total':
        data = chg.data['total']
    if options['spin_channel'] == 'up':
        data = chg.spin_data[Spin.up]
    if options['spin_channel'] == 'down':
        data = chg.spin_data[Spin.down]
    if options['spin_channel'] == 'polarization':
        data = chg.spin_data[Spin.up] - chg.spin_data[Spin.down]
    
    data /= structure.volume
    
    cs = plane_to_cs(plane)
    posinplane = numpy.dot(cs,numpy.dot(structure.frac_coords - center, rprim).T).T
    print('Atom positions in plane:')
    for i,pos in enumerate(posinplane):
        print('  %s: %f %f (z = %f)' % (structure[i].specie, pos[0], pos[1], pos[2]))
    
    # Interpolate
    x,y,z = interpolate_plane(data,rprim,plane=plane,center=center,dim=repeat,res=res)
    print('Cell density:\n  minimum: %f\n  maximum: %f' % (data.min(),data.max()))
    print('In-plane density:\n  minimum: %f\n  maximum: %f' % (z.min(),z.max()))
    print('Note: PAW calculations may contain negative values for the pseudo-density in the core regions.')
    
    suffix = u'.gz' if options['gzip'] == True else u''
    
    if options['format'] == 'xyz':
        out = numpy.array([x.flatten(),y.flatten(),z.flatten()]).T
        with open('chg_%s_in_plane.csv%s' % (options['spin_channel'],suffix), 'w+') as f:
            numpy.savetxt(f, out, delimiter = u',', header = u' X, Y, Z')
            
    if options['format'] == 'matrix':
        with open('chg_%s_in_plane.mat%s' % (options['spin_channel'],suffix), 'w+') as f:
            numpy.savetxt(f, z, delimiter = u' ', header = u' Z matrix')