#!/usr/bin/env python

from ase.io.vasp import read_vasp

structure = read_vasp('POSCAR')

print '      ',
for atom1 in structure:
	print '%s      ' % atom1.symbol,
print
for (i,atom1) in enumerate(structure):
	print '%s ' % atom1.symbol,
	for (j,atom2) in enumerate(structure):
		if (i >= j):
			print '        ',
		else:
			print '%8.3f' % structure.get_distance(i,j,True),
	print
