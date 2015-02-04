#!/usr/bin/env python
from __future__ import print_function, division

import sys
import numpy
import itertools
import scipy.io.netcdf as netcdf

def levicivita():
    """Generates the nonzero components and value of the Levi-Civita symbol. """
    yield +1.,(0,1,2) # levicivita = +1 for (x,y,z), etc.
    yield +1.,(1,2,0)
    yield +1.,(2,0,1)
    yield -1.,(1,0,2)
    yield -1.,(0,2,1)
    yield -1.,(2,1,0)

def pauli(s1,s2,component):
    """Gives the complex components of the pauli matrices."""
    # Convert spins to 0 or 1
    s1 = bool(s1)
    s2 = bool(s2)
    if component == 0:
        return complex(0.) if s1 == s2 else complex(1.)
    elif component == 1:
        return complex(0.) if s1 == s2 else complex(0., 2. * s1 - 1.)
    elif component == 2:
        return complex(1. - 2. * s1) if s1 == s2 else complex(0.)
    else:
        # Invalid cartesian component: just return zero.
        return complex(0.)

def gVectors(shape):
    """
    Construct the G-vectors that correspond to the grid shape.
    
    Ex: `gGrid([2,2,3])` will return 12 G-vectors with the FFT ordering.
    """
    assert numpy.size(shape) == 3
    # FFT components in each direction
    # 0, 1, 2 ... N/2, N/2-1, ... -1
    x = numpy.fft.fftfreq(shape[0],1./shape[0])
    y = numpy.fft.fftfreq(shape[1],1./shape[1])
    z = numpy.fft.fftfreq(shape[2],1./shape[2])
    # Create a grid of points
    x,y,z = numpy.meshgrid(x,y,z)
    # Flatten into a list of vectors: i.e. vec[:,n] is n-th Fourier component
    grid = numpy.array([x.flatten(), y.flatten(), z.flatten()]).T
    # 0 0 0
    # 0 0 1
    # ...
    # -1 -1 -1
    # Return as integer array
    return numpy.array(grid, dtype=int)

class WaveFunction:
    def __init__(self, ncfile):
        self._netcdf = netcdf.netcdf_file(ncfile, 'r')
        #self.nband = self._netcdf.dimensions['max_number_of_states']
        ## Number of k-points
        #self.nkpt = self._netcdf.dimensions['number_of_kpoints']
        ## Number of spin components
        #self.nspin = self._netcdf.dimensions['number_of_spins']
        ## Number of G-vectors
        #self.nG = [0,0,0]
        #self.nG[0] = self._netcdf.dimensions['number_of_grid_points_vector1']
        #self.nG[1] = self._netcdf.dimensions['number_of_grid_points_vector2']
        #self.nG[2] = self._netcdf.dimensions['number_of_grid_points_vector3']
        ## Real or complex wave functions
        #self._rorc = self._netcdf.dimensions['real_or_complex_wavefunctions']
        self.prim = self._netcdf.variables['primitive_vectors'][:,:]
        self.rprim = 2*numpy.pi*numpy.linalg.inv(self.prim).T
        self.kpts = ncfile.variables['reduced_coordinates_of_kpoints'][:,:]

    def get_wavefunction(self, iband, ikpt, ispin, ispinor = 0):
        # Read dimensions
        [nspin,nkpt,nband,nspinor,nG1,nG2,nG3,rorc] = \
            netcdf.variables['real_space_wavefunctions'].shape
        # Valid arguments?
        assert iband < nband
        assert ikpt < nkpt
        assert ispin < nspin
        assert ispinor < nspinor
        # Read wave function components
        wfk = self._netcdf.variables['real_space_wavefunctions'][ispin,ikpt,iband,ispinor,:,:,:,:]
        # Convert to complex
        if rorc == 2: wfk.dtype = '>c16'
        # Reshape to G-vector dimensions
        wfk = numpy.reshape(wfk,[nG1,nG2,nG3])
        # Fourier transform
        return numpy.fft.ifftn(wfk[:,:,:])

def matrix_element(wf, iband1, ikpt1, ispin1, iband2, ikpt2, ispin2, potin):
    """
    Calculate the matrix element between two wave functions using the provided potential.
    
    input:
    i*1, i*2     : wave function quantum numbers for band, kpt and spin.
    pot          : F.T. of the potential
    
    global:
    rprim (primitive vectors), rkpt (all direct k-points)
    
    returns: matrix element for the spin orbit interaction.
    """
    
    if ikpt1 != ikpt2: return 0.
    
    # Get the wave functions on G-grid
    wfk1 = wf.get_wavefunction(iband1, ikpt1, ispin1)
    wfk2 = wf.get_wavefunction(iband2, ikpt2, ispin2)

    # Potential on G-grid (make a copy since we may resize it)
    potential = potin.copy()

    # The k-point for the second (ket) wavefunction.
    kpt = wf.kpts[ikpt2]
    
    # Find maximum G grid in each direction.
    maxGrid = numpy.max([wfk1.shape, wfk2.shape, potential.shape], axis = 0)

    # rotate G vectors in cartesian coordinates
    gvmax = numpy.dot(wf.rprim, gVectors(maxGrid).T).T
    
    # Pad the wave functions and potential with zeros.
    wfk1.resize(maxGrid) # G'
    wfk2.resize(maxGrid) # G''
    potential.resize(maxGrid) # G'' - G'
    
    # Calculate column vectors for wfk1, Ej and Pk
    wf1 = wfk1.reshape([len(gvmax),1])
    Pk = (kpt + gvmax) * wfk2.reshape([len(gvmax),1])
    Ej = gvmax * potential.reshape([len(gvmax),1])

    me = []
    # Loop over nonzero terms of the Levi-Civita symbol (6 terms) to sum up the matrix element
    for  levi,(i,j,k) in levicivita():
        # i,j,k are the cartesian components
        
        # Check if this term is automatically zero
        psigma = pauli(ispin1, ispin2, i)
        if abs(psigma) < 1.E-12: continue

        # Convolution: view Ej and Pk on the grid, then use FFT convolution
        EjPk = numpy.fft.ifftn(numpy.fft.fftn(Ej[:,j].reshape(maxGrid)) * \
                               numpy.fft.fftn(Pk[:,k].reshape(maxGrid)))

        # Dot product: view convolution as column vector and dot with wave function.
        me.append(-levi * psigma * numpy.dot(wf1.T.conj(), EjPk.reshape(len(gvmax)))[0])
    
    # Add up all the nonzero ijk terms
    return sum(me)