# -*- coding: utf-8 -*-
"""

"""
from __future__ import print_function
import numpy 
from numpy import loadtxt, genfromtxt
import pylab


def read_levels(fname):
    """parse the  data of

    :param fname: The name of the ascii file containing the levels information.
     (see below for an example of the content of the file). The file can be
     found at http://www.physast.uga.edu/ugamop/energies/H2Xvjlevels
    :return: A numpy matrix of shape (vmax + 1, jmax + 1) where vmax and jmax
      are the maximum vibrational and rotational transition quantum numbers
      mentioned in the file "fname".

    to numpy arrays and provide utilities to get information from the data.
      # Rovibrational energy levels for H_2
      # v   J        -binding energy
      #           a.u.            cm-1
        0   0  0.16456626E+00  0.36118114E+05
        0   1  0.16402653E+00  0.35999657E+05
        0   2  0.16295202E+00  0.35763830E+05
        0   3  0.16135249E+00  0.35412773E+05
        0   4  0.15924216E+00  0.34949609E+05
        0   5  0.15663934E+00  0.34378356E+05
        0   6  0.15356587E+00  0.33703808E+05
                   .
                   .
                   .
                   .
       13   7  0.23025618E-03  0.50535384E+02
       14   0  0.65673931E-03  0.14413760E+03
       14   1  0.57867630E-03  0.12700475E+03
       14   2  0.42998862E-03  0.94371580E+02
       14   3  0.22670691E-03  0.49756409E+02
    """

    v, j, a, energies = loadtxt(fname, unpack=True)

    v, j = numpy.int32(v), numpy.int32(j)

    # conversion of the energies from cm-1 -> Joule
    energies *= 1.98630e-23

    nv_max, nj_max = v.max(), j.max()

    enh2 = numpy.zeros((nv_max+1, nj_max+1), 'f8')

    # having v=0 as 0 for the energies, the bottom of the well
    # is at -2179.0 cm-1 + convention into Joule
    # bottom = 2179.0*1.98630e-23

    # filling the matrix enh2 with the values read from the text file
    for i in range(len(v)):
        vi, ji = v[i], j[i]
        enh2[vi, ji] = energies[i]

#        depth = enh2[0,0] + bottom
#        enh2[vi,ji]  =  enh2[vi,ji] - depth

    return enh2


def read_levels_lique(fname):
    """
    read the energy levels from the file "H2Xvjlevels_francios_mod.cs" that is
    recovered from ABC. The file has (should have) the following format.

      ----------------------------------------------------------------------
      Channel   a       v       j       k       Energy (eV)
      ----------------------------------------------------------------------
       92       2       0       0       0         0.27034
       93       2       0       1       0         0.28504

    :param fname: The paths to the ascii file containing the levels
     information. The file should have the following format:
    :return: numpy array array of N rows, the v and j levels and their
     corresponding energy. The elements of the array are stored in increasing
     order in the energy.

    .. code-block:: python

        levels = read_levels('/path/to/H2Xvjlevels_francios_mod.cs')

        print('{:3}{:3}{:10}'.format('v', 'j', 'E(eV)'))
        for level in levels:
           print('{:3}{:3}{:10}'.format(level['v'], level['j'], level['E']))
    """
    data_read = numpy.loadtxt(fname, skiprows=3)

    # get the sorting indices list from the last column (energy)
    inds = numpy.argsort(data_read[:, -1])

    v, j, energies = data_read[inds, 2], data_read[inds, 3], data_read[inds, 5]

    data = numpy.zeros(inds.shape, dtype=[('v', 'i4'),
                                          ('j', 'i4'),
                                          ('E', 'f8')])

    data['v'] = v
    data['j'] = j
    data['E'] = energies

    return data
