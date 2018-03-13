# -*- coding: utf-8 -*-

#    read_energy_levels.py is part of Frigus.

#    Frigus: software to compure the energy exchange in a multi-level system
#    Copyright (C) 2016-2018 Mher V. Kazandjian and Carla Maria Coppola

#    Frigus is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, version 3 of the License.    
#
#    Frigus is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with Frigus.  If not, see <http://www.gnu.org/licenses/>.

"""
module that implements energy levels readers
"""
from __future__ import print_function
import numpy
from numpy import loadtxt

from astropy import units as u
from astropy import constants

from frigus.species import EnergyLevelsMolecule, EnergyLevelsOnDegreeOfFreedom
from frigus.utils import linear_2d_index


def read_levels(fname):
    """
    Parse the  data of

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
    energies *= 1.98630e-23  # .. todo:: do this conversion with astropy

    nv_max, nj_max = v.max(), j.max()

    enh2 = numpy.zeros((nv_max+1, nj_max+1), 'f8')

    # having v=0 as 0 for the energies, the bottom of the well
    # is at -2179.0 cm-1 + convention into Joule
    # bottom = 2179.0*1.98630e-23

    # filling the matrix enh2 with the values read from the text file
    for i in range(len(v)):
        vi, ji = v[i], j[i]
        enh2[vi, ji] = energies[i]

        # depth = enh2[0,0] + bottom
        # enh2[vi,ji]  =  enh2[vi,ji] - depth

    return enh2


def read_levels_lique(fname):
    """
    Read the energy levels from the file "H2Xvjlevels_francois_mod.cs" that is
    recovered from ABC. The file has (should have) the following format.

      ----------------------------------------------------------------------
      Channel   a       v       j       k       Energy (eV)
      ----------------------------------------------------------------------
       92       2       0       0       0         0.27034
       93       2       0       1       0         0.28504

    :param fname: The paths to the ascii file containing the levels
     information.
    :return: numpy array array of N rows, the v and j levels and their
     corresponding energy. The elements of the array are stored in increasing
     order in the energy.

    .. code-block:: python

        levels = read_levels_lique('/path/to/H2Xvjlevels_francois_mod.cs')

        print('{:3}{:3}{:10}'.format('v', 'j', 'E(eV)'))
        for level in levels:
           print('{:3}{:3}{:10}'.format(level['v'], level['j'], level['E']))
    """
    data_read = numpy.loadtxt(fname, skiprows=3)

    # get the sorting indices list from the last column (energy)
    inds = numpy.argsort(data_read[:, -1])

    v, j, energies = data_read[inds, 2], data_read[inds, 3], data_read[inds, 5]

    energy_levels = EnergyLevelsMolecule(inds.size, energy_unit=u.eV)

    energy_levels.data['v'] = v
    energy_levels.data['j'] = j
    energy_levels.data['g'] = 2*j + 1
    energy_levels.data['E'] = energies * u.eV
    energy_levels.data['label'] = linear_2d_index(v, j)

    return energy_levels


def read_levels_wrathmall_and_flower(fname):
    """
    read the energy levels from the file "H2Xvjlevels_flower.cs". The file has
    (should have) the following format.

      lev no.   v           J       Energy (K)

      1       0       0       0
      2       0       2       509.85
      3       0       4       1681.678

    :param fname: The paths to the ascii file containing the levels
     information.
    :return: numpy array array of N rows, the v and j levels and their
     corresponding energy. The elements of the array are stored in increasing
     order in the energy.

    .. code-block:: python

        levels = read_levels_wrathmall_and_flower('/path/to/H2Xvjlevels_flower.cs')

        print('{:3}{:3}{:10}'.format('v', 'j', 'E(eV)'))
        for level in levels:
           print('{:3}{:3}{:10}'.format(level['v'], level['j'], level['E']))
    """
    data_read = numpy.loadtxt(fname, skiprows=2)

    # get the sorting indices list from the last column (energy)
    inds = numpy.argsort(data_read[:, -1])

    v, j, energies = data_read[inds, 1], data_read[inds, 2], data_read[inds, 3]

    energy_levels = EnergyLevelsMolecule(inds.size, energy_unit=u.K)

    energy_levels.data['v'] = v
    energy_levels.data['j'] = j
    energy_levels.data['g'] = 2*j + 1
    energy_levels.data['E'] = (energies*u.K).to(
        u.eV,
        equivalencies=u.temperature_energy()
    )
    energy_levels.data['label'] = linear_2d_index(v, j)

    return energy_levels


def read_levels_lipovka(fname):
    """
    read the energy levels from the file "flower_roueff_data.dat" that is
    recovered from the website:
                                                                  v(?)      j (?)     j(?)    E (K? or cm-1?)
      0.000     0.000     0.000     0.000     0.000     0.000     0.000     0.000     0.000     0.000
      0.000     0.000     0.000     0.000     0.000     1.000     0.000     1.000     1.000   128.400

    :param fname: The paths to the ascii file containing the levels
     information.
    :return: EnergyLevelsMolecule object
    """
    v, j, energies = [], [], []
    with open(fname) as fobj:
        linesold = fobj.readlines()
        lines = []
        for line_no, line in enumerate(linesold):
            if line.isspace() is False and line_no >= 2 and line_no <= 12:
                line_split = line.split()
                v.append(0)
                j.append(int(float(line_split[-2])))
                energies.append(float(line_split[-1]))

    v = numpy.array(v, 'i')
    j = numpy.array(j, 'i')
    energies = numpy.array(energies, 'f8')

    # assume the data is sorted in increasing energy values
    energy_levels = EnergyLevelsMolecule(energies.size, energy_unit=u.K)

    energy_levels.data['v'] = v
    energy_levels.data['j'] = j
    energy_levels.data['g'] = 2*j + 1
    energy_levels.data['E'] = (energies * u.K * constants.k_B).to(u.eV)
    energy_levels.data['label'] = linear_2d_index(v, j)

    return energy_levels


def read_levels_two_levels_test_1(fname):
    """
    Read the energy levels of the two_levels_test_1 system.

    :param fname: The paths to the ascii file containing the levels
     information.
    :return: a EnergyLevelsMolecule object holding the level data
    """

    data_read = numpy.loadtxt(fname, skiprows=3)

    # get the sorting indices list from the last column (energy)
    inds = numpy.argsort(data_read[:, -1])

    j, g, energies = data_read[inds, :].T

    energy_levels = EnergyLevelsOnDegreeOfFreedom(
        n_levels=inds.size,
        energy_unit=u.eV
    )

    energy_levels.data['j'] = j
    energy_levels.data['g'] = g
    energy_levels.data['E'] = energies*u.eV
    energy_levels.data['label'] = j

    return energy_levels
