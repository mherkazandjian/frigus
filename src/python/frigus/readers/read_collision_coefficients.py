# -*- coding: utf-8 -*-

#    read_collision_coefficients.py is part of Frigus.

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
Read the collisional rates (collision coefficients)
"""
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

import numpy
from numpy import (loadtxt, int32, zeros, unique, void,
                   ascontiguousarray, dtype, hstack)

from astropy import units as u


def unique_level_pairs(vj):
    """
    From a list of levels find the list of unique levels and return them.

    :param iterable vj: the list of v levels where each item vj[x] is a level.
    The shape of vj should be (2,n) where n is the number of levels.
    :return: (ndarray) The unique levels. The shape of this array is
    (2,n_unique) where n_unique is the number of unique transitions.

    .. code-block: python

        vj = array([[0, 0, 0, 7, 4, 8, 4, 6, 9, 8, 9, 6, 2, 0, 5, 5, 9, 8, 1],
                    [4, 4, 3, 5, 3, 6, 3, 8, 2, 5, 7, 8, 8, 6, 6, 9, 1, 0, 9]])
        vj_unique = unique_level_pairs(vj)
    """
    assert vj.shape[0] == 2

    a = vj.T
    new_dtype = dtype((void, a.dtype.itemsize * 2))
    b = ascontiguousarray(a).view(new_dtype)
    _, idx = unique(b, return_index=True)
    unique_a = a[idx]
    return unique_a.T


def read_collision_coefficients_lique_and_wrathmall(fname):
    """
    Parse the collisional data by François.

    These are the coefficient rates K_ij where i > j (so these fill the lower
    triangular K matrix).

    Supplementary Information for manuscript:
    "Revisited study of the ro-vibrational excitation of H$_2$ by H: Towards a
    revision of the cooling of astrophysical media" by François LIQUE

    The table contains the H2-H collisional rate coefficients

    v, v' : initial and final vibrational state
    j, j' : initial and final rotational state

    v  j  v' j'    k(cm3 s-1) (T) , T= 100 to 5000K by steps of 100K

    0  1  0  0     0.3098E-21  0.1521E-18  0.1791E-16  0.3164E-15.....
    0  2  0  0     0.6561E-13  0.7861E-13  0.9070E-13  0.1266E-12.....

    :param string fname: The path to the ascii data.
    :return: a tuple of 3 elements.
      The first element is a 5D array that holds all the rate coefficients.
      K[v, j, v', j', T_index] ( in m3/s)

      The second element is a 1D temperature array (T) corresponding to the
      rate coefficients in the 5D array. This array has the same size as the
      last dimension of K (i.e K[0,0,0,0,:].size = T.size

      The third element is a tuple of 4 elements:

         - The first element is a 2D array of shape (2, n_transitions) (ini).
           The columns of this array (v, j = ini) are the initial v and j of
           the transitions for a certain T section. i.e. the number of non-zero
           elements in K for a certain temperature is equal to the number of
           elements in the v or j columns.
           v.size = j.size = where(K[..., 0] > 0)[0].size

         - The second element is the same of the first element but for the
           final transitions (v', j' = fin)

         - The third element is a 2D array of shape (2, n_unique_transitions)
           that are the unique levels involved in all the transitions.

         - The last element is an array of shape (T.size, n_transitions)
           which are the collisional coefficient rates with non-zero values
           for each value of temperature in the T array.
    """
    def find_temperature_array():
        """
        Search the header of the data file and return the range of
        temperatures used. The temperature range is assumed to be on the
        10th line of the header
        """
        _t_vals = None
        with open(fname) as fobj:
            for line_num, line in enumerate(fobj):
                if line_num == 8:
                    _t_vals = loadtxt(StringIO(line), delimiter=',') * u.Kelvin
                    break

        assert _t_vals is not None
        return _t_vals

    t_values = find_temperature_array()

    # read the data from the original ascii file
    data_read = loadtxt(fname, unpack=True, skiprows=10)
    (v, j, vp, jp), cr = int32(data_read[0:4]), data_read[4:]
    ini = zeros((2, v.size), 'i')
    fin = zeros((2, v.size), 'i')

    # declare the array where the data will be stored
    nv, nj, nvp, njp = v.max()+1, j.max()+1, vp.max()+1, jp.max()+1
    nv_max, nj_max = int(max(nv, nvp)), int(max(nj, njp))
    data = zeros((t_values.size, nv_max, nj_max, nv_max, nj_max), 'f8')

    # copy the read data into the container array
    for i, cri in enumerate(cr.T):
        data[:, v[i], j[i], vp[i], jp[i]] = cri
        ini[:, i] = v[i], j[i]
        fin[:, i] = vp[i], jp[i]

    # find the unique levels from from the transitions
    unique_levels = unique_level_pairs(hstack((unique_level_pairs(ini),
                                               unique_level_pairs(fin))))

    # set the units of the data to be returned
    data_with_units = data * (u.cm**3 / u.second)
    cr_with_units = cr * (u.cm**3 / u.second)

    # convert the units to m^3/s
    data_with_units = data_with_units.to(u.m**3 / u.second)
    cr_with_units = cr_with_units.to(u.m**3 / u.second)

    return data_with_units, t_values, (ini, fin, unique_levels, cr_with_units)


def read_collision_coefficients_lipovka(fname):
    """
    Parse the collisional data used by lipovka.

    These are the coefficient rates K_ij where i > j (so these fill the lower
    triangular K matrix).

    The table contains the HD-H collisional rate coefficients

    v, v' : initial and final vibrational state
    j, j' : initial and final rotational state

    the data are stored in blocks for each temperature. The colomns are the
    initial ro-vibrational (v,j) level and the rows are the final ones (v', j')
    http://massey.dur.ac.uk/drf/HD_H/

    :param string fname: The path to the ascii data.
    :return: a tuple of 3 elements.
      The first element is a 5D array that holds all the rate coefficients.
      K[T_index, v, j, v', j'] ( in m3/s)

      The second element is a 1D temperature array (T) corresponding to the
      rate coefficients in the 5D array. This array has the same size as the
      last dimension of K (i.e K[:,0,0,0,0].size = T.size

      The third element is a tuple of 4 elements:

         - The first element is a 2D array of shape (2, n_transitions) (ini).
           The columns of this array (v, j = ini) are the initial v and j of
           the transitions for a certain T section. i.e. the number of non-zero
           elements in K for a certain temperature is equal to the number of
           elements in the v or j columns. It is assumed that 'ini' are the 
           same for each temperature value. 
           v.size = j.size = where(K[..., 0] > 0)[0].size

         - The second element is the same of the first element but for the
           final transitions (v', j' = fin)

         - The third element is a 2D array of shape (2, n_unique_transitions)
           that are the unique levels involved in all the transitions.

         - The last element is an array of shape (T.size, n_transitions)
           which are the collisional coefficient rates with non-zero values
           for each value of temperature in the T array. Each row corresponds
           to the ini, fin collision rates of the first and second element
           in the tuple mentioned earlier.
    """
    class Reader(object):
        """
        Parse the  data by Flower and Roueff contained in the file
        flower_roueff_data.dat downloaded from http://massey.dur.ac.uk/drf/HD_H

        A block of data (for a certain temperature) is defined as everything
        between:

        100
        (0,0) (0,1) (0,2) (0,3) (0,4) (0,5) (0,6) (0,7) (0,8) (0, 9)
        1.0D-09 1.5D-11 5.4D-13 4.0D-14 7.8D-15 6.1D-16 2.4D-17 8.7D-18 2.2D-18 2.3D-19
        ...
        ...

        .. code-block:: python

            reader = read_cr_lipovka.Reader('path_to_the_file')

            # print all the collision rates for the transition (v,j) -> (vp,jp)
            # for
            # all the temperatures
            print reader.data[v, j, vp, jp, :]
        """
        def __init__(self, fname, tiny=1e-70):
            """Constructor"""

            self.fname = fname
            self.data = None
            """
            the collision rates for all the transitions for all the tempratures
            """

            self.ini = None   #: the initial v,j of all the transitions
            self.fin = None   #: the final v,j of all the transitions
            self.tkin = None
            """the temperatures at which the collisional data are given"""

            self.read_data()

        def read_data(self):
            """Read all the data to temporary storage in self.data"""
            with open(self.fname) as fobj:
                linesold = fobj.readlines()
                lines = []
                for line_no, line in enumerate(linesold):
                    if line.isspace() is False and line_no >= 12:
                        lines.append(line)
                raw_data = ''.join(lines)

            cr, ini, fin, tkin = self.parse_data(raw_data)
            self.data = cr
            self.ini = ini
            self.fin = fin
            self.tkin = tkin

        def parse_data(self, raw_data):
            """
            Parses the read data into blocks, one block for each temperature
            """
            # split the raw data into blocks
            blocks = raw_data.split('T =')[1:]

            tkin = zeros(len(blocks), 'f8')
            cr_tkin = []
            ini_tkin = []
            fin_tkin = []

            for ib, block in enumerate(blocks):
                ini, fin, t, cr = self.parse_block(block)
                tkin[ib] = t
                cr_tkin.append(cr)
                ini_tkin.append(ini)
                fin_tkin.append(fin)

            return cr_tkin, ini_tkin, fin_tkin, tkin

        def parse_block(self, block):
            """
            Parse a block of data and return the temperature, the levels and
            the collision rates
            """
            lines = iter(filter(lambda x: len(x) > 0, block.split('\n')))

            # get the temperature
            t_kin = next(lines)
            assert 'K' in t_kin
            tkin = t_kin.replace('K', '')
            # print T, tkin

            # get the header (levels)
            levels = next(lines).replace(
                ' ', ''
            ).replace(
                ')(', ':'
            )[1:-1].split(':')

            v, j = [], []
            for level in levels:
                v.append(int32(level.split(',')[0]))
                j.append(int32(level.split(',')[1]))
            v, j = numpy.array(v), numpy.array(j)

            ini = numpy.repeat(numpy.vstack((v, j)), len(v), axis=1)
            fin = numpy.tile([v, j], len(v))

            cr = zeros((v.size, j.size, v.size, j.size), 'f8')
            for final_state_ind, line in enumerate(lines):
                vp, jp = v[final_state_ind], j[final_state_ind]
                if len(line.strip()) > 0:
                    data = numpy.float64(line.replace('D', 'E').strip().split())
                    for initial_state_ind, datum in enumerate(data):
                        vi, ji = v[initial_state_ind], j[initial_state_ind]
                        cr[vi, ji, vp, jp] = datum

            return ini, fin, tkin, cr

    reader = Reader(fname)

    t_values = reader.tkin

    # read the data from the original ascii file
    # data_read = loadtxt(fname, unpack=True, skiprows=10)
    # (v, j, vp, jp), cr = int32(data_read[0:4]), data_read[4:]

    ini = reader.ini[0]
    fin = reader.fin[0]

    v, j = ini[0, :], ini[1, :]
    vp, jp = fin[0, :], fin[1, :]

    n_transitions_per_t_value = v.size

    # declare the array where the data will be stored in a tensor
    nv, nj, nvp, njp = v.max()+1, j.max()+1, vp.max()+1, jp.max()+1
    nv_max, nj_max = max(nv, nvp), max(nj, njp)
    data = zeros((t_values.size, nv_max, nj_max, nv_max, nj_max), 'f8')

    # declare the array where the data will be stored in a matrix
    cr = numpy.zeros((t_values.size, n_transitions_per_t_value), 'f8')

    # copy the read data into the container array
    for i, cri in enumerate(reader.data):
        data[i, v, j, vp, jp] = cri[v, j, vp, jp]
        cr[i, :] = cri[v, j, vp, jp]

    # find the unique levels from from the transitions
    unique_levels = unique_level_pairs(hstack((unique_level_pairs(ini),
                                               unique_level_pairs(fin))))

    # set the units of the data to be returned
    data_with_units = data * (u.cm**3 / u.second)
    cr_with_units = cr * (u.cm**3 / u.second)
    t_values = t_values * u.K

    # convert the units to m^3/s
    data_with_units = data_with_units.to(u.m**3 / u.second)
    cr_with_units = cr_with_units.to(u.m**3 / u.second)

    return data_with_units, t_values, (ini, fin, unique_levels, cr_with_units)
