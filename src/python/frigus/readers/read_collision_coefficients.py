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
from astropy.constants import k_B


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
      K[T_index, v, j, v', j'] ( in m3/s)

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

         - The third element is a 2D array of shape (2, n_unique_levels)
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
    unique_levels = unique_level_pairs(
        hstack((unique_level_pairs(ini),
                unique_level_pairs(fin)))
    )

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


def read_collision_coefficients_esposito_h2_he(fname):
    """
    Parse the collisional data by Fabrizio Esposito.

    These are the coefficient rates K_ij and K_ji (so these fill the lower and
    upper triangular K matrix). However, only the downwards coefficients are
    accurate; the upwards must be found with detailed balance.

    As in the format used in the convention tvjwk (file HeH2_tvjwk.res)

    The table contains the H2-He collisional rate coefficients
    index v j v' j'

    index : counter for the transition, starting from 0
    v, v' : initial and final vibrational state
    j, j' : initial and final rotational state

    T [K]      k(cm3 s-1)
    T from 100 K to 1000K by steps of 100K, from 1000 K to 10000 K by step of
    500 K, from 10000 K to 20000 K by step of 1000 K

    Each transition is a block like:

    0 0 0 0 0
    100 2.8454E-009
    200 3.6106E-009
    300 3.8642E-009
    400 3.9323E-009
    ........
    ........
    18000 6.0465E-009
    19000 6.1032E-009
    20000 6.1530E-009


    followed by 2 empty lines.

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

    class Reader(object):
        """
        Parse the  data by Fabrizio Esposito

        """
        def __init__(self, fname, tiny=1e-70):
            """Constructor"""

            self.fname = fname
            self.data = None
            self.cr = None
            """
            the collision rates for all the transitions for all the tempratures
            """
            self.ini = None   #: the initial v,j of all the transitions
            self.fin = None   #: the final v,j of all the transitions
            self.tkin = None
            """the temperatures at which the collisional data are given"""

            self.read_data()


        # def find_temperature_array(self):
        #     """
        #     Parse the first block in the data file and return the range of
        #     temperatures used.
        #     """
        #     _t_vals = None
        #     with open(self.fname) as fobj:
        #         _t_vals = []
        #         for line_num, line in enumerate(fobj):
        #             if 1 <= line_num <= 38:
        #                 temp = float(line.split(sep=' ')[0])
        #                 _t_vals.append(temp)
        #
        #     assert _t_vals is not None
        #     return numpy.array(_t_vals) * u.Kelvin
        def read_data(self):
            with open(self.fname) as fobj:
                linesold = fobj.readlines()
                lines = []
                _v = []
                _j = []
                _vp = []
                _jp = []
                _ini = []
                _fin = []
                _t_vals = []
                for line_no, line in enumerate(linesold):
                    if line.isspace() is False:
                        if len(line.split()) > 2:
                            lines.append('block ' + line)
                            (
                            _index_read, _v_read, _j_read, _vp_read, _jp_read) = \
                                line.split()
                            _v.append(int(_v_read))
                            _j.append(int(_j_read))
                            _vp.append(int(_vp_read))
                            _jp.append(int(_jp_read))
                            _ini.append([int(_v_read), int(_j_read)])
                            _fin.append([int(_vp_read),int(_jp_read)])
                        else:
                            lines.append(line)
                            if 1 <= line_no <= 38:
                                temp = float(line.split(sep=' ')[0])
                                _t_vals.append(temp)

                raw_data = ''.join(lines)
                _v = numpy.array(_v)
                _j = numpy.array(_j)
                _vp = numpy.array(_vp)
                _jp = numpy.array(_jp)
                _ini = numpy.array(_ini)
                _fin = numpy.array(_fin)
                _t_vals = numpy.array(_t_vals)

            # declare the array where the data will be stored in a tensor
            nv, nj, nvp, njp = _v.max() + 1, _j.max() + 1, _vp.max() + 1, _jp.max() + 1
            nv_max, nj_max = max(nv, nvp), max(nj, njp)
            data = zeros((_t_vals.size, nv_max, nj_max, nv_max, nj_max), 'f8')
            cr = zeros((_t_vals.size, len(_v)), 'f8')

            self.ini = _ini
            self.fin = _fin
            self.tkin = _t_vals

            data, cr = self.parse_data(cr, data, raw_data)

            self.data = data

        def parse_block(self, block):
            """
            Parse a block of data and return the transition, the levels and
            the collision rates.
            It returns an array of floats corresponding to the reaction rates
            for each temperature of the temperature array
            """
            lines = list(filter(lambda x: len(x) > 0, block.split('\n')))

            cr_block = zeros(self.tkin.size, 'f8')

            for ntemp, line in enumerate(lines[1:]):
                cr_block[ntemp] = line.split()[1]

            return cr_block

        def parse_data(self, cr, data, raw_data):
            """
            Parses the read data into blocks, one block for each transition
            """
            # split the raw data into blocks
            blocks = raw_data.split('block')[1:]
            blocks = numpy.array(blocks)

            for ib, block in enumerate(blocks):
                cr_block = self.parse_block(block)
                # assign the reaction rates for all the temperatures
                # for the ib-th transition
                cr[:, ib] = cr_block[:]

                # assign the data for the multidimensional transition tensor
                v_ini, j_ini = self.ini[ib]
                v_fin, j_fin = self.fin[ib]
                data[:, v_ini, j_ini, v_fin, j_fin] = cr[:, ib]

            self.data = data
            self.cr = cr

            return data, cr

    reader = Reader(fname)

    ini = reader.ini.T
    fin = reader.fin.T

    # find the unique levels from from the transitions
    unique_levels = unique_level_pairs(
        hstack(
            (unique_level_pairs(ini),
            unique_level_pairs(fin))
        )
    )

    # set the units of the data to be returned
    data_with_units = reader.data * (u.cm**3 / u.second)
    cr_with_units = reader.cr * (u.cm**3 / u.second)
    t_values = reader.tkin * u.K

    # convert the units to m^3/s
    data_with_units = data_with_units.to(u.m**3 / u.second)
    cr_with_units = cr_with_units.to(u.m**3 / u.second)

    return data_with_units, t_values, (ini, fin, unique_levels, cr_with_units)

def compute_collision_coefficients_gerlich_h2_hp(fname1, fname2):
    """
    Parse the coefficients k0 and Delta_E0 provided by Gerlich in

    https://aip.scitation.org/doi/pdf/10.1063/1.457980

    to compute the reaction rates for the process:
    H+ + H2(j) -> H+ + H2(j')
    The data are fully provided for the transitions up to j = 7.

    These are the coefficient rates K_ij; K_ji (upwards transitions)
    must be found with detailed balance.

    The thermal reactions can be calculated as:

    K_j_j'(T) = k0_j_j' * exp(-Delta_E0_j_j'/k/T)

    In the read file there are 2 tables:
     - first table: k0 (units: 10**-10 cm3/s)
     - second table: Delta_E0 (units: meV)

    # indices: (initial,final for j, v=0 for all)
    #   j'  0     1       2      3      4      5      6     7
    # j
    # 0
    # 1
    # 2
    # 3
    # 4
    # 5
    # 6
    # 7

    The temperature range in which these fits are valid is:

    t_kin in [10, 500] K

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

    class Reader(object):
        """
        Parse the data by Gerlic
        """
        def __init__(self, fname, tiny=1e-70):
            """Constructor"""

            self.fname1 = fname1
            self.fname2 = fname2
            self.data = None
            self.cr = None
            """
            the collision rates for all the transitions for all the temperatures
            """
            self.ini = None   #: the initial v,j of all the transitions
            self.fin = None   #: the final v,j of all the transitions
            self.tkin = None
            """the temperatures at which the collisional data are given"""

            self.read_data()

        def read_data(self):
            k0 = numpy.loadtxt(self.fname1) * 1e-10 * u.cm ** 3 / u.s
            delta_E0 = numpy.loadtxt(self.fname2) * u.meV
            _t_vals = numpy.logspace(1., 2.7, 22)

            K_gerlich = []
            for itemp, t_kin in enumerate(_t_vals):
                t_kin *= u.K
                K_gerlich_array = k0 * numpy.exp(
                    -1. * delta_E0.to(u.J) / k_B / t_kin
                )
                K_gerlich.append(K_gerlich_array)

            n_transitions = len(K_gerlich[:][0])**2

            ini = numpy.vstack(
                (
                    numpy.zeros(n_transitions),  # all the v's are zero
                    numpy.repeat([0, 1, 2, 3, 4, 5, 6, 7], len(K_gerlich[:][0]))
                )
            ).astype('i')

            fin = numpy.vstack(
                (
                    numpy.zeros(n_transitions),  # all the v's are zero
                    numpy.tile([0, 1, 2, 3, 4, 5, 6, 7], len(K_gerlich[:][0]))
                )
            ).astype('i')

            # DEBUG
            # for i, f in zip(ini[1], fin[1]):
            #     print(i, f, delta_E0[i, f], k0[i, f])

            # declare the array where the data will be stored in a tensor
            nv, nj, nvp, njp = \
                0 + 1, len(K_gerlich[:][0]), 0 + 1, len(K_gerlich[:][0])

            nv_max, nj_max = max(nv, nvp), max(nj, njp)
            data = zeros((_t_vals.size, nv_max, nj_max, nv_max, nj_max), 'f8')
            # cr = zeros((_t_vals.size, len(_v)), 'f8')
            cr = zeros((_t_vals.size, n_transitions), 'f8')

            for itemp, _ in enumerate(_t_vals):
                for j in numpy.arange(len(K_gerlich[:][0])):
                    for jp in numpy.arange(len(K_gerlich[:][0])):
                        data[itemp, 0, j, 0, jp] = K_gerlich[itemp][j][jp].value
                        # assign the reaction rates for all the temperatures
                        # for the ib-th transition
                        cr[itemp, j*len(K_gerlich[:][0])+jp] =\
                            K_gerlich[itemp][j][jp].value

            self.ini = ini
            self.fin = fin
            self.tkin = _t_vals
            self.data = data
            self.cr = cr

            return data, cr

    reader = Reader(fname1, fname2)

    ini = reader.ini
    fin = reader.fin

    # find the unique levels from from the transitions
    unique_levels = unique_level_pairs(
        hstack(
            (unique_level_pairs(ini),
             unique_level_pairs(fin))
        )
    )

    # set the units of the data to be returned
    data_with_units = reader.data * (u.cm**3 / u.second)
    cr_with_units = reader.cr * (u.cm**3 / u.second)
    t_values = reader.tkin * u.K

    # convert the units to m^3/s
    data_with_units = data_with_units.to(u.m**3 / u.second)
    cr_with_units = cr_with_units.to(u.m**3 / u.second)

    return data_with_units, t_values, (ini, fin, unique_levels, cr_with_units)
