# -*- coding: utf-8 -*-

#    read_einstein_coefficient.py is part of Frigus.

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
read the Einstein coefficients
"""
import os
import numpy
from numpy import zeros, testing, array

from astropy import units as u

from frigus import utils

DATADIR = utils.datadir_path()


def read_einstein_simbotin():
    """
    Read the data provided by Simbotin from multiple files and returns the A
    matrix for transitions (v', j') -> (v'', j'') with the limitation that
    delta j i.e |j'' - j'| = 0 or 2

    0   j''=j'
                  14              13              12              11              10               9         .....     1
    1    0   0.6094641D-13   0.1367911D-12   0.1860857D-12   0.1826816D-12   0.9029585D-13   0.8476571D-14   .....

    The data is read from the following files:
          Read/j2jdown
          Read/j2j
          Read/j2jup

    :return: A 4D matrix holding the A coefficients. A[v', j', v'', j'']

    .. code-block:: python

        # get the data as a 4D matrix
        A = read_einstein_simbotin()

        # get the Einstein coefficient for the transition (v'=3, j'=9) to
        # (v''=0, j''=18)
        print(A[3, 9, 0, 18])
    """

    # lists that store the read data. These are the nonzero entries of A
    # that is also returned
    vp_nnz, jp_nnz, vpp_nnz, jpp_nnz, A_nnz = [], [], [], [], []

    def read_j2j_x(fname, A, delta_j=None, skip_rows=None):
        """
        Read the Einstein coefficients from the file "fname" and modify the
        corresponding entries of 'A'.

        .. warning:: this function modifies the values of A

        in the snippet below:
        - the first row is initial vibrational level v'
        - the first column is the initial rotational level j'
        - the second column is the final vibrational level v''
        - the final rational level can be obtained from j' (the first column)
          by adding 'delta_j'

           j'   v''| v'      14              13              12              11              10               9               8               7               6               5               4               3               2               1               0
           2    0  |     0.2237387D-13   0.7602406D-13   0.1977334D-12   0.5213964D-12   0.1451901D-11   0.4314544D-11   0.1378568D-10   0.4779910D-10   0.1820544D-09   0.7742611D-09   0.3770442D-08   0.2158490D-07   0.1274956D-06   0.2526477D-06   0.2941861D-10
           2    1  |     0.9136254D-12   0.2921039D-11   0.6942250D-11   0.1660788D-10   0.4217889D-10   0.1154684D-09   0.3432847D-09   0.1116306D-08   0.4011464D-08   0.1621802D-07   0.7435754D-07   0.3193196D-06   0.3682392D-06   0.2785271D-10

        :param fname: path to the file containing the Einstein coefficient data
        :param array_like A: The 4D numpy array that will be populated with the read data
        :param int delta_j: the increment that is added to  j' to obtain j'' of
         a transition.
        :param int skip_rows: The number of rows to skip in parsing the data
         file. (i.e the number of lines of the header).
        """

        # opening the original ascii file and discard empty lines
        with open(fname) as fobj:
            lines = list(filter(lambda x: x.strip() != '', fobj.readlines()))

        lines = lines[skip_rows:]
        while True:

            if len(lines) == 0:
                break

            line = lines.pop(0)

            vp_all = list(map(int, line.split()))

            while True:
                # processing a block and break the loop only when the last line
                # of the block is identified as having 3 tokens
                line = lines.pop(0)
                tokens = line.split()
                jp, vpp = int(tokens[0]), int(tokens[1])
                for i, A_i in enumerate(tokens[2:]):
                    A_tmp = numpy.float64(A_i.replace('D', 'E'))
                    vp, jpp = vp_all[i], jp + delta_j
                    A[vp, jp, vpp, jpp] = A_tmp

                    if vp == vpp and jp == jpp:
                        raise ValueError('v -> v, j -> j transition. This is '
                                         'not possible')
                    vp_nnz.append(vp)
                    jp_nnz.append(jp)
                    vpp_nnz.append(vpp)
                    jpp_nnz.append(jpp)
                    A_nnz.append(A_tmp)

                # We determine that the end of a block has reached by checking
                # that the first 10 chars of next line are empty
                if len(lines) >= 1:
                    next_line = lines[0]
                    if next_line[0:10].strip() == '':
                        break
                else:
                    break

    # get the vmax and the jmax from the file j2jdown (the info is found only
    # in that file and not in j2j nor in j2jup)
    # The greatest rotational number for H2, starting from j = 0
    # (tot. rotational numbers: jmax=32). The greatest vibrational number for
    # H2, starting from v = 0 (tot. vibrational numbers: vmax=15)
    with open(os.path.join(DATADIR, 'j2jdown')) as fobj:
        lines = fobj.readlines()
    jmax = numpy.array(lines[1].split(), dtype='i').max()
    vmax = numpy.array(lines[2].split(), dtype='i').max()

    # define the A matrix (vmax and jmax assume zero indexing that is how they
    # they are provided in the data files)
    A = zeros((vmax + 1, jmax + 1, vmax + 1, jmax + 1), 'f8')

    read_j2j_x(os.path.join(DATADIR, 'j2jdown'), A, delta_j=-2, skip_rows=3)
    read_j2j_x(os.path.join(DATADIR, 'j2j'),     A, delta_j=0, skip_rows=1)
    read_j2j_x(os.path.join(DATADIR, 'j2jup'),   A, delta_j=2, skip_rows=1)

    testing.assert_approx_equal(A.sum(), 9.3724e-4, significant=4)
    testing.assert_approx_equal(A[7, 2, 4, 2], 4.133e-7, significant=4)

    # setting units to the Einstein coefficients that will be returned
    A_with_units = A / u.second
    A_nnz_with_units = A_nnz / u.second

    return A_with_units, (
        array(vp_nnz, 'i4'),
        array(jp_nnz, 'i4'),
        array(vpp_nnz, 'i4'),
        array(jpp_nnz, 'i4'),
        A_nnz_with_units
    )


def read_einstein_coppola():
    """
    Read the data provided by Coppola for transitions (v', j') -> (v'', j'')
    with the limitation that delta j i.e |j'' - j'| = 1

    :return: A 4D matrix holding the A coefficients. A[v', j', v'', j'']

    .. code-block:: python

        .. todo:: add example
    """

    # lists that store the read data. These are the nonzero entries of A
    # that is also returned
    vp_nnz, jp_nnz, vpp_nnz, jpp_nnz, A_nnz = [], [], [], [], []

    data = numpy.loadtxt(
        os.path.join(DATADIR, 'lipovka', 'hd_einstein_coeffs.dat'),
        skiprows=2
    )

    vp_nnz = data[:, 0].astype(numpy.int32)
    jp_nnz = data[:, 1].astype(numpy.int32)
    vpp_nnz = data[:, 2].astype(numpy.int32)
    jpp_nnz = data[:, 3].astype(numpy.int32)
    A_nnz = data[:, 5].astype(numpy.float64)

    jmax = jp_nnz.max()
    vmax = vp_nnz.max()

    # define the A matrix (vmax and jmax assume zero indexing that is how they
    # they are provided in the data files)
    A = zeros((vmax + 1, jmax + 1, vmax + 1, jmax + 1), 'f8')

    A[vp_nnz, jp_nnz, vpp_nnz, jpp_nnz] = A_nnz

    # .. todo:: add these tests
    # testing.assert_approx_equal(A.sum(), 9.3724e-4, significant=4)
    # testing.assert_approx_equal(A[7, 2, 4, 2], 4.133e-7, significant=4)

    # setting units to the Einstein coefficients that will be returned
    A_with_units = A / u.second
    A_nnz_with_units = A_nnz / u.second

    return A_with_units, (
        vp_nnz,
        jp_nnz,
        vpp_nnz,
        jpp_nnz,
        A_nnz_with_units)
